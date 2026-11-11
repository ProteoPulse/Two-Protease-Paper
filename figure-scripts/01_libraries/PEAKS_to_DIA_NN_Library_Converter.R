# --------------------------------------------------------------------------------------------------
# Two-Protease-Paper | Convert PEAKS 13.1 library exports to DIA-NN library format
# Manuscript: Methods, 'Data Analysis' (PEAKS libraries -> DIA-NN).
# Author: Simion Kreimer.
# Run: build_unimod_tsv() once (needs unimod.xml from unimod.org), then
#      convert_peaks_to_diann('peaks.csv', 'proteome.fasta', 'unimod_full.tsv').
# Packages: data.table, stringr, tools, xml2 (installed on demand).
# Released as used for the manuscript: only file paths, line endings and this header were changed;
# the analysis code is untouched. Paths come from environment variables, see ../../README.md.
# --------------------------------------------------------------------------------------------------

# =========================
# PEAKS -> DIA-NN Converter
# with file browser, full UniMod TSV import,
# FASTA protein mapping, robust column resolver,
# PEAKS-rounded PTM matching, fragment loss parsing,
# PTM summary + unmatched/ambiguous logging
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(tools)
})

# ---------------------------
# File chooser helpers
# ---------------------------
choose_file_one <- function(prompt = "Select a file") {
  cat(prompt, "...\n")
  if (.Platform$OS.type == "windows") {
    p <- choose.files()
    if (length(p) == 0 || p == "") stop("No file selected.")
    return(p[1])
  } else {
    return(file.choose())
  }
}

choose_save_path <- function(default_name = "unimod_full.tsv", prompt = "Choose where to save") {
  cat(prompt, "...\n")
  if (.Platform$OS.type == "windows") {
    p <- choose.files(default = default_name)
    if (length(p) == 0 || p == "") stop("No save path selected.")
    return(p[1])
  } else {
    return(file.path(getwd(), default_name))
  }
}

# ---------------------------
# Column resolver (robust to header variants)
# ---------------------------
resolve_col <- function(dt, candidates, required = TRUE) {
  nms <- names(dt)
  nms_l <- tolower(nms)
  cand_l <- tolower(candidates)
  
  # exact (case-insensitive)
  hit <- nms[match(cand_l, nms_l)]
  hit <- hit[!is.na(hit)]
  if (length(hit) > 0) return(hit[1])
  
  # fuzzy contains
  for (cand in cand_l) {
    idx <- which(grepl(cand, nms_l, fixed = TRUE))
    if (length(idx) > 0) return(nms[idx[1]])
  }
  
  if (required) {
    stop(
      "Could not find a required column. Tried: ", paste(candidates, collapse = ", "),
      "\nAvailable columns:\n", paste(nms, collapse = ", ")
    )
  }
  NA_character_
}

# ---------------------------
# Build FULL UniMod TSV locally from downloaded UniMod XML
# (run once; then the converter will import the TSV via file browser)
# ---------------------------
build_unimod_tsv <- function(unimod_xml_path = NULL, out_tsv_path = NULL) {
  if (!requireNamespace("xml2", quietly = TRUE)) install.packages("xml2")
  suppressPackageStartupMessages({
    library(xml2)
    library(data.table)
  })
  
  if (is.null(unimod_xml_path)) unimod_xml_path <- choose_file_one("Select UniMod XML (unimod.xml)")
  if (is.null(out_tsv_path)) out_tsv_path <- choose_save_path("unimod_full.tsv", "Choose where to save the UniMod TSV")
  
  doc <- xml2::read_xml(unimod_xml_path)
  
  mod_nodes <- xml2::xml_find_all(doc, ".//*[local-name()='mod']")
  cat("Found <mod> nodes:", length(mod_nodes), "\n")
  if (length(mod_nodes) == 0) stop("No <mod> entries found. Wrong file?")
  
  unimod_list <- lapply(mod_nodes, function(m) {
    id <- suppressWarnings(as.integer(xml2::xml_attr(m, "record_id")))
    
    # UniMod unimod.xml typically stores name in attributes on <mod>
    nm <- xml2::xml_attr(m, "title")
    if (is.na(nm) || nm == "") nm <- xml2::xml_attr(m, "full_name")
    if (is.na(nm) || nm == "") nm <- xml2::xml_text(m)
    
    # mono mass stored as attribute on <delta mono_mass="...">
    delta_node <- xml2::xml_find_first(m, ".//*[local-name()='delta']")
    mass_attr <- xml2::xml_attr(delta_node, "mono_mass")
    if (is.na(mass_attr) || mass_attr == "") mass_attr <- xml2::xml_attr(delta_node, "delta_mono_mass")
    mass <- suppressWarnings(as.numeric(mass_attr))
    
    data.table(Name = nm, Id = id, MonoMass = mass)
  })
  
  unimod_df <- rbindlist(unimod_list, fill = TRUE)
  unimod_df <- unimod_df[!is.na(Id) & !is.na(MonoMass) & !is.na(Name) & nzchar(Name)]
  
  if (nrow(unimod_df) == 0) {
    stop("Parsed 0 UniMod rows. The XML schema is unexpected or name/mass fields weren't found.")
  }
  
  setkey(unimod_df, Id)
  fwrite(unimod_df, out_tsv_path, sep = "\t", quote = FALSE, na = "NA")
  
  cat("✔ Saved full UniMod table to:\n", out_tsv_path, "\n")
  cat("Rows written:", nrow(unimod_df), "\n")
  cat("First 5 rows:\n")
  print(head(unimod_df, 5))
  
  invisible(out_tsv_path)
}

# ---------------------------
# FASTA parsing and peptide->protein mapping
# ---------------------------
read_fasta_to_table <- function(fasta_file) {
  lines <- readLines(fasta_file, warn = FALSE)
  headers <- grep("^>", lines)
  if (length(headers) == 0) stop("No FASTA headers found (lines starting with '>').")
  
  starts <- headers + 1
  stops <- c(headers[-1] - 1, length(lines))
  
  rbindlist(lapply(seq_along(headers), function(i) {
    seq <- paste(lines[starts[i]:stops[i]], collapse = "")
    header <- sub("^>", "", lines[headers[i]])
    
    # Extract accession between first two pipes if present: sp|O00178|GTPB1_HUMAN -> O00178
    acc <- if (grepl("\\|", header)) sub(".*\\|(.*?)\\|.*", "\\1", header) else sub(" .*", "", header)
    
    gene <- if (grepl("GN=", header)) sub(".*GN=([^ ]+).*", "\\1", header) else NA_character_
    
    data.table(Protein = acc, Gene = gene, Sequence = seq)
  }))
}

match_peptides_to_proteins <- function(peptide_seqs, fasta_df) {
  rbindlist(lapply(peptide_seqs, function(peptide) {
    hits <- fasta_df[grepl(peptide, Sequence, fixed = TRUE)]
    data.table(
      PeptideSequence = peptide,
      ProteinId = if (nrow(hits)) hits$Protein[1] else NA_character_,
      AllMappedProteins = if (nrow(hits)) paste(unique(hits$Protein), collapse = ";") else NA_character_,
      GeneName = if (nrow(hits)) hits$Gene[1] else NA_character_,
      AllMappedGenes = if (nrow(hits)) paste(unique(na.omit(hits$Gene)), collapse = ";") else NA_character_,
      Proteotypic = (nrow(hits) == 1)
    )
  }))
}

# ---------------------------
# Mod parsing: PEAKS Modifications -> [UniMod:ID] insertion
# Supports negative deltas and PEAKS rounding to 2 decimals
# ---------------------------
unmatched_mods <- character()
ambiguous_mods <- character()
cam_offsite_examples <- character()
# Helper: evaluate mapping quality using UniMod:4->C and UniMod:35->M
score_unimod_localization <- function(modseq) {
  # capture residue immediately before [UniMod:ID]
  m <- stringr::str_match_all(modseq, "([A-Z])\\[UniMod:(\\d+)\\]")[[1]]
  if (nrow(m) == 0) return(list(n4 = 0, ok4 = 0, n35 = 0, ok35 = 0, score = 0))
  
  aa  <- m[, 2]
  uid <- as.integer(m[, 3])
  
  n4  <- sum(uid == 4,  na.rm = TRUE)
  ok4 <- sum(uid == 4  & aa == "C", na.rm = TRUE)
  
  n35  <- sum(uid == 35, na.rm = TRUE)
  ok35 <- sum(uid == 35 & aa == "M", na.rm = TRUE)
  
  # score: fraction correct among those two mods (avoid div-by-zero)
  denom <- (n4 + n35)
  score <- if (denom > 0) (ok4 + ok35) / denom else 0
  
  list(n4 = n4, ok4 = ok4, n35 = n35, ok35 = ok35, score = score)
}

map_to_unimod <- function(seq, mods_str, unimod_df,
                          index_base = "0",     # PEAKS example shows 0-based
                          mass_tol = 0.03) {
  
  if (is.na(seq) || seq == "" || is.na(mods_str) || mods_str == "") return(seq)
  
  # normalize unicode minus
  mods_str <- gsub("\u2212|\u2013|\u2014", "-", mods_str)
  
  seq_split <- strsplit(seq, "")[[1]]
  prefix <- ""
  suffix <- ""
  
  # split PEAKS tokens (your example uses ;)
  tokens <- unlist(strsplit(mods_str, ";", fixed = TRUE))
  tokens <- trimws(tokens)
  tokens <- tokens[tokens != ""]
  
  # helper: parse one token into (pos, name, context, delta)
  parse_token <- function(tok) {
    # Style A (your example): 3-Carbamidomethylation-(57.02)
    m1 <- stringr::str_match(tok, "^(\\d+)-(.+?)-\\(([+-]?[0-9.]+)\\)$")
    if (!any(is.na(m1))) {
      return(list(
        pos = as.integer(m1[2]),
        name = trimws(m1[3]),
        context = "",
        delta = as.numeric(m1[4])
      ))
    }
    
    # Style B: 3-Name (context)-(57.02)  OR  3-Name(context)-(57.02)
    m2 <- stringr::str_match(tok, "^(\\d+)-(.+?)\\s*\\(([^)]*)\\)-\\(([+-]?[0-9.]+)\\)$")
    if (!any(is.na(m2))) {
      return(list(
        pos = as.integer(m2[2]),
        name = trimws(m2[3]),
        context = m2[4],
        delta = as.numeric(m2[5])
      ))
    }
    
    # If neither matches, give up on this token
    return(NULL)
  }
  
  # choose unimod id by (1) name rules, then (2) mass tolerance
  choose_unimod_id <- function(mod_name_l, delta_raw) {
    # Force carbamid* -> UniMod:4
    if (grepl("carbamid", mod_name_l, fixed = TRUE)) return(4L)
    
    # amidation: common +/− ambiguity; UniMod amidation is negative, but we mainly need correct ID
    if (grepl("amidation", mod_name_l, fixed = TRUE)) {
      # If your unimod table contains amidation, the name match will usually find it; otherwise mass match.
      # Keep delta as-is; just let name match win if present.
    }
    
    # mass tolerance candidates
    cand <- unimod_df[abs(MonoMass - delta_raw) <= mass_tol]
    if (nrow(cand) == 0) return(NA_integer_)
    
    # prefer name hit among candidates
    nh <- cand[grepl(mod_name_l, tolower(Name), fixed = TRUE)]
    if (nrow(nh) >= 1) return(nh$Id[1])
    return(cand$Id[1])
  }
  
  for (tok in tokens) {
    info <- parse_token(tok)
    if (is.null(info)) {
      unmatched_mods <<- unique(c(unmatched_mods, paste0("Unparsed token: ", tok)))
      next
    }
    
    pos_raw <- info$pos
    mod_name_raw <- info$name
    context <- info$context
    delta_raw <- info$delta
    if (is.na(pos_raw) || is.na(delta_raw)) next
    
    mod_name_l <- tolower(mod_name_raw)
    
    # Terminal only if explicitly stated (some PEAKS exports include this)
    is_nterm <- grepl("N-term", context, ignore.case = TRUE)
    is_cterm <- grepl("C-term", context, ignore.case = TRUE)
    
    # Force amidation to C-term behavior if it appears without context
    if (grepl("amidation", mod_name_l, fixed = TRUE)) {
      is_cterm <- TRUE
      is_nterm <- FALSE
    }
    
    unimod_id <- choose_unimod_id(mod_name_l, delta_raw)
    
    if (is.na(unimod_id)) {
      unmatched_mods <<- unique(c(unmatched_mods, paste0(mod_name_raw, " (", round(delta_raw, 2), ")")))
      next
    }
    
    tag <- paste0("[UniMod:", unimod_id, "]")
    
    # CRITICAL FIX: UniMod:4 must NOT be treated as N-term/C-term
    if (unimod_id == 4L) { is_nterm <- FALSE; is_cterm <- FALSE }
    
    if (is_nterm) { prefix <- paste0(prefix, tag); next }
    if (is_cterm) { suffix <- paste0(suffix, tag); next }
    
    # Residue placement (after residue)
    idx <- if (index_base == "0") pos_raw + 1L else pos_raw
    
    if (idx >= 1 && idx <= length(seq_split)) {
      # Optional: log offsite CAM placements (you said some exist)
      if (unimod_id == 4L && seq_split[idx] != "C") {
        cam_offsite_examples <<- unique(c(cam_offsite_examples,
                                          paste0(seq, " | pos=", pos_raw, " -> ", seq_split[idx], " | token=", tok)))
      }
      
      seq_split[idx] <- paste0(seq_split[idx], tag)
    } else {
      unmatched_mods <<- unique(c(unmatched_mods, paste0(mod_name_raw, " (", round(delta_raw, 2), ") [pos=", pos_raw, " out of range]")))
    }
  }
  
  paste0(prefix, paste(seq_split, collapse = ""), suffix)
}

# ---------------------------
# Fragment loss parsing: b2-H2O, y7-NH3, etc.
# ---------------------------
parse_fragment_annotation <- function(ann) {
  frag_type <- str_extract(ann, "^[abcxyz]")
  series <- suppressWarnings(as.integer(str_extract(ann, "\\d+")))
  
  loss <- "noloss"
  if (grepl("-", ann, fixed = TRUE)) {
    loss_raw <- sub("^[^-]+-", "", ann)
    loss <- if (!is.na(loss_raw) && nzchar(loss_raw)) loss_raw else "loss"
  }
  
  list(FragmentType = frag_type, FragmentSeriesNumber = series, FragmentLossType = loss)
}

# ---------------------------
# MAIN: Convert PEAKS -> DIA-NN
# Uses file dialogs by default
# ---------------------------
convert_peaks_to_diann <- function(peaks_filepath = NULL,
                                   fasta_filepath = NULL,
                                   unimod_tsv_filepath = NULL) {
  
  # pick files if not provided
  if (is.null(peaks_filepath)) peaks_filepath <- choose_file_one("Select PEAKS spectral library file")
  if (is.null(fasta_filepath)) fasta_filepath <- choose_file_one("Select FASTA protein database")
  if (is.null(unimod_tsv_filepath)) unimod_tsv_filepath <- choose_file_one("Select FULL UniMod TSV (build once with build_unimod_tsv())")
  
  # output file
  file_base <- file_path_sans_ext(basename(peaks_filepath))
  output_file <- file.path(dirname(peaks_filepath), paste0(file_base, "_DIANN.tsv"))
  
  # load inputs
  peaks <- fread(peaks_filepath)
  
  # Resolve PEAKS column variants
  col_prec_mz   <- resolve_col(peaks, c("m/z", "mz", "precursor m/z", "precursor mz", "precursormz"))
  col_z         <- resolve_col(peaks, c("z", "charge", "charge state", "precursor charge", "precursorcharge"))
  col_rt        <- resolve_col(peaks, c("rt", "retention time", "retention time (min)", "rt (min)", "rt[min]", "rt [min]", "rt(min)", "rt_sec", "rt (sec)", "rt [sec]"))
  col_seq       <- resolve_col(peaks, c("sequence (backbone)", "sequence", "peptide sequence", "peptidesequence", "peptide"))
  col_mods      <- resolve_col(peaks, c("modifications", "modification", "mods"), required = FALSE)
  col_peakslist <- resolve_col(peaks, c("peaks list", "peakslist", "fragment ions", "fragments", "peakslist_in", "peaks list "), required = TRUE)
  
  # Standardize to internal names
  setnames(peaks,
           old = c(col_prec_mz, col_z, col_rt, col_seq, col_peakslist),
           new = c("PrecursorMz_in", "PrecursorCharge_in", "RT_in", "Seq_in", "PeaksList_in"))
  
  # Mods column optional: if not found, create blank
  if (is.na(col_mods)) {
    peaks[, Mods_in := ""]
  } else {
    setnames(peaks, old = col_mods, new = "Mods_in")
  }
  
  # Load UniMod TSV
  unimod_df <- fread(unimod_tsv_filepath)
  if (!all(c("Name", "Id", "MonoMass") %in% names(unimod_df))) {
    stop("UniMod TSV must have columns: Name, Id, MonoMass")
  }
  unimod_df[, `:=`(Id = as.integer(Id), MonoMass = as.numeric(MonoMass))]
  unimod_df <- unimod_df[!is.na(Id) & !is.na(MonoMass) & !is.na(Name) & nzchar(Name)]
  setkey(unimod_df, Id)
  
  # FASTA mapping
  fasta_df <- read_fasta_to_table(fasta_filepath)
  pep_map <- match_peptides_to_proteins(unique(peaks$Seq_in), fasta_df)
  peaks <- merge(peaks, pep_map, by.x = "Seq_in", by.y = "PeptideSequence", all.x = TRUE)
  
  # ModifiedPeptideSequence
  peaks[, ModifiedPeptideSequence := mapply(
    map_to_unimod,
    Seq_in,
    Mods_in,
    MoreArgs = list(unimod_df = unimod_df, index_base = "0", mass_tol = 0.03)
  )]
  
  
  chk <- rbindlist(lapply(peaks$ModifiedPeptideSequence, function(s) as.data.table(score_unimod_localization(s))), fill = TRUE)
  cat(sprintf("Unimod:4 on C: %d/%d | Unimod:35 on M: %d/%d\n",
              sum(chk$ok4), sum(chk$n4), sum(chk$ok35), sum(chk$n35)))
  
  # PTM usage summary
  ptm_ids <- str_extract_all(peaks$ModifiedPeptideSequence, "UniMod:\\d+")
  ptm_ids_flat <- unlist(ptm_ids)
  ptm_counts <- sort(table(ptm_ids_flat), decreasing = TRUE)
  
  # Expand fragments into DIA-NN format
  entries <- vector("list", 0L)
  k <- 0L
  
  for (i in seq_len(nrow(peaks))) {
    e <- peaks[i]
    pl <- e$PeaksList_in
    if (is.na(pl) || pl == "") next
    
    frags <- strsplit(pl, ";", fixed = TRUE)[[1]]
    for (frag in frags) {
      parts <- strsplit(frag, ":", fixed = TRUE)[[1]]
      if (length(parts) < 3) next
      
      prod_mz <- suppressWarnings(as.numeric(parts[1]))
      lib_int <- suppressWarnings(as.numeric(parts[2]))
      ann <- parts[3]
      
      fp <- parse_fragment_annotation(ann)
      
      k <- k + 1L
      entries[[k]] <- data.table(
        PrecursorMz = e$PrecursorMz_in,
        ProductMz = prod_mz,
        Annotation = ann,
        ProteinId = e$ProteinId,
        GeneName = e$GeneName,
        PeptideSequence = e$Seq_in,
        ModifiedPeptideSequence = e$ModifiedPeptideSequence,
        PrecursorCharge = e$PrecursorCharge_in,
        LibraryIntensity = lib_int,
        NormalizedRetentionTime = e$RT_in,
        PrecursorIonMobility = NA,
        FragmentType = fp$FragmentType,
        FragmentCharge = 1,
        FragmentSeriesNumber = fp$FragmentSeriesNumber,
        FragmentLossType = fp$FragmentLossType,
        AverageExperimentalRetentionTime = e$RT_in,
        AllMappedProteins = e$AllMappedProteins,
        AllMappedGenes = e$AllMappedGenes,
        Proteotypic = e$Proteotypic
      )
    }
  }
  
  out <- rbindlist(entries, fill = TRUE)
  fwrite(out, output_file, sep = "\t", quote = FALSE, na = "NA")
  cat("✔ DIA-NN spectral library written to:\n", output_file, "\n")
  
  # Report PTMs used
  if (length(ptm_counts) > 0) {
    cat("\n📌 PTMs used (from ModifiedPeptideSequence):\n")
    ptm_dt <- data.table(UniModTag = names(ptm_counts), Count = as.integer(ptm_counts))
    ptm_dt[, Id := as.integer(sub("UniMod:", "", UniModTag))]
    ptm_dt <- merge(ptm_dt, unimod_df[, .(Id, Name, MonoMass)], by = "Id", all.x = TRUE)
    setorder(ptm_dt, -Count)
    
    # Print top 50
    to_print <- head(ptm_dt, 50)
    for (r in seq_len(nrow(to_print))) {
      cat(sprintf(" - %s [UniMod:%d] (Δ=%.6f): %d\n",
                  to_print$Name[r], to_print$Id[r], to_print$MonoMass[r], to_print$Count[r]))
    }
    if (nrow(ptm_dt) > 50) cat(" ... (showing top 50)\n")
  } else {
    cat("\n📌 PTMs used: none detected in ModifiedPeptideSequence.\n")
  }
  
  # Unmatched + ambiguous logs
  if (length(unmatched_mods) > 0) {
    cat("\n⚠️ Unmatched modifications (PEAKS -> UniMod):\n")
    for (m in unmatched_mods) cat(" - ", m, "\n", sep = "")
  } else {
    cat("\n✅ No unmatched modifications.\n")
  }
  
  if (length(ambiguous_mods) > 0) {
    cat("\n⚠️ Ambiguous modifications (mass matched multiple UniMod entries):\n")
    for (m in head(ambiguous_mods, 50)) cat(" - ", m, "\n", sep = "")
    if (length(ambiguous_mods) > 50) cat(" ... (showing first 50)\n")
  }
  if (length(cam_offsite_examples) > 0) {
    cat("\n[CAM offsite examples] showing first 20:\n")
    print(head(cam_offsite_examples, 20))
  }
  
  invisible(output_file)
}

# ---------------------------
# HOW TO RUN
# ---------------------------
# 1) Download unimod.xml from UniMod and build the TSV once:
#    build_unimod_tsv()
#
# 2) Convert with dialogs:
#    convert_peaks_to_diann()
#
# Optional: pass paths directly:
# convert_peaks_to_diann("peaks.csv", "proteome.fasta", "unimod_full.tsv")
