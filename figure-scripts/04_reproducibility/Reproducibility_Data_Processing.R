# --------------------------------------------------------------------------------------------------
# Two-Protease-Paper | Replicate CV analysis of the reproducibility cohort
# Manuscript: Fig 4A (inter-day CV violins, well-quantified counts per sample, measured data).
# Authors: Simion Kreimer (original); Mario M. Alba (manuscript configuration).
# Inputs: TWO_PROTEASE_DATA/Reproducibility_Search/{Trypsin,Krakatoa} Library/report.{pr,pg}_matrix.tsv;
#         TWO_PROTEASE_FASTA for PTM residue numbering.
# Outputs: TWO_PROTEASE_OUT/04_reproducibility/DIANN_Replicate_CV_QC_<timestamp>/
# Run: Rscript figure-scripts/04_reproducibility/Reproducibility_Data_Processing.R (prompts take defaults).
# Packages: data.table, ggplot2, grid.
# Released as used for the manuscript: only file paths, line endings and this header were changed;
# the analysis code is untouched. Paths come from environment variables, see ../../README.md.
# --------------------------------------------------------------------------------------------------

# ============================================================
# DIA-NN replicate CV QC after combined search results
# Purpose:
#   1) Read one or more DIA-NN search result folders.
#   2) Build protein and peptide matrices for each search strategy.
#   3) Build an "All Searches" combined matrix from all strategies.
#   4) Parse sample names as Day_SampleID_Replicate from sample columns.
#   5) Calculate feature-level intra-day and inter-day CVs.
#   6) Export low-CV features detected in all expected replicates.
#   7) Export well-quantified PTM peptides with protein/gene context.
#   8) Plot CV distributions as violin plots, separately for each search strategy
#      and for "All Searches".
#   9) Summarize search-strategy overlap of well-quantified proteins.
#
# Expected sample name:
#   1_2_3
#   where:
#     1 = preparation day
#     2 = sample ID
#     3 = replicate analysis number
#
# If your raw column names are full paths, the script uses the file basename and
# removes .raw/.d/.mzML before parsing.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

# ==============================
# USER SETTINGS
# ==============================
pg_file <- "report.pg_matrix.tsv"
pr_file <- "report.pr_matrix.tsv"

ALL_SEARCHES_LABEL <- "All Searches"

# Peptide intensity collapsing from pr_matrix rows to modified peptide rows.
# "sum" is usually best if multiple precursors map to the same modified peptide.
peptide_collapse_method <- "sum"  # "sum" or "mean"

zero_as_na <- TRUE

# If TRUE, columns that do not parse as Day_SampleID_Replicate are written to
# sample_name_parse_failures.csv and excluded from CV analysis instead of stopping.
drop_unparsed_sample_columns <- TRUE

# Parser behavior:
#   "last_numeric_triplet"  = use the last 3 underscore/dash/dot separated numeric tokens
#   "first_numeric_triplet" = use the first 3 numeric tokens
#   "strict_first3"         = require first 3 tokens to be Day/SampleID/Replicate
sample_parse_mode <- "last_numeric_triplet" 

# CV thresholds used to define "well quantified".
intra_cv_threshold <- 20
inter_cv_threshold <- 20

# PTMs/modified peptides can be noisier; use a separate inter-day CV threshold.
ptm_inter_cv_threshold <- 40

# PTM plot: only show modified sites that pass in all samples.
# Set to NULL to use all unique SampleID values detected in the experiment.
ptm_required_samples_for_plot <- 8

# Required for exact full-protein PTM residue numbering unless DIA-NN already exports
# a protein-position/site column. Set manually or choose interactively.
# Example: fasta_file <- "C:/path/to/proteome.fasta"
# TWO_PROTEASE_FASTA = UniProtKB human reviewed proteome (see the repository README).
fasta_file <- Sys.getenv("TWO_PROTEASE_FASTA", "uniprotkb_Human_AND_reviewed_true_AND_m_2025_10_28.fasta")
fasta_db <- list(seq_by_accession = list(), seq_by_header = list())

# Minimum number of replicate injections required to calculate intra-day CV.
min_replicates_for_intra_cv <- 2

# Minimum number of preparation days required to calculate inter-day CV.
min_days_for_inter_cv <- 2

# If TRUE, "well quantified intra-day" requires detection in every replicate
# that exists for that SampleID + Day.
require_all_replicates_for_well_quantified <- TRUE

# If TRUE, "well quantified across days" requires a feature to pass intra-day QC
# on every day available for that SampleID, and pass inter-day CV.
require_all_days_for_consensus <- TRUE

# Save high-resolution PNG and PDF plots.
PUB_DPI <- 600

# Maximum number of sample/day panels shown in a single violin plot before
# switching x-axis text smaller.
large_plot_label_threshold <- 16

# Well-quantified protein replicate intensity plots.
max_wq_proteins_per_intensity_plot <- Inf

# Sample fingerprint plots:
# Each sample gets a panel, each well-quantified protein is shown as a box/whisker
# across replicate measurements, ordered by global average intensity.
max_wq_proteins_per_fingerprint_plot <- Inf
show_protein_names_on_fingerprint_xaxis <- FALSE
wq_intensity_log_base <- 2

# ==============================
# FOLDER PICKER
# ==============================
pick_parent_folder <- function() {
  if (interactive() && .Platform$OS.type == "windows") {
    p <- utils::choose.dir(caption = "Select parent folder containing DIA-NN search result folders")
    if (is.na(p) || !nzchar(p)) stop("No folder selected.")
    return(p)
  }
  stop("Folder picker is Windows interactive only. Set base_dir manually near the top of the script.")
}

# Comment this line and set base_dir manually if not running interactively on Windows.
# base_dir <- pick_parent_folder()                      # interactive picker (disabled)
# 2-protease paper: replicate inputs live in Results/Reproducibility_Search (its immediate
# children = "Krakatoa Library" / "Trypsin Library"); discovery is non-recursive, so base_dir
# points there. Output is decoupled below to land directly in the Results root as a timestamped folder.
# Paths come from environment variables (see the repository README): TWO_PROTEASE_DATA = the DIA-NN
# search results (MassIVE MSV000103457), TWO_PROTEASE_OUT = where reports are written.
base_dir       <- file.path(Sys.getenv("TWO_PROTEASE_DATA", "data"), "Reproducibility_Search")
results_root   <- file.path(Sys.getenv("TWO_PROTEASE_OUT", "output"), "04_reproducibility")

# ---- Paste/Rscript-safe prompts (2-protease): take all default answers --------
# Every readline() returns "" so each prompt uses its DEFAULT branch and can never
# consume a following pasted line. Set ASSUME_DEFAULTS <- FALSE to answer prompts
# yourself when running via source() in an interactive console.
ASSUME_DEFAULTS <- TRUE
if (ASSUME_DEFAULTS) readline <- function(prompt = "") { if (nzchar(prompt)) cat(prompt); "" }
# base_dir <- "C:/path/to/parent/folder"

timestamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
out_dir <- file.path(results_root, paste0("DIANN_Replicate_CV_QC_", timestamp))   # output -> <TWO_PROTEASE_OUT>/04_reproducibility/<timestamp>
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "run_log.txt")
sink(log_file, split = TRUE)
cat("DIA-NN replicate CV QC started:", format(Sys.time()), "\n")
cat("Base dir:", base_dir, "\n")
cat("Output dir:", out_dir, "\n\n")
on.exit({
  cat("\nRun finished:", format(Sys.time()), "\n")
  sink()
}, add = TRUE)

# ==============================
# GENERAL HELPERS
# ==============================
pub_theme <- function(base_size = 17) {
  theme_classic(base_size = base_size) +
    theme(panel.spacing = unit(0.15, "lines"),
      plot.title = element_text(face = "bold", color = "black", size = base_size + 4),
      plot.subtitle = element_text(face = "bold", color = "black", size = base_size + 1),
      axis.title = element_text(face = "bold", color = "black", size = base_size + 2),
      axis.text = element_text(color = "black", face = "bold", size = base_size),
      strip.text = element_text(face = "bold", color = "black", size = base_size + 1),
      legend.title = element_text(face = "bold", color = "black", size = base_size),
      legend.text = element_text(face = "bold", color = "black", size = base_size - 1),
      panel.grid = element_blank(),
      strip.background = element_blank()
    )
}

save_pub <- function(plot_obj, out_png, width = 12, height = 7, dpi = PUB_DPI) {
  # Keep dimensions within practical bounds and use limitsize=FALSE as a safeguard.
  width <- min(max(width, 4), 48)
  height <- min(max(height, 4), 48)

  ggsave(
    out_png, plot_obj,
    width = width, height = height, dpi = dpi, units = "in",
    limitsize = FALSE
  )
  out_pdf <- sub("\\.png$", ".pdf", out_png, ignore.case = TRUE)
  ggsave(
    out_pdf, plot_obj,
    width = width, height = height, device = cairo_pdf, units = "in",
    limitsize = FALSE
  )
}

sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

fmt_int <- function(x) formatC(as.integer(x), format = "d", big.mark = ",")

natural_sort_key <- function(x) {
  # Produces numeric-aware ordering: S1, S2, ..., S9, S10 rather than S1, S10, S2.
  x <- as.character(x)
  prefix <- sub("([0-9]+).*", "", x)
  num <- suppressWarnings(as.integer(sub("^[^0-9]*([0-9]+).*", "\\1", x)))
  num[is.na(num)] <- .Machine$integer.max
  data.table(original = x, prefix = prefix, num = num)[order(prefix, num, original), original]
}

sampleday_order <- function(dt) {
  tmp <- unique(as.data.table(dt)[, c("Day", "SampleID"), with = FALSE])
  tmp[, DayNum := suppressWarnings(as.integer(Day))]
  tmp[, SampleNum := suppressWarnings(as.integer(SampleID))]
  tmp[is.na(DayNum), DayNum := .Machine$integer.max]
  tmp[is.na(SampleNum), SampleNum := .Machine$integer.max]
  tmp[, SampleDay := paste0("D", Day, "_S", SampleID)]
  tmp[order(DayNum, SampleNum, Day, SampleID), SampleDay]
}

add_sample_coverage <- function(dt, sample_meta, feature_cols = c("Search", "Level", "FeatureID")) {
  dt <- as.data.table(dt)
  sm <- as.data.table(sample_meta)
  total_samples <- uniqueN(sm$SampleID)
  total_sample_days <- uniqueN(paste(sm$SampleID, sm$Day, sep = "___"))

  # Inter-day CV tables are sample-level and do not have a Day column.
  if ("Day" %in% names(dt)) {
    cov <- dt[, list(
      SamplesDetected = uniqueN(SampleID),
      SampleDaysDetected = uniqueN(paste(SampleID, Day, sep = "___")),
      SamplesDetectedList = paste(natural_sort_key(unique(paste0("S", SampleID))), collapse = ";"),
      SampleDaysDetectedList = paste(sampleday_order(.SD), collapse = ";")
    ), by = feature_cols]
  } else {
    cov <- dt[, list(
      SamplesDetected = uniqueN(SampleID),
      SampleDaysDetected = NA_integer_,
      SamplesDetectedList = paste(natural_sort_key(unique(paste0("S", SampleID))), collapse = ";"),
      SampleDaysDetectedList = NA_character_
    ), by = feature_cols]
  }

  cov[, TotalSamplesInExperiment := total_samples]
  cov[, TotalSampleDaysInExperiment := total_sample_days]
  cov[, PercentSamplesDetected := 100 * SamplesDetected / TotalSamplesInExperiment]
  cov[, PercentSampleDaysDetected := ifelse(is.na(SampleDaysDetected), NA_real_, 100 * SampleDaysDetected / TotalSampleDaysInExperiment)]
  merge(dt, cov, by = feature_cols, all.x = TRUE)
}

identify_well_quantified_interday <- function(inter_cv, threshold, sample_meta) {
  dt <- copy(as.data.table(inter_cv))
  dt <- dt[
    MeasuredDays >= min_days_for_inter_cv &
      !is.na(CV) &
      CV <= threshold &
      CompleteDays %in% TRUE
  ]
  dt[, WellQuantifiedInterDay := TRUE]
  dt <- add_sample_coverage(dt, sample_meta)
  dt[]
}


get_basename_noext <- function(x) {
  x2 <- gsub("\\\\", "/", as.character(x))
  b <- basename(x2)
  b <- sub("\\.raw$", "", b, ignore.case = TRUE)
  b <- sub("\\.d$", "", b, ignore.case = TRUE)
  b <- sub("\\.mzML$", "", b, ignore.case = TRUE)
  b <- sub("\\.mzML\\.gz$", "", b, ignore.case = TRUE)
  b
}


repair_duplicate_names <- function(dt) {
  # DIA-NN can contain duplicate annotation column names after export/merge.
  # data.table allows them, but setnames/joins can behave unpredictably.
  nm <- names(dt)
  if (anyDuplicated(nm)) {
    names(dt) <- make.unique(nm, sep = "_dup")
  }
  dt
}

safe_get_first_column <- function(dt, candidates) {
  nms <- names(dt)
  for (cand in candidates) {
    hit <- which(nms == cand)
    if (length(hit) > 0) return(nms[hit[1]])
    # Also accept make.unique-renamed duplicate columns such as Protein.Group_dup1.
    hit2 <- grep(paste0("^", gsub("([.])", "\\\\\\1", cand), "(_dup[0-9]+)?$"), nms, value = TRUE)
    if (length(hit2) > 0) return(hit2[1])
  }
  NA_character_
}


sanitize_ids <- function(ids) {
  ids <- as.character(ids)
  ids <- trimws(ids)
  ids[ids %in% c("", "NA", "NaN")] <- NA_character_
  ids
}

normalize_colname <- function(x) tolower(gsub("[^a-z0-9]+", "", x))

is_non_sample_col <- function(colname) {
  raw <- tolower(trimws(as.character(colname)))
  nrm <- normalize_colname(colname)
  if (nrm %in% c(
    "proteingroup", "genes", "proteinids", "proteinnames", "firstprotein",
    "modifiedsequence", "modifiedpeptidesequence", "strippedsequence",
    "precursorid", "precursorcharge", "precursorionmobility",
    "proteotypic", "nsequences", "nsequence", "nseq", "nseqs",
    "qvalue", "pgqvalue", "pep", "pgpep", "libqvalue",
    "decoy", "decoyevidence", "rt", "rtstart", "rtstop",
    "ms1area", "ms1profilecorr"
  )) return(TRUE)
  if (grepl("^precursor\\.?\\s*charge$", raw)) return(TRUE)
  if (grepl("^proteotypic$", raw)) return(TRUE)
  if (grepl("^n\\.?\\s*sequences$", raw)) return(TRUE)
  FALSE
}

numeric_sample_cols <- function(dt, exclude = character()) {
  cols <- setdiff(names(dt), exclude)
  cols <- cols[!vapply(cols, is_non_sample_col, logical(1))]
  keep <- character()
  for (cname in cols) {
    v <- suppressWarnings(as.numeric(dt[[cname]]))
    if (sum(!is.na(v)) > 0) keep <- c(keep, cname)
  }
  unique(keep)
}

make_unique_search_names <- function(paths) {
  nms <- basename(paths)
  if (anyDuplicated(nms)) nms <- make.unique(nms, sep = "_")
  names(paths) <- nms
  paths
}

discover_search_folders <- function(parent_dir, pg_file, pr_file) {
  parent_pg <- file.path(parent_dir, pg_file)
  parent_pr <- file.path(parent_dir, pr_file)
  if (file.exists(parent_pg) && file.exists(parent_pr)) {
    out <- parent_dir
    names(out) <- basename(parent_dir)
    return(out)
  }
  subs <- list.dirs(parent_dir, full.names = TRUE, recursive = FALSE)
  ok <- vapply(subs, function(d) {
    file.exists(file.path(d, pg_file)) && file.exists(file.path(d, pr_file))
  }, logical(1))
  found <- subs[ok]
  if (length(found) == 0) {
    stop("No search folders found. Expected each folder to contain ", pg_file, " and ", pr_file)
  }
  found <- found[order(basename(found))]
  make_unique_search_names(found)
}

prompt_search_short_names <- function(search_paths_named) {
  cat("\nSearch strategies found:\n")
  print(data.table(Search = names(search_paths_named), Folder = unname(search_paths_named)))
  cat("\nShort names for search strategies. Press Enter to keep defaults.\n")
  new_names <- character(length(search_paths_named))
  for (i in seq_along(search_paths_named)) {
    default <- names(search_paths_named)[i]
    ans <- readline(sprintf("  %s -> short name [%s]: ", default, default))
    ans <- trimws(ans)
    if (!nzchar(ans)) ans <- default
    new_names[i] <- ans
  }
  names(search_paths_named) <- make.unique(new_names, sep = "_")
  search_paths_named
}

# ==============================
# SAMPLE METADATA
# ==============================
parse_sample_name <- function(sample_col) {
  clean <- get_basename_noext(sample_col)

  # Normalize common separators to underscores, then split.
  norm <- gsub("[^A-Za-z0-9]+", "_", clean)
  toks <- strsplit(norm, "_", fixed = TRUE)[[1]]
  toks <- toks[nzchar(toks)]

  day <- sample_id <- rep <- NA_character_
  parse_ok <- FALSE

  if (sample_parse_mode == "strict_first3") {
    if (length(toks) >= 3) {
      day <- toks[1]
      sample_id <- toks[2]
      rep <- toks[3]
      parse_ok <- grepl("^[0-9]+$", day) && grepl("^[0-9]+$", sample_id) && grepl("^[0-9]+$", rep)
    }
  } else {
    numeric_idx <- which(grepl("^[0-9]+$", toks))

    if (length(numeric_idx) >= 3) {
      if (sample_parse_mode == "first_numeric_triplet") {
        use_idx <- numeric_idx[1:3]
      } else {
        # Default: last numeric triplet. This is safest when file names contain
        # dates or run numbers before the true Day_SampleID_Replicate key.
        use_idx <- tail(numeric_idx, 3)
      }
      day <- toks[use_idx[1]]
      sample_id <- toks[use_idx[2]]
      rep <- toks[use_idx[3]]
      parse_ok <- TRUE
    }
  }

  data.table(
    Sample = sample_col,
    SampleClean = clean,
    ParsedTokens = paste(toks, collapse = "|"),
    Day = as.character(day),
    SampleID = as.character(sample_id),
    Replicate = as.character(rep),
    ParseOK = parse_ok
  )
}

build_sample_meta <- function(samples) {
  sm <- rbindlist(lapply(samples, parse_sample_name), fill = TRUE)
  sm <- unique(sm, by = "Sample")
  sm[, SampleDay := paste0("D", Day, "_S", SampleID)]
  sm[, SampleDayLabel := paste0("Day ", Day, " / Sample ", SampleID)]
  sm[, ReplicateLabel := paste0("R", Replicate)]
  sm[]
}

validate_sample_meta <- function(sample_meta) {
  bad <- sample_meta[is.na(Day) | is.na(SampleID) | is.na(Replicate) | !ParseOK]
  if (nrow(bad) > 0) {
    fwrite(bad, file.path(out_dir, "sample_name_parse_failures.csv"))
    msg <- paste0(
      "Some columns could not be parsed as Day_SampleID_Replicate. ",
      "They were written to sample_name_parse_failures.csv."
    )
    if (isTRUE(drop_unparsed_sample_columns)) {
      cat("WARNING: ", msg, "
", sep = "")
      cat("Dropping unparsed columns from CV analysis:
")
      print(bad[, c("Sample", "SampleClean", "ParsedTokens"), with = FALSE])
    } else {
      stop(msg)
    }
  }

  good <- sample_meta[!(is.na(Day) | is.na(SampleID) | is.na(Replicate) | !ParseOK)]
  dup <- good[, list(N = .N), by = c("Day", "SampleID", "Replicate")][N > 1]
  if (nrow(dup) > 0) {
    fwrite(dup, file.path(out_dir, "duplicate_day_sample_replicate_keys.csv"))
    cat("WARNING: Duplicate Day/SampleID/Replicate keys found. See duplicate_day_sample_replicate_keys.csv
")
  }
  invisible(TRUE)
}


# ==============================
# MATRIX IMPORT
# ==============================
matrix_from_dt <- function(dt, id_col, num_cols, zero_as_na = TRUE, collapse_method = "sum") {
  ids <- sanitize_ids(dt[[id_col]])
  keep <- !is.na(ids) & nzchar(ids)
  dt <- dt[keep]
  ids <- ids[keep]

  for (cname in num_cols) dt[[cname]] <- suppressWarnings(as.numeric(dt[[cname]]))
  m <- as.matrix(dt[, ..num_cols])
  mode(m) <- "numeric"
  if (zero_as_na) m[m == 0] <- NA_real_

  if (anyDuplicated(ids)) {
    agg <- data.table(FeatureID = ids)
    for (j in seq_along(num_cols)) agg[[num_cols[j]]] <- m[, j]
    if (collapse_method == "sum") {
      agg <- agg[, lapply(.SD, function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)),
                 by = "FeatureID", .SDcols = num_cols]
    } else if (collapse_method == "mean") {
      agg <- agg[, lapply(.SD, function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)),
                 by = "FeatureID", .SDcols = num_cols]
    } else {
      stop("Unknown collapse method: ", collapse_method)
    }
    rn <- agg[["FeatureID"]]
    agg[, FeatureID := NULL]
    m2 <- as.matrix(agg)
    mode(m2) <- "numeric"
    rownames(m2) <- rn
    colnames(m2) <- num_cols
    return(m2)
  }

  rownames(m) <- ids
  colnames(m) <- num_cols
  m
}

pick_first_existing <- function(candidates, names_vec) {
  hit <- candidates[candidates %in% names_vec]
  if (length(hit) == 0) NA_character_ else hit[1]
}




get_first_existing_col <- function(dt, candidates) {
  for (cc in candidates) {
    if (cc %in% names(dt)) return(cc)
  }
  NA_character_
}

extract_first_gene_value <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  vapply(x, function(xx) {
    if (!nzchar(xx) || xx == "NA") return(NA_character_)
    parts <- unlist(strsplit(xx, ";|,|\\s+", perl = TRUE))
    parts <- trimws(parts)
    parts <- parts[nzchar(parts) & parts != "NA"]
    if (length(parts) == 0) return(NA_character_)
    parts[1]
  }, character(1))
}

force_genes_column <- function(dt) {
  dt <- as.data.table(dt)

  # Create canonical Genes column if missing.
  if (!("Genes" %in% names(dt))) {
    dt[, Genes := NA_character_]
  } else {
    dt[, Genes := as.character(Genes)]
  }

  # Fill from Gene if present.
  if ("Gene" %in% names(dt)) {
    bad <- is.na(dt[["Genes"]]) | !nzchar(dt[["Genes"]]) | dt[["Genes"]] == "NA"
    if (any(bad)) dt[bad, Genes := as.character(get("Gene"))]
  }

  # Fill from repaired duplicate Genes columns if present.
  dup_gene_cols <- grep("^Genes(_dup[0-9]+)?$", names(dt), value = TRUE)
  dup_gene_cols <- setdiff(dup_gene_cols, "Genes")
  if (length(dup_gene_cols) > 0) {
    for (cc in dup_gene_cols) {
      bad <- is.na(dt[["Genes"]]) | !nzchar(dt[["Genes"]]) | dt[["Genes"]] == "NA"
      if (any(bad)) dt[bad, Genes := as.character(get(cc))]
    }
  }

  dt[, Genes := extract_first_gene_value(Genes)]
  dt[]
}

safe_annotation_table <- function(dt, id_col, ann_cols, level_name) {
  dt <- repair_duplicate_names(as.data.table(dt))
  dt <- force_genes_column(dt)

  id_idx <- match(id_col, names(dt))
  if (is.na(id_idx)) stop("ID column not found after duplicate-name repair: ", id_col)

  # Always preserve Genes explicitly.
  ann_cols <- unique(c("Genes", ann_cols))
  ann_cols <- ann_cols[ann_cols %in% names(dt)]
  ann_cols <- setdiff(ann_cols, id_col)

  keep_names <- unique(c(names(dt)[id_idx], ann_cols))
  ann <- dt[, ..keep_names]
  setnames(ann, names(dt)[id_idx], "FeatureID")

  ann[, FeatureID := sanitize_ids(FeatureID)]
  ann <- ann[!is.na(FeatureID) & nzchar(FeatureID)]

  ann <- ann[, lapply(.SD, function(x) {
    vals <- unique(na.omit(as.character(x)))
    vals <- vals[nzchar(vals) & vals != "NA"]
    paste(vals, collapse = ";")
  }), by = "FeatureID"]

  ann <- force_genes_column(ann)
  ann[, Level := level_name]
  ann[]
}


read_protein_matrix <- function(path) {
  dt <- fread(path, sep = "\t", header = TRUE, data.table = TRUE, showProgress = FALSE)
  dt <- repair_duplicate_names(dt)
  id_col <- safe_get_first_column(dt, c("Protein.Group", "Genes", names(dt)[1]))
  if (is.na(id_col)) stop("Could not infer protein ID column in ", path)

  ann_cols <- intersect(c("Genes", "Protein.Group", "Protein.Ids", "Protein.Names", "First.Protein.Description"), names(dt))
  num_cols <- numeric_sample_cols(dt, exclude = unique(c(id_col, ann_cols)))
  if (length(num_cols) == 0) stop("No numeric sample columns found in ", path)

  m <- matrix_from_dt(dt, id_col = id_col, num_cols = num_cols, zero_as_na = zero_as_na, collapse_method = "sum")

  ann <- safe_annotation_table(dt, id_col, ann_cols, "Protein")

  list(matrix = m, annotation = ann, sample_cols = num_cols)
}

read_peptide_matrix <- function(path) {
  dt <- fread(path, sep = "\t", header = TRUE, data.table = TRUE, showProgress = FALSE)
  dt <- repair_duplicate_names(dt)

  pep_col <- safe_get_first_column(dt, c(
    "Modified.Sequence",
    "ModifiedPeptideSequence",
    "Modified sequence",
    "Modified.Peptide.Sequence",
    "ModifiedPeptide.Sequence",
    "Stripped.Sequence"
  ))
  if (is.na(pep_col)) stop("Could not find modified peptide column in ", path)

  ann_cols <- intersect(c("Genes", "Protein.Group", "Protein.Ids", "Protein.Names",
                          "Stripped.Sequence", "Modified.Sequence", "ModifiedPeptideSequence",
                          "Precursor.Id", "PrecursorCharge", "Precursor.Charge"), names(dt))
  num_cols <- numeric_sample_cols(dt, exclude = unique(c(pep_col, ann_cols)))
  if (length(num_cols) == 0) stop("No numeric sample columns found in ", path)

  dt[[pep_col]] <- sanitize_ids(dt[[pep_col]])
  dt <- dt[!is.na(dt[[pep_col]]) & nzchar(dt[[pep_col]])]
  for (cname in num_cols) dt[[cname]] <- suppressWarnings(as.numeric(dt[[cname]]))
  if (zero_as_na) {
    for (cname in num_cols) dt[[cname]][dt[[cname]] == 0] <- NA_real_
  }

  if (peptide_collapse_method == "sum") {
    pep <- dt[, lapply(.SD, function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)),
              by = pep_col, .SDcols = num_cols]
  } else if (peptide_collapse_method == "mean") {
    pep <- dt[, lapply(.SD, function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)),
              by = pep_col, .SDcols = num_cols]
  } else {
    stop("Unknown peptide_collapse_method: ", peptide_collapse_method)
  }

  rn <- sanitize_ids(pep[[pep_col]])
  pep[[pep_col]] <- NULL
  m <- as.matrix(pep)
  mode(m) <- "numeric"
  rownames(m) <- rn
  colnames(m) <- num_cols

  ann <- safe_annotation_table(dt, pep_col, ann_cols, "Peptide")

  list(matrix = m, annotation = ann, sample_cols = num_cols)
}

align_matrix_samples <- function(m, all_samples) {
  missing <- setdiff(all_samples, colnames(m))
  if (length(missing) > 0) {
    add <- matrix(NA_real_, nrow = nrow(m), ncol = length(missing))
    rownames(add) <- rownames(m)
    colnames(add) <- missing
    m <- cbind(m, add)
  }
  m[, all_samples, drop = FALSE]
}

align_matrix_ids <- function(m, all_ids) {
  missing <- setdiff(all_ids, rownames(m))
  if (length(missing) > 0) {
    add <- matrix(NA_real_, nrow = length(missing), ncol = ncol(m))
    rownames(add) <- missing
    colnames(add) <- colnames(m)
    m <- rbind(m, add)
  }
  m[all_ids, , drop = FALSE]
}

build_all_searches_matrix <- function(mats_named) {
  all_ids <- sort(unique(unlist(lapply(mats_named, rownames))))
  all_samples <- sort(unique(unlist(lapply(mats_named, colnames))))

  mats <- lapply(mats_named, function(m) {
    m2 <- align_matrix_samples(m, all_samples)
    align_matrix_ids(m2, all_ids)
  })

  sum_m <- matrix(0, nrow = length(all_ids), ncol = length(all_samples))
  cnt_m <- matrix(0, nrow = length(all_ids), ncol = length(all_samples))
  rownames(sum_m) <- all_ids
  colnames(sum_m) <- all_samples
  rownames(cnt_m) <- all_ids
  colnames(cnt_m) <- all_samples

  for (nm in names(mats)) {
    m <- mats[[nm]]
    ok <- !is.na(m)
    sum_m[ok] <- sum_m[ok] + m[ok]
    cnt_m[ok] <- cnt_m[ok] + 1
  }

  out <- sum_m
  out[cnt_m == 0] <- NA_real_
  out[cnt_m > 0] <- sum_m[cnt_m > 0] / cnt_m[cnt_m > 0]
  out
}

combine_annotations <- function(ann_list) {
  ann <- rbindlist(ann_list, fill = TRUE)
  if (nrow(ann) == 0) return(data.table())
  ann <- ann[, lapply(.SD, function(x) paste(unique(na.omit(as.character(x))), collapse = ";")),
             by = c("Level", "FeatureID")]
  ann[]
}

# ==============================
# CV CALCULATIONS
# ==============================
cv_percent <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 2) return(NA_real_)
  mn <- mean(x)
  if (!is.finite(mn) || mn <= 0) return(NA_real_)
  100 * stats::sd(x) / mn
}

matrix_to_long <- function(m, search_name, level_name, sample_meta) {
  dt <- as.data.table(m, keep.rownames = "FeatureID")
  long <- melt(
    dt,
    id.vars = "FeatureID",
    variable.name = "Sample",
    value.name = "Intensity",
    variable.factor = FALSE
  )
  long[, Search := search_name]
  long[, Level := level_name]
  long <- merge(long, sample_meta, by = "Sample", all.x = TRUE)
  long[]
}

calc_intra_day_cv <- function(m, search_name, level_name, sample_meta) {
  long <- matrix_to_long(m, search_name, level_name, sample_meta)
  expected <- unique(long[, c("SampleID", "Day", "Replicate"), with = FALSE])
  expected_counts <- expected[, list(ExpectedReplicates = uniqueN(Replicate)), by = c("SampleID", "Day")]

  out <- long[, list(
    MeanIntensity = if (all(is.na(Intensity))) NA_real_ else mean(Intensity, na.rm = TRUE),
    SDIntensity = if (sum(is.finite(Intensity) & Intensity > 0) < 2) NA_real_ else stats::sd(Intensity[is.finite(Intensity) & Intensity > 0]),
    CV = cv_percent(Intensity),
    MeasuredReplicates = sum(is.finite(Intensity) & Intensity > 0),
    TotalReplicateRows = .N,
    ReplicatesMeasured = paste(sort(unique(Replicate[is.finite(Intensity) & Intensity > 0])), collapse = ";")
  ), by = c("Search", "Level", "FeatureID", "SampleID", "Day")]

  out <- merge(out, expected_counts, by = c("SampleID", "Day"), all.x = TRUE)
  out[, CompleteReplicates := MeasuredReplicates == ExpectedReplicates]
  out[, CV_Type := "Intra-day"]
  setcolorder(out, c("CV_Type", "Search", "Level", "FeatureID", "SampleID", "Day",
                     "CV", "MeanIntensity", "SDIntensity", "MeasuredReplicates",
                     "ExpectedReplicates", "CompleteReplicates", "ReplicatesMeasured"))
  out[]
}

calc_inter_day_cv <- function(m, search_name, level_name, sample_meta) {
  long <- matrix_to_long(m, search_name, level_name, sample_meta)

  day_means <- long[, list(
    DayMeanIntensity = if (all(is.na(Intensity))) NA_real_ else mean(Intensity, na.rm = TRUE),
    MeasuredReplicates = sum(is.finite(Intensity) & Intensity > 0),
    ExpectedReplicates = uniqueN(Replicate)
  ), by = c("Search", "Level", "FeatureID", "SampleID", "Day")]
  day_means[, CompleteReplicates := MeasuredReplicates == ExpectedReplicates]

  expected_days <- unique(sample_meta[, c("SampleID", "Day"), with = FALSE])
  expected_day_counts <- expected_days[, list(ExpectedDays = uniqueN(Day)), by = "SampleID"]

  out <- day_means[, list(
    MeanIntensity = if (all(is.na(DayMeanIntensity))) NA_real_ else mean(DayMeanIntensity, na.rm = TRUE),
    SDIntensity = if (sum(is.finite(DayMeanIntensity) & DayMeanIntensity > 0) < 2) NA_real_ else stats::sd(DayMeanIntensity[is.finite(DayMeanIntensity) & DayMeanIntensity > 0]),
    CV = cv_percent(DayMeanIntensity),
    MeasuredDays = sum(is.finite(DayMeanIntensity) & DayMeanIntensity > 0),
    CompleteReplicateDays = sum(CompleteReplicates %in% TRUE),
    DaysMeasured = paste(sort(unique(Day[is.finite(DayMeanIntensity) & DayMeanIntensity > 0])), collapse = ";")
  ), by = c("Search", "Level", "FeatureID", "SampleID")]

  out <- merge(out, expected_day_counts, by = "SampleID", all.x = TRUE)
  out[, CompleteDays := MeasuredDays == ExpectedDays]
  out[, CV_Type := "Inter-day"]
  setcolorder(out, c("CV_Type", "Search", "Level", "FeatureID", "SampleID",
                     "CV", "MeanIntensity", "SDIntensity", "MeasuredDays",
                     "ExpectedDays", "CompleteDays", "CompleteReplicateDays", "DaysMeasured"))
  out[]
}

identify_well_quantified_intra <- function(intra_cv) {
  dt <- copy(intra_cv)
  dt <- dt[MeasuredReplicates >= min_replicates_for_intra_cv]
  dt <- dt[!is.na(CV) & CV <= intra_cv_threshold]
  if (require_all_replicates_for_well_quantified) {
    dt <- dt[CompleteReplicates %in% TRUE]
  }
  dt[, WellQuantifiedIntra := TRUE]
  dt[]
}

identify_well_quantified_consensus <- function(wq_intra, inter_cv) {
  if (nrow(wq_intra) == 0 || nrow(inter_cv) == 0) return(data.table())

  expected_days <- unique(inter_cv[, c("Search", "Level", "FeatureID", "SampleID", "ExpectedDays"), with = FALSE])
  pass_days <- wq_intra[, list(PassingDays = uniqueN(Day),
                               DaysPassing = paste(sort(unique(Day)), collapse = ";")),
                        by = c("Search", "Level", "FeatureID", "SampleID")]
  pass_days <- merge(pass_days, expected_days, by = c("Search", "Level", "FeatureID", "SampleID"), all.x = TRUE)

  inter_pass <- inter_cv[
    MeasuredDays >= min_days_for_inter_cv & !is.na(CV) & CV <= inter_cv_threshold
  ]
  if (require_all_days_for_consensus) {
    inter_pass <- inter_pass[CompleteDays %in% TRUE]
    pass_days <- pass_days[PassingDays == ExpectedDays]
  }

  out <- merge(
    pass_days,
    inter_pass[, c("Search", "Level", "FeatureID", "SampleID", "CV", "MeanIntensity", "MeasuredDays", "ExpectedDays", "CompleteDays"), with = FALSE],
    by = c("Search", "Level", "FeatureID", "SampleID"),
    all = FALSE,
    suffixes = c("_IntraSummary", "_InterDay")
  )
  setnames(out, "CV", "InterDayCV")
  out[, WellQuantifiedConsensus := TRUE]
  out[]
}


standardize_protein_gene_columns <- function(dt) {
  dt <- as.data.table(dt)

  # Create robust Protein column from any available DIA-NN-style protein annotation.
  if (!("Protein" %in% names(dt))) {
    protein_candidates <- c("Protein.Group", "Protein.Ids", "Protein.Names", "First.Protein.Description")
    hit <- protein_candidates[protein_candidates %in% names(dt)]
    if (length(hit) > 0) {
      dt[, Protein := as.character(get(hit[1]))]
    } else {
      dt[, Protein := NA_character_]
    }
  }

  # Create robust Gene column.
  if (!("Gene" %in% names(dt))) {
    if ("Genes" %in% names(dt)) {
      dt[, Gene := as.character(Genes)]
    } else {
      dt[, Gene := NA_character_]
    }
  }

  dt[]
}


add_annotations <- function(dt, ann) {
  dt <- as.data.table(dt)
  ann <- as.data.table(ann)
  if (nrow(dt) == 0 || nrow(ann) == 0) return(dt)

  # Remove stale annotation columns before adding canonical ones back.
  stale <- intersect(c("Genes", "Gene", "Protein", "Protein.Group", "Protein.Ids", "Protein.Names"), names(dt))
  if (length(stale) > 0) dt[, (stale) := NULL]

  out <- merge(dt, ann, by = c("Level", "FeatureID"), all.x = TRUE)
  out <- standardize_protein_gene_columns(out)
  out[]
}

# ==============================
# PTM EXTRACTION
# ==============================
unimod_name_map <- c(
  "1" = "Acetyl",
  "4" = "Carbamidomethyl",
  "5" = "Carbamyl",
  "7" = "Deamidated",
  "21" = "Phospho",
  "27" = "Glu->pyro-Glu",
  "28" = "Gln->pyro-Glu",
  "35" = "Oxidation",
  "36" = "Dimethyl",
  "37" = "Trimethyl",
  "40" = "Sulfo",
  "121" = "GlyGly",
  "259" = "Label:13C(6)15N(2)",
  "267" = "Label:13C(6)15N(4)",
  "737" = "TMT6plex",
  "2016" = "TMTpro"
)

extract_unimod_ids <- function(x) {
  x <- as.character(x)
  hits <- gregexpr("UniMod:([0-9]+)", x, perl = TRUE)
  reg <- regmatches(x, hits)[[1]]
  if (length(reg) == 0 || identical(reg, character(0))) return(character(0))
  unique(sub("UniMod:", "", reg))
}

extract_mass_mods <- function(x) {
  x <- as.character(x)
  hits <- gregexpr("\\([+-]?[0-9]+\\.?[0-9]*\\)", x, perl = TRUE)
  reg <- regmatches(x, hits)[[1]]
  if (length(reg) == 0 || identical(reg, character(0))) return(character(0))
  unique(gsub("[()]", "", reg))
}

describe_mods <- function(feature_id) {
  ids <- extract_unimod_ids(feature_id)
  if (length(ids) > 0) {
    names <- unimod_name_map[ids]
    names[is.na(names)] <- paste0("UniMod:", ids[is.na(names)])
    return(paste(unique(names), collapse = ";"))
  }
  mass <- extract_mass_mods(feature_id)
  if (length(mass) > 0) return(paste0("MassDelta:", paste(mass, collapse = ";MassDelta:")))
  "Unmodified_or_unknown"
}

make_ptm_table <- function(wq_peptides_annotated) {
  if (nrow(wq_peptides_annotated) == 0) return(data.table())
  dt <- copy(wq_peptides_annotated)
  dt <- standardize_protein_gene_columns(dt)
  dt[, UniModIDs := vapply(FeatureID, function(x) paste(extract_unimod_ids(x), collapse = ";"), character(1))]
  dt[, PTM := vapply(FeatureID, describe_mods, character(1))]
  dt[]
}

# ==============================
# PLOTS
# ==============================
plot_cv_violin_one_search <- function(cv_dt, search_name, level_name, cv_type, out_png) {
  dt <- cv_dt[Search == search_name & Level == level_name & CV_Type == cv_type & is.finite(CV)]
  if (nrow(dt) == 0) return(invisible(NULL))

  # Use identical search-strategy color mapping as the well_quantified_inter_day_protein count heatmap.
  search_levels <- get_search_levels(cv_dt$Search)
  pal <- search_palette(search_levels)
  dt[, Search := factor(as.character(Search), levels = search_levels)]

  if (cv_type == "Intra-day") {
    dt[, XLabel := paste0("D", Day, "_S", SampleID)]
    x_lab <- "Day / Sample"
  } else {
    dt[, XLabel := paste0("S", SampleID)]
    x_lab <- "Sample"
  }

  if (cv_type == "Intra-day") {
    x_levels <- sampleday_order(dt)
  } else {
    tmp_order <- unique(dt[, .(SampleID)])
    tmp_order[, SampleNum := suppressWarnings(as.integer(SampleID))]
    tmp_order[is.na(SampleNum), SampleNum := .Machine$integer.max]
    x_levels <- paste0("S", tmp_order[order(SampleNum, SampleID), SampleID])
  }
  dt[, XLabel := factor(XLabel, levels = x_levels)]

  base_size <- if (length(x_levels) > large_plot_label_threshold) 10 else 13

  p <- ggplot(dt, aes(x = XLabel, y = CV, fill = Search)) +
    geom_violin(trim = TRUE, scale = "width", color = "black", linewidth = 0.3, alpha = 0.80) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.25, alpha = 0.90) +
    geom_hline(yintercept = if (cv_type == "Intra-day") intra_cv_threshold else inter_cv_threshold,
               linetype = "dashed", linewidth = 0.5) +
    scale_fill_manual(values = pal, drop = FALSE) +
    labs(
      title = paste(search_name, "-", level_name, cv_type, "CV distribution"),
      x = x_lab,
      y = "CV (%)",
      fill = "Search"
    ) +
    pub_theme(base_size = base_size) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  save_pub(p, out_png, width = max(9, length(x_levels) * 0.45), height = 6.5)
  invisible(p)
}

plot_all_cv_violins <- function(intra_cv_all, inter_cv_all, searches_with_all, out_dir) {
  plot_dir <- file.path(out_dir, "CV_Violin_Plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  for (s in searches_with_all) {
    for (lvl in c("Protein", "Peptide")) {
      plot_cv_violin_one_search(
        intra_cv_all, s, lvl, "Intra-day",
        file.path(plot_dir, paste0("IntraDay_CV_Violin_", sanitize_filename(lvl), "_", sanitize_filename(s), ".png"))
      )
      plot_cv_violin_one_search(
        inter_cv_all, s, lvl, "Inter-day",
        file.path(plot_dir, paste0("InterDay_CV_Violin_", sanitize_filename(lvl), "_", sanitize_filename(s), ".png"))
      )
    }
  }
}

plot_wq_count_heatmap <- function(wq_protein, out_png) {
  dt0 <- as.data.table(wq_protein)
  if (nrow(dt0) == 0) return(invisible(NULL))

  if ("Day" %in% names(dt0)) {
    dt <- dt0[, list(WellQuantifiedProteins = uniqueN(FeatureID)), by = c("Search", "SampleID", "Day")]
    dt[, XLabel := paste0("D", Day, "_S", SampleID)]
    x_levels <- sampleday_order(dt)
    x_title <- "Day / Sample"
  } else {
    dt <- dt0[, list(WellQuantifiedProteins = uniqueN(FeatureID)), by = c("Search", "SampleID")]
    tmp <- unique(dt[, .(SampleID)])
    tmp[, SampleNum := suppressWarnings(as.integer(SampleID))]
    tmp[is.na(SampleNum), SampleNum := .Machine$integer.max]
    x_levels <- paste0("S", tmp[order(SampleNum, SampleID), SampleID])
    dt[, XLabel := paste0("S", SampleID)]
    x_title <- "Sample"
  }

  dt[, XLabel := factor(XLabel, levels = x_levels)]
  search_levels <- get_search_levels(dt$Search)
  dt[, Search := factor(as.character(Search), levels = search_levels)]
  pal <- search_palette(search_levels)

  p <- ggplot(dt, aes(x = XLabel, y = Search, fill = Search, alpha = WellQuantifiedProteins)) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = fmt_int(WellQuantifiedProteins)), size = 6.4, color = "black", fontface = "bold") +
    scale_fill_manual(values = pal) +
    scale_alpha_continuous(range = c(0.55, 1.0), guide = "none") +
    labs(
      title = "Inter-day low-CV proteins detected across all replicate sets",
      x = x_title,
      y = "Search strategy",
      fill = "Search"
    ) +
    pub_theme(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
          axis.text.y = element_text(face = "bold", color = "black"),
          legend.position = "bottom")

  save_pub(p, out_png, width = max(10, length(unique(dt$XLabel)) * 0.50), height = max(4.8, length(unique(dt$Search)) * 0.58))
  invisible(p)
}

make_pairwise_overlap <- function(dt_sets, id_col = "FeatureID", group_col = "Search") {
  groups <- sort(unique(dt_sets[[group_col]]))
  out <- list()
  k <- 1
  for (a in groups) {
    set_a <- unique(dt_sets[get(group_col) == a][[id_col]])
    for (b in groups) {
      set_b <- unique(dt_sets[get(group_col) == b][[id_col]])
      inter_n <- length(intersect(set_a, set_b))
      union_n <- length(union(set_a, set_b))
      pct_jaccard <- if (union_n == 0) NA_real_ else 100 * inter_n / union_n
      pct_of_a <- if (length(set_a) == 0) NA_real_ else 100 * inter_n / length(set_a)
      out[[k]] <- data.table(
        RowSearch = a,
        ColSearch = b,
        Intersection = inter_n,
        Union = union_n,
        JaccardPercent = pct_jaccard,
        PercentOfRowSearch = pct_of_a,
        RowSearchTotal = length(set_a),
        ColSearchTotal = length(set_b)
      )
      k <- k + 1
    }
  }
  rbindlist(out)
}

plot_pairwise_overlap_heatmap <- function(overlap_dt, out_png, value_col = "JaccardPercent",
                                          title = "Search-strategy overlap") {
  if (nrow(overlap_dt) == 0) return(invisible(NULL))
  dt <- copy(overlap_dt)
  dt[, RowSearch := factor(RowSearch, levels = rev(sort(unique(as.character(RowSearch)))))]
  dt[, ColSearch := factor(ColSearch, levels = sort(unique(as.character(ColSearch))))]
  dt[, Label := paste0(fmt_int(Intersection), "\n", sprintf("%.1f%%", get(value_col)))]

  p <- ggplot(dt, aes(x = ColSearch, y = RowSearch, fill = get(value_col))) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = Label), size = 3.3, color = "black", lineheight = 0.9) +
    scale_fill_gradient(low = "white", high = "grey30", na.value = "white") +
    labs(title = title, x = NULL, y = NULL, fill = "%") +
    pub_theme(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  n <- length(unique(dt$ColSearch))
  save_pub(p, out_png, width = max(6, n * 1.1), height = max(5, n * 0.8))
  invisible(p)
}

make_membership_patterns <- function(wq_intra_protein_no_all) {
  if (nrow(wq_intra_protein_no_all) == 0) return(data.table())
  keys <- unique(wq_intra_protein_no_all[, c("SampleID", "FeatureID"), with = FALSE])
  searches <- sort(unique(wq_intra_protein_no_all$Search))

  for (s in searches) {
    present <- unique(wq_intra_protein_no_all[Search == s, c("SampleID", "FeatureID"), with = FALSE])
    present[, value := TRUE]
    setnames(present, "value", s)
    keys <- merge(keys, present, by = c("SampleID", "FeatureID"), all.x = TRUE)
    keys[is.na(get(s)), (s) := FALSE]
  }

  keys[, SearchPattern := apply(.SD, 1, function(row) {
    hit <- searches[as.logical(row)]
    if (length(hit) == 0) "None" else paste(hit, collapse = " + ")
  }), .SDcols = searches]

  keys[]
}

plot_membership_pattern_heatmap <- function(pattern_dt, out_png) {
  if (nrow(pattern_dt) == 0) return(invisible(NULL))
  dt <- pattern_dt[, list(Proteins = uniqueN(FeatureID)), by = c("SampleID", "SearchPattern")]
  dt[, XLabel := paste0("S", SampleID)]
  tmp <- unique(dt[, .(SampleID)])
  tmp[, SampleNum := suppressWarnings(as.integer(SampleID))]
  tmp[is.na(SampleNum), SampleNum := .Machine$integer.max]
  dt[, XLabel := factor(XLabel, levels = paste0("S", tmp[order(SampleNum, SampleID), SampleID]))]
  dt[, SearchPattern := factor(SearchPattern, levels = unique(dt[order(SearchPattern), SearchPattern]))]

  p <- ggplot(dt, aes(x = XLabel, y = SearchPattern, fill = Proteins)) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = fmt_int(Proteins)), size = 3.0, color = "black") +
    scale_fill_gradient(low = "white", high = "grey30") +
    labs(
      title = "Inter-day low-CV protein search-strategy membership patterns",
      x = "Sample",
      y = "Search strategy pattern",
      fill = "Proteins"
    ) +
    pub_theme(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_pub(p, out_png,
           width = max(10, length(unique(dt$XLabel)) * 0.45),
           height = max(5, length(unique(dt$SearchPattern)) * 0.35))
  invisible(p)
}



# ==============================
# CONSISTENT COLOR + LABEL HELPERS
# ==============================

get_search_levels <- function(x) {
  x <- as.character(x)
  # Keep current script order, with All Searches last if present.
  lv <- unique(x)
  if (ALL_SEARCHES_LABEL %in% lv) {
    lv <- c(setdiff(lv, ALL_SEARCHES_LABEL), ALL_SEARCHES_LABEL)
  }
  lv
}


search_palette <- function(searches) {
  searches <- get_search_levels(searches)
  base_cols <- c(
    "#0072B2", # blue
    "#D55E00", # vermillion
    "#009E73", # green
    "#CC79A7", # purple/pink
    "#E69F00", # orange
    "#56B4E9", # sky blue
    "#F0E442", # yellow
    "#999999", # grey
    "#332288", "#88CCEE", "#44AA99", "#117733",
    "#DDCC77", "#CC6677", "#882255", "#AA4499"
  )
  cols <- rep(base_cols, length.out = length(searches))
  names(cols) <- searches
  if (ALL_SEARCHES_LABEL %in% searches) cols[ALL_SEARCHES_LABEL] <- "#8DA0CB"
  cols
}

get_feature_display_name <- function(dt) {
  dt <- as.data.table(dt)

  gene <- if ("Gene" %in% names(dt)) as.character(dt$Gene) else rep(NA_character_, nrow(dt))
  protein <- if ("Protein" %in% names(dt)) as.character(dt$Protein) else rep(NA_character_, nrow(dt))
  fid <- as.character(dt$FeatureID)

  gene[is.na(gene) | !nzchar(gene)] <- NA_character_
  protein[is.na(protein) | !nzchar(protein)] <- NA_character_

  label <- ifelse(!is.na(gene), paste0(gene, " / ", fid),
                  ifelse(!is.na(protein), paste0(protein, " / ", fid), fid))
  label
}


get_gene_display_name <- function(dt) {
  dt <- as.data.table(dt)
  dt <- force_genes_column(dt)
  fid <- if ("FeatureID" %in% names(dt)) as.character(dt$FeatureID) else rep(NA_character_, nrow(dt))
  gene <- as.character(dt$Genes)
  bad <- is.na(gene) | !nzchar(gene) | gene == "NA"
  gene[bad] <- fid[bad]
  gene
}






excluded_ptm_regex <- paste(
  c(
    "^Unmodified_or_unknown$",
    "Carbamidomethyl",
    "pyro-Glu",
    "Pyro",
    "Deamidated",
    "Deamidation",
    "Amidated",
    "Amidation",
    "Oxidation",
    "MassDelta:57\\.0",       # common carbamidomethyl-style mass
    "MassDelta:15\\.9",       # oxidation-ish
    "MassDelta:0\\.98",       # deamidation-ish
    "UniMod:2",
    "UniMod:23",
    "UniMod:30"
  ),
  collapse = "|"
)

filter_reportable_ptms <- function(wq_ptm) {
  dt <- as.data.table(wq_ptm)
  if (nrow(dt) == 0) return(dt)
  dt <- dt[!is.na(PTM) & nzchar(PTM)]
  dt <- dt[!grepl(excluded_ptm_regex, PTM, ignore.case = TRUE)]
  dt[]
}


paginate_dt <- function(dt, rows_per_page = 80) {
  dt <- as.data.table(dt)
  if (nrow(dt) == 0) return(list(dt))
  split(dt, ceiling(seq_len(nrow(dt)) / rows_per_page))
}

page_path <- function(out_png, page_i, n_pages) {
  if (n_pages <= 1) return(out_png)
  sub("\\.png$", paste0("_page", sprintf("%02d", page_i), ".png"), out_png, ignore.case = TRUE)
}




# ==============================
# FASTA + FULL-PROTEIN PTM SITE MAPPING
# ==============================
pick_fasta_file <- function() {
  if (!is.null(fasta_file) && nzchar(as.character(fasta_file)) && file.exists(fasta_file)) {
    return(fasta_file)
  }

  if (interactive() && .Platform$OS.type == "windows") {
    ans <- trimws(readline("Select FASTA for full-protein PTM site mapping? [Y/n]: "))
    if (!tolower(ans) %in% c("n", "no")) {
      p <- utils::choose.files(
        caption = "Select FASTA file for protein sequence mapping",
        filters = matrix(c("FASTA files", "*.fasta;*.fa;*.faa;*.fsa", "All files", "*.*"), ncol = 2, byrow = TRUE)
      )
      if (length(p) > 0 && file.exists(p[1])) return(p[1])
    }
  }

  NA_character_
}

read_fasta_sequences <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) {
    return(list(seq_by_accession = list(), seq_by_header = list()))
  }

  lines <- readLines(path, warn = FALSE)
  headers <- grep("^>", lines)
  if (length(headers) == 0) {
    return(list(seq_by_accession = list(), seq_by_header = list()))
  }

  seq_by_accession <- list()
  seq_by_header <- list()

  for (i in seq_along(headers)) {
    start <- headers[i]
    end <- if (i < length(headers)) headers[i + 1] - 1L else length(lines)
    header <- sub("^>", "", lines[start])
    seq_lines <- lines[(start + 1L):end]
    seq <- paste0(seq_lines[!grepl("^>", seq_lines)], collapse = "")
    seq <- toupper(gsub("[^A-Z]", "", seq))
    if (!nzchar(seq)) next

    accs <- extract_accessions_from_protein_string(header)
    if (length(accs) == 0) accs <- header

    seq_by_header[[header]] <- seq
    for (acc in accs) {
      if (!acc %in% names(seq_by_accession)) seq_by_accession[[acc]] <- seq
    }
  }

  list(seq_by_accession = seq_by_accession, seq_by_header = seq_by_header)
}

extract_accessions_from_protein_string <- function(x) {
  x <- as.character(x)
  if (length(x) == 0 || is.na(x[1])) return(character(0))
  pieces <- unlist(strsplit(x, ";|,| ", perl = TRUE))
  pieces <- trimws(pieces)
  pieces <- pieces[nzchar(pieces) & pieces != "NA"]

  out <- character()
  for (p in pieces) {
    # UniProt style sp|P12345|GENE or tr|A0A...|NAME
    if (grepl("\\|", p)) {
      pp <- strsplit(p, "\\|", perl = TRUE)[[1]]
      pp <- pp[nzchar(pp)]
      if (length(pp) >= 2) out <- c(out, pp[2])
      out <- c(out, pp)
    } else {
      out <- c(out, p)
    }
  }

  unique(out[nzchar(out)])
}

get_candidate_protein_accessions <- function(row_dt) {
  row_dt <- as.data.table(row_dt)
  candidates <- character()

  cols <- c("Protein", "Protein.Group", "Protein.Ids", "Protein.Names", "Genes")
  cols <- cols[cols %in% names(row_dt)]
  for (cc in cols) {
    candidates <- c(candidates, extract_accessions_from_protein_string(row_dt[[cc]][1]))
  }

  unique(candidates[nzchar(candidates)])
}

strip_modified_sequence_to_peptide <- function(modseq) {
  s <- as.character(modseq)[1]
  if (is.na(s) || !nzchar(s)) return("")

  # Remove bracketed/parenthesized modification annotations, then keep amino-acid letters.
  s <- gsub("\\([^)]*\\)", "", s, perl = TRUE)
  s <- gsub("\\[[^]]*\\]", "", s, perl = TRUE)
  s <- toupper(gsub("[^A-Z]", "", s))
  s
}

find_peptide_start_in_fasta <- function(peptide, row_dt, fasta_db) {
  peptide <- toupper(as.character(peptide)[1])
  if (is.na(peptide) || !nzchar(peptide)) return(NA_integer_)
  if (length(fasta_db$seq_by_accession) == 0 && length(fasta_db$seq_by_header) == 0) return(NA_integer_)

  accs <- get_candidate_protein_accessions(row_dt)

  # First try candidate protein accessions from DIA-NN annotation.
  for (acc in accs) {
    if (acc %in% names(fasta_db$seq_by_accession)) {
      seq <- fasta_db$seq_by_accession[[acc]]
      pos <- regexpr(peptide, seq, fixed = TRUE)[1]
      if (pos > 0) return(as.integer(pos))
    }
  }

  # Fallback: search all FASTA entries. This is slower but robust.
  for (seq in fasta_db$seq_by_accession) {
    pos <- regexpr(peptide, seq, fixed = TRUE)[1]
    if (pos > 0) return(as.integer(pos))
  }

  NA_integer_
}


extract_unimod_positions_from_modified_sequence <- function(modseq) {
  # Returns a data.table with peptide-local positions, residues, and UniMod IDs.
  # This function uses only explicit data.table column references downstream.
  s <- as.character(modseq)[1]
  if (is.na(s) || !nzchar(s)) {
    return(data.table(PeptidePosition = integer(), Residue = character(), UniModID = character()))
  }

  chars <- strsplit(s, "", fixed = TRUE)[[1]]
  aa_pos <- 0L
  last_aa <- NA_character_
  pos_vec <- integer()
  res_vec <- character()
  id_vec <- character()

  i <- 1L
  while (i <= length(chars)) {
    ch <- chars[i]

    if (grepl("^[A-Z]$", ch)) {
      aa_pos <- aa_pos + 1L
      last_aa <- ch
      i <- i + 1L
      next
    }

    if (ch %in% c("(", "[")) {
      close_ch <- if (ch == "(") ")" else "]"
      j <- i + 1L
      token_chars <- character()
      while (j <= length(chars) && chars[j] != close_ch) {
        token_chars <- c(token_chars, chars[j])
        j <- j + 1L
      }
      token <- paste0(token_chars, collapse = "")
      ids <- regmatches(token, gregexpr("UniMod:[0-9]+", token, perl = TRUE))[[1]]

      if (length(ids) > 0 && !identical(ids, character(0)) && !is.na(last_aa)) {
        for (id in ids) {
          pos_vec <- c(pos_vec, aa_pos)
          res_vec <- c(res_vec, last_aa)
          id_vec <- c(id_vec, id)
        }
      }

      i <- max(j + 1L, i + 1L)
      next
    }

    i <- i + 1L
  }

  if (length(pos_vec) == 0) {
    return(data.table(PeptidePosition = integer(), Residue = character(), UniModID = character()))
  }

  data.table(
    PeptidePosition = as.integer(pos_vec),
    Residue = as.character(res_vec),
    UniModID = as.character(id_vec)
  )
}

find_position_column <- function(dt) {
  candidates <- c(
    "Protein.Position", "ProteinPosition", "Protein.Positions", "ProteinPositions",
    "Position", "Positions", "PTM.Site", "PTMSite", "Site", "Sites",
    "Global.Position", "GlobalPosition", "Protein.Start", "ProteinStart",
    "Start.Position", "StartPosition"
  )
  hit <- candidates[candidates %in% names(dt)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

make_modsite_labels_for_row <- function(feature_id, row_dt) {
  row_dt <- as.data.table(row_dt)
  sites <- extract_unimod_positions_from_modified_sequence(feature_id)

  if (nrow(sites) == 0) {
    ids <- extract_unimod_ids(feature_id)
    if (length(ids) == 0) return("")
    return(paste(paste0("UNMAPPED U:", ids), collapse = ";"))
  }

  pos_col <- find_position_column(row_dt)
  protein_start <- NA_integer_

  # 1) Prefer explicit DIA-NN protein position/site column if present.
  if (!is.na(pos_col) && pos_col %in% names(row_dt)) {
    raw_split <- unlist(strsplit(as.character(row_dt[[pos_col]][1]), ";|,", perl = TRUE))
    raw_pos <- suppressWarnings(as.integer(raw_split[1]))
    if (length(raw_pos) > 0 && is.finite(raw_pos)) protein_start <- raw_pos
  }

  # 2) Otherwise map stripped modified peptide to the full FASTA protein sequence.
  if (!is.finite(protein_start)) {
    peptide <- strip_modified_sequence_to_peptide(feature_id)
    protein_start <- find_peptide_start_in_fasta(peptide, row_dt, fasta_db)
  }

  pep_pos <- as.integer(sites[["PeptidePosition"]])
  residues <- as.character(sites[["Residue"]])
  unimods <- as.character(sites[["UniModID"]])

  if (length(pep_pos) == 0 || length(residues) == 0 || length(unimods) == 0) return("")

  if (is.finite(protein_start)) {
    protein_pos <- protein_start + pep_pos - 1L
    labels <- paste0(residues, protein_pos, " ", sub("UniMod:", "U:", unimods, fixed = TRUE))
  } else {
    # No peptide-local fallback in final label; keep it visibly unmapped.
    labels <- paste0("UNMAPPED_", residues, " ", sub("UniMod:", "U:", unimods, fixed = TRUE))
  }

  paste(unique(labels), collapse = ";")
}

summarize_ptms_by_protein <- function(wq_ptm) {
  dt <- as.data.table(wq_ptm)
  if (nrow(dt) == 0) {
    return(data.table(Search = character(), GeneLabel = character(),
                      PTM_Sites = character(), PTM_Count = integer(), PTM_Features = character()))
  }

  dt <- filter_reportable_ptms(dt)
  if (nrow(dt) == 0) {
    return(data.table(Search = character(), GeneLabel = character(),
                      PTM_Sites = character(), PTM_Count = integer(), PTM_Features = character()))
  }

  dt <- standardize_protein_gene_columns(dt)
  dt[, GeneLabel := get_gene_display_name(.SD), .SDcols = names(dt)]

  labels <- character(nrow(dt))
  for (ii in seq_len(nrow(dt))) {
    labels[ii] <- make_modsite_labels_for_row(dt[["FeatureID"]][ii], dt[ii])
  }
  dt[, ModSiteLabel := labels]
  dt[, ModSiteLabel := gsub("UniMod:", "U:", ModSiteLabel, fixed = TRUE)]
  dt <- dt[nzchar(ModSiteLabel)]

  if (nrow(dt) == 0) {
    return(data.table(Search = character(), GeneLabel = character(),
                      PTM_Sites = character(), PTM_Count = integer(), PTM_Features = character()))
  }

  out <- dt[, list(
    PTM_Count = uniqueN(ModSiteLabel),
    PTM_Sites = paste(sort(unique(ModSiteLabel)), collapse = "; "),
    PTM_Features = paste(sort(unique(FeatureID)), collapse = " | ")
  ), by = c("Search", "GeneLabel")]

  out[]
}

attach_ptms_to_protein_coverage <- function(protein_cov, wq_ptm) {
  pc <- as.data.table(protein_cov)
  if (nrow(pc) == 0) return(pc)
  pc <- standardize_protein_gene_columns(pc)
  pc[, GeneLabel := get_gene_display_name(.SD), .SDcols = names(pc)]
  pc[, PTM_Count := 0L]
  pc[, PTM_Sites := ""]
  pc[, PTM_Features := ""]
  pc[, TileLabel := GeneLabel]
  pc[]
}






make_ptm_all_sample_plot_table <- function(wq_ptm, sample_meta) {
  dt <- as.data.table(wq_ptm)
  if (nrow(dt) == 0) return(data.table())

  dt <- filter_reportable_ptms(dt)
  if (nrow(dt) == 0) return(data.table())

  dt <- standardize_protein_gene_columns(dt)
  dt[, GeneLabel := get_gene_display_name(.SD), .SDcols = names(dt)]

  # Ensure modsite labels exist at row level.
  labels <- character(nrow(dt))
  for (ii in seq_len(nrow(dt))) {
    labels[ii] <- make_modsite_labels_for_row(dt[["FeatureID"]][ii], dt[ii])
  }
  dt[, ModSiteLabel := labels]
  dt <- dt[nzchar(ModSiteLabel) & !grepl("^UNMAPPED", ModSiteLabel)]
  dt <- dt[!grepl("UniMod:(2|23|30)([^0-9]|$)", ModSiteLabel, perl = TRUE)]
  dt <- dt[!is_peptide_like_label(GeneLabel)]
  if (nrow(dt) == 0) return(data.table())

  required_n <- ptm_required_samples_for_plot
  if (is.null(required_n) || !is.finite(required_n)) {
    required_n <- uniqueN(sample_meta$SampleID)
  }

  out <- dt[, list(
    SamplesPassing = uniqueN(SampleID),
    SamplesPassingList = paste(natural_sort_key(unique(paste0("S", SampleID))), collapse = ";"),
    MedianInterDayCV = median(CV, na.rm = TRUE),
    MeanInterDayCV = mean(CV, na.rm = TRUE),
    Features = paste(sort(unique(FeatureID)), collapse = " | ")
  ), by = c("Search", "GeneLabel", "ModSiteLabel", "PTM")]

  out <- out[SamplesPassing >= required_n]
  out[, RequiredSamples := required_n]
  setorder(out, Search, GeneLabel, ModSiteLabel)
  out[]
}



is_peptide_like_label <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  grepl("^[ACDEFGHIKLMNPQRSTVWY]{6,}$", x) |
    grepl("UniMod", x, ignore.case = TRUE) |
    grepl("\\(|\\)|\\[|\\]", x) |
    grepl("^[A-Z]{8,}[0-9]*$", x)
}


add_ptm_search_uniqueness <- function(ptm_site_dt) {
  dt <- as.data.table(ptm_site_dt)
  if (nrow(dt) == 0) return(dt)

  # Unique means this gene + full-protein mod-site label is observed in exactly one
  # non-All-Searches strategy. All Searches is ignored for uniqueness calls.
  no_all <- dt[Search != ALL_SEARCHES_LABEL]
  if (nrow(no_all) == 0) {
    dt[, UniqueToSearch := FALSE]
    dt[, UniqueSearch := NA_character_]
    return(dt[])
  }

  u <- no_all[, list(
    N_Searches = uniqueN(Search),
    SearchList = paste(sort(unique(Search)), collapse = ";")
  ), by = c("GeneLabel", "ModSiteLabel")]
  u[, UniqueToSearch := N_Searches == 1L]
  u[, UniqueSearch := ifelse(UniqueToSearch, SearchList, NA_character_)]

  dt <- merge(dt, u[, .(GeneLabel, ModSiteLabel, UniqueToSearch, UniqueSearch)],
              by = c("GeneLabel", "ModSiteLabel"), all.x = TRUE)
  dt[is.na(UniqueToSearch), UniqueToSearch := FALSE]
  dt[]
}



extract_modsite_position <- function(x) {
  x <- as.character(x)
  # First number after a residue letter, e.g. S123 U:21 -> 123
  m <- regexpr("[A-Z]([0-9]+)", x, perl = TRUE)
  out <- rep(NA_integer_, length(x))
  hit <- m > 0
  if (any(hit)) {
    token <- regmatches(x, m)
    out[hit] <- suppressWarnings(as.integer(gsub("^[A-Z]", "", token)))
  }
  out
}


plot_ptm_sites_all_samples_by_protein <- function(ptm_site_dt, out_png) {
  dt_all <- as.data.table(ptm_site_dt)
  if (nrow(dt_all) == 0) return(invisible(NULL))

  # Exclude entries without true gene annotation and peptide-only labels.
  dt_all <- dt_all[!is.na(GeneLabel) & nzchar(GeneLabel) & GeneLabel != "NA"]
  dt_all <- dt_all[!grepl("^no_gene$|^no gene$|^UNMAPPED|^pep[0-9]", GeneLabel, ignore.case = TRUE)]
  dt_all <- dt_all[!is_peptide_like_label(GeneLabel)]
  dt_all <- dt_all[!grepl("^pep[0-9]|UNMAPPED", ModSiteLabel, ignore.case = TRUE)]
  dt_all <- dt_all[!is.na(ModSiteLabel) & nzchar(ModSiteLabel)]
  dt_all <- dt_all[!grepl("UniMod:(2|23|30)([^0-9]|$)", ModSiteLabel, perl = TRUE)]
  dt_all <- dt_all[!grepl("UniMod:(2|23|30)([^0-9]|$)", PTM, perl = TRUE)]

  # Exclude compact-plot ambiguous multi-site labels from the plot, but they remain in CSV.
  dt_all <- dt_all[!grepl(";", ModSiteLabel, fixed = TRUE)]

  if (nrow(dt_all) == 0) return(invisible(NULL))

  search_levels <- get_search_levels(dt_all$Search)

  for (s in search_levels) {
    dt <- dt_all[Search == s]
    if (nrow(dt) == 0) next

    # Protein on X-axis, PTM cells stacked above each protein.
    gene_order <- dt[, .N, by = GeneLabel][order(-N, GeneLabel), GeneLabel]
    dt[, GeneLabel := factor(GeneLabel, levels = gene_order)]
    dt[, ModSitePosition := extract_modsite_position(ModSiteLabel)]
    dt[is.na(ModSitePosition), ModSitePosition := .Machine$integer.max]
    dt <- dt[order(GeneLabel, ModSitePosition, ModSiteLabel)]
    dt[, SiteStack := seq_len(.N), by = GeneLabel]

    # Red border for sites unique to this search result.
    dt[, BorderColor := ifelse(UniqueToSearch %in% TRUE & UniqueSearch == s, "red3", "white")]
    dt[, BorderWidth := ifelse(UniqueToSearch %in% TRUE & UniqueSearch == s, 1.25, 0.20)]

    label_size <- ifelse(length(gene_order) > 80, 1.9,
                         ifelse(length(gene_order) > 45, 2.2,
                                ifelse(length(gene_order) > 25, 2.6, 3.1)))

    p <- ggplot(dt, aes(x = GeneLabel, y = SiteStack, fill = MedianInterDayCV)) +
      geom_tile(aes(color = BorderColor, linewidth = BorderWidth), width = 0.88, height = 0.90) +
      geom_text(aes(label = ModSiteLabel), size = label_size, color = "black", fontface = "bold", lineheight = 0.78) +
      scale_color_identity() +
      scale_linewidth_identity() +
      scale_y_continuous(expand = expansion(mult = c(0.005, 0.025)), breaks = NULL) +
      scale_fill_gradientn(colors = coverage_palette(), name = "Median
inter-day CV") +
      labs(
        title = paste0("Well-quantified PTM sites present in all ", unique(dt$RequiredSamples), " samples"),
        subtitle = paste0("Search strategy: ", s, "   |   Stacked N-to-C by residue position; thick red border = unique to this search"),
        x = "Gene / protein",
        y = "PTM sites stacked per protein"
      ) +
      pub_theme(base_size = 16) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()
      )

    out_s <- sub("\\.png$", paste0("_", sanitize_filename(s), ".png"), out_png, ignore.case = TRUE)
    plot_width <- max(14, min(48, length(gene_order) * 0.95))
    plot_height <- max(6.5, min(30, max(dt$SiteStack, na.rm = TRUE) * 0.55))
    save_pub(p, out_s, width = plot_width, height = plot_height)
  }

  invisible(NULL)
}





# ==============================
# PROTEIN SAMPLE COVERAGE OUTPUTS
# ==============================
coverage_palette <- function() {
  # Light presentation-friendly gradient so black text remains readable.
  c("#E8F5E9", "#C8E6C9", "#A5D6A7", "#FFF59D", "#FFE082", "#FFCC80")
}

make_protein_sample_coverage <- function(wq_inter_protein, sample_meta) {
  dt <- as.data.table(wq_inter_protein)
  if (nrow(dt) == 0) return(data.table())

  # One row per Search / Protein / SampleID that meets:
  #   CompleteDays TRUE, inter-day CV below threshold, and all replicate sets present.
  pass <- unique(dt[, c(
    "Search", "Level", "FeatureID", "SampleID", "CV", "MeanIntensity",
    "MeasuredDays", "ExpectedDays", "CompleteDays",
    "Protein", "Gene"
  ), with = FALSE])

  total_samples <- uniqueN(sample_meta$SampleID)

  cov <- pass[, list(
    SamplesPassing = uniqueN(SampleID),
    TotalSamplesInExperiment = total_samples,
    PercentSamplesPassing = 100 * uniqueN(SampleID) / total_samples,
    SamplesPassingList = paste(natural_sort_key(unique(paste0("S", SampleID))), collapse = ";"),
    MedianInterDayCV = median(CV, na.rm = TRUE),
    MeanInterDayCV = mean(CV, na.rm = TRUE),
    MinInterDayCV = min(CV, na.rm = TRUE),
    MaxInterDayCV = max(CV, na.rm = TRUE),
    MedianIntensity = median(MeanIntensity, na.rm = TRUE)
  ), by = c("Search", "Level", "FeatureID", "Protein", "Gene")]

  setorder(cov, Search, -SamplesPassing, MedianInterDayCV, FeatureID)
  cov[]
}

make_protein_sample_matrix <- function(wq_inter_protein, sample_meta) {
  dt <- as.data.table(wq_inter_protein)
  if (nrow(dt) == 0) return(data.table())

  all_samples <- unique(sample_meta[, c("SampleID"), with = FALSE])
  all_samples[, SampleNum := suppressWarnings(as.integer(SampleID))]
  all_samples[is.na(SampleNum), SampleNum := .Machine$integer.max]
  all_samples <- all_samples[order(SampleNum, SampleID)]

  proteins <- unique(dt[, c("Search", "FeatureID", "Protein", "Gene"), with = FALSE])
  proteins[, tmp_join_key := 1]
  all_samples[, tmp_join_key := 1]
  grid <- merge(
    proteins,
    all_samples[, .(SampleID, tmp_join_key)],
    by = "tmp_join_key",
    allow.cartesian = TRUE
  )
  grid[, tmp_join_key := NULL]
  proteins[, tmp_join_key := NULL]
  all_samples[, tmp_join_key := NULL]

  pass <- unique(dt[, c("Search", "FeatureID", "SampleID", "CV"), with = FALSE])
  pass[, Pass := TRUE]

  mat <- merge(grid, pass, by = c("Search", "FeatureID", "SampleID"), all.x = TRUE)
  mat[is.na(Pass), Pass := FALSE]
  mat[]
}

plot_protein_coverage_distribution <- function(protein_cov, out_png) {
  dt <- as.data.table(protein_cov)
  if (nrow(dt) == 0) return(invisible(NULL))

  dt[, Search := factor(Search, levels = unique(Search))]

  p <- ggplot(dt, aes(x = Search, y = SamplesPassing, fill = Search)) +
    geom_violin(trim = TRUE, scale = "width", color = "black", linewidth = 0.25, alpha = 0.85) +
    geom_boxplot(width = 0.13, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.25) +
    stat_summary(fun = median, geom = "point", shape = 23, size = 2.6, fill = "white", color = "black") +
    scale_fill_manual(values = rep(coverage_palette(), length.out = length(unique(dt$Search)))) +
    labs(
      title = "Protein reproducibility: samples passing inter-day CV <20%",
      x = "Search strategy",
      y = "Number of samples where protein passes criteria"
    ) +
    pub_theme(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

  save_pub(p, out_png, width = max(9, length(unique(dt$Search)) * 1.1), height = 6.5)
  invisible(p)
}

plot_top_protein_coverage_bar <- function(protein_cov, out_png, top_n = Inf, rows_per_page = 80) {
  dt <- as.data.table(protein_cov)
  if (nrow(dt) == 0) return(invisible(NULL))

  if (ALL_SEARCHES_LABEL %in% dt$Search) {
    dt_plot <- dt[Search == ALL_SEARCHES_LABEL]
  } else {
    dt_plot <- dt
  }

  dt_plot <- dt_plot[order(-SamplesPassing, MedianInterDayCV)]
  if (is.finite(top_n)) dt_plot <- head(dt_plot, top_n)

  if (!("GeneLabel" %in% names(dt_plot))) {
    dt_plot[, GeneLabel := get_gene_display_name(.SD), .SDcols = names(dt_plot)]
  }
  dt_plot[, TileLabel := GeneLabel]

  dt_plot[, SampleCount := SamplesPassing]
  dt_plot[, RankInCount := seq_len(.N), by = SampleCount]

  total_samples <- max(dt_plot$TotalSamplesInExperiment, na.rm = TRUE)
  max_rank <- max(dt_plot$RankInCount, na.rm = TRUE)
  x_breaks <- seq_len(total_samples)

  label_size <- ifelse(max_rank > 90, 2.8,
                       ifelse(max_rank > 55, 3.2,
                              ifelse(max_rank > 30, 3.8, 4.5)))

  p <- ggplot(dt_plot, aes(x = SampleCount, y = RankInCount, fill = MedianInterDayCV)) +
    geom_tile(color = NA, linewidth = 0, width = 1.0, height = 0.96) +
    geom_text(aes(label = TileLabel), size = label_size, color = "black", fontface = "bold", lineheight = 0.90) +
    scale_x_continuous(breaks = x_breaks, limits = c(0.5, total_samples + 0.5), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_gradientn(colors = coverage_palette(), name = "Median
inter-day CV") +
    labs(
      title = "Reproducible proteins by sample coverage",
      subtitle = "Each tile = one passing protein; text = gene symbol",
      x = "Number of samples passing CV <20% + complete replicate detection",
      y = "Proteins within sample-count category"
    ) +
    pub_theme(base_size = 16) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey82", linewidth = 0.30),
      panel.grid.major.y = element_blank()
    )

  plot_width <- max(9, min(16, total_samples * 1.0))
  plot_height <- max(5.2, min(36, max_rank * 0.28))
  save_pub(p, out_png, width = plot_width, height = plot_height)
  invisible(p)
}

plot_protein_sample_pass_heatmap <- function(protein_matrix, protein_cov, out_png, top_n = 50) {
  mat <- as.data.table(protein_matrix)
  cov <- as.data.table(protein_cov)
  if (nrow(mat) == 0 || nrow(cov) == 0) return(invisible(NULL))

  # Use All Searches if available, otherwise the first search.
  search_use <- if (ALL_SEARCHES_LABEL %in% cov$Search) ALL_SEARCHES_LABEL else cov$Search[1]

  top <- cov[Search == search_use][order(-SamplesPassing, MedianInterDayCV)]
  top <- head(top, top_n)
  dt <- mat[Search == search_use & FeatureID %in% top$FeatureID]
  if (nrow(dt) == 0) return(invisible(NULL))

  dt <- merge(dt, top[, .(FeatureID, SamplesPassing, MedianInterDayCV)], by = "FeatureID", all.x = TRUE)
  dt[, ProteinLabel := get_feature_display_name(.SD), .SDcols = names(dt)]
  dt[, ProteinLabel := factor(ProteinLabel, levels = rev(unique(get_feature_display_name(top))))]

  sample_order <- unique(dt[, .(SampleID)])
  sample_order[, SampleNum := suppressWarnings(as.integer(SampleID))]
  sample_order[is.na(SampleNum), SampleNum := .Machine$integer.max]
  sample_levels <- paste0("S", sample_order[order(SampleNum, SampleID), SampleID])
  dt[, SampleLabel := factor(paste0("S", SampleID), levels = sample_levels)]

  dt[, PassLabel := ifelse(Pass, "Pass", "Fail / not detected")]

  p <- ggplot(dt, aes(x = SampleLabel, y = ProteinLabel, fill = PassLabel)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_manual(values = c("Pass" = "#22A884", "Fail / not detected" = "#FDE725")) +
    labs(
      title = paste0("Protein pass/fail matrix: ", search_use),
      subtitle = "Pass = inter-day CV <20% with complete replicate detection",
      x = "Sample",
      y = NULL,
      fill = NULL
    ) +
    pub_theme(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 7))

  save_pub(p, out_png, width = max(9, length(sample_levels) * 0.45), height = max(7, min(top_n, nrow(top)) * 0.18))
  invisible(p)
}

make_ptm_sample_coverage <- function(wq_ptm, sample_meta) {
  dt <- as.data.table(wq_ptm)
  if (nrow(dt) == 0) return(data.table())

  total_samples <- uniqueN(sample_meta$SampleID)

  cov <- dt[, list(
    SamplesPassing = uniqueN(SampleID),
    TotalSamplesInExperiment = total_samples,
    PercentSamplesPassing = 100 * uniqueN(SampleID) / total_samples,
    SamplesPassingList = paste(natural_sort_key(unique(paste0("S", SampleID))), collapse = ";"),
    MedianInterDayCV = median(CV, na.rm = TRUE),
    MeanInterDayCV = mean(CV, na.rm = TRUE),
    Proteins = uniqueN(Protein, na.rm = TRUE),
    Genes = uniqueN(Gene, na.rm = TRUE)
  ), by = c("Search", "FeatureID", "PTM", "Protein", "Gene")]

  setorder(cov, Search, -SamplesPassing, MedianInterDayCV, PTM, FeatureID)
  cov[]
}

plot_ptm_coverage_bar <- function(ptm_cov, out_png, top_n = Inf, rows_per_page = 80) {
  # PTMs are now integrated into the protein sample-coverage plot.
  invisible(NULL)
}







# ==============================
# WELL-QUANTIFIED PROTEIN INTENSITY PLOTS
# ==============================
make_wq_protein_intensity_long <- function(protein_mats, wq_inter_protein, sample_meta) {
  wq <- as.data.table(wq_inter_protein)
  if (nrow(wq) == 0) return(data.table())
  wq <- standardize_protein_gene_columns(wq)
  wq[, GeneLabel := get_gene_display_name(.SD), .SDcols = names(wq)]

  # Add coverage/stat columns here so this function works even before
  # protein_sample_coverage has been created.
  cov_cols <- wq[, list(
    SamplesPassing = uniqueN(SampleID),
    MedianInterDayCV = median(CV, na.rm = TRUE)
  ), by = c("Search", "FeatureID")]

  wq <- merge(wq, cov_cols, by = c("Search", "FeatureID"), all.x = TRUE, suffixes = c("", "_calc"))

  # If columns already existed and merge created *_calc, fill from calculated versions.
  if ("SamplesPassing_calc" %in% names(wq)) {
    if (!("SamplesPassing" %in% names(wq))) wq[, SamplesPassing := SamplesPassing_calc]
    wq[is.na(SamplesPassing), SamplesPassing := SamplesPassing_calc]
    wq[, SamplesPassing_calc := NULL]
  }
  if ("MedianInterDayCV_calc" %in% names(wq)) {
    if (!("MedianInterDayCV" %in% names(wq))) wq[, MedianInterDayCV := MedianInterDayCV_calc]
    wq[is.na(MedianInterDayCV), MedianInterDayCV := MedianInterDayCV_calc]
    wq[, MedianInterDayCV_calc := NULL]
  }

  out <- list()
  k <- 1L

  for (s in names(protein_mats)) {
    if (!s %in% wq$Search) next

    m <- protein_mats[[s]]
    ids <- unique(wq[Search == s, FeatureID])
    ids <- intersect(ids, rownames(m))
    if (length(ids) == 0) next

    mm <- m[ids, , drop = FALSE]
    dt <- as.data.table(mm, keep.rownames = "FeatureID")
    long <- melt(
      dt,
      id.vars = "FeatureID",
      variable.name = "Sample",
      value.name = "Intensity",
      variable.factor = FALSE
    )
    long <- long[is.finite(Intensity) & Intensity > 0]
    if (nrow(long) == 0) next

    long[, Search := s]
    long <- merge(long, sample_meta[, .(Sample, Day, SampleID, Replicate)], by = "Sample", all.x = TRUE)

    wq_meta <- unique(wq[Search == s, .(Search, FeatureID, GeneLabel, SamplesPassing, MedianInterDayCV)])
    long <- merge(long, wq_meta, by = c("Search", "FeatureID"), all.x = TRUE)

    long[, LogIntensity := if (wq_intensity_log_base == 10) log10(Intensity) else log2(Intensity)]
    long[, SampleLabel := paste0("S", SampleID)]
    long[, ReplicateLabel := paste0("D", Day, "_R", Replicate)]

    out[[k]] <- long
    k <- k + 1L
  }

  rbindlist(out, fill = TRUE)
}

plot_wq_protein_intensities_one_search <- function(wq_int_long, search_name, out_png) {
  dt <- as.data.table(wq_int_long)
  dt <- dt[Search == search_name & is.finite(LogIntensity)]
  if (nrow(dt) == 0) return(invisible(NULL))

  # Order samples numerically so S10 follows S9.
  sample_order <- unique(dt[, .(SampleID)])
  sample_order[, SampleNum := suppressWarnings(as.integer(SampleID))]
  sample_order[is.na(SampleNum), SampleNum := .Machine$integer.max]
  sample_order <- sample_order[order(SampleNum, SampleID)]
  sample_levels <- paste0("S", sample_order$SampleID)
  dt[, SampleLabel := factor(paste0("S", SampleID), levels = sample_levels)]

  # Order proteins by global average intensity, highest to lowest.
  prot_order <- dt[, .(
    MeanLogIntensity = mean(LogIntensity, na.rm = TRUE),
    GeneLabel = GeneLabel[which.max(!is.na(GeneLabel) & nzchar(GeneLabel))]
  ), by = FeatureID]
  prot_order[is.na(GeneLabel) | !nzchar(GeneLabel), GeneLabel := FeatureID]
  prot_order <- prot_order[order(-MeanLogIntensity, GeneLabel)]

  if (is.finite(max_wq_proteins_per_fingerprint_plot)) {
    prot_order <- head(prot_order, max_wq_proteins_per_fingerprint_plot)
    dt <- dt[FeatureID %in% prot_order$FeatureID]
  }

  prot_order[, ProteinLabel := make.unique(GeneLabel)]
  prot_order[, ProteinIndex := seq_len(.N)]
  dt <- merge(dt, prot_order[, .(FeatureID, ProteinLabel, ProteinIndex, MeanLogIntensity)], by = "FeatureID", all.x = TRUE)

  # Summarize each protein within each sample using mean ± SD of replicate measurements.
  sum_dt <- dt[, .(
    MeanLogIntensity = mean(LogIntensity, na.rm = TRUE),
    SDLogIntensity = if (.N >= 2) sd(LogIntensity, na.rm = TRUE) else 0,
    NReplicateMeasurements = .N
  ), by = .(Search, SampleLabel, FeatureID, ProteinLabel, ProteinIndex)]

  sum_dt[, ymin := MeanLogIntensity - SDLogIntensity]
  sum_dt[, ymax := MeanLogIntensity + SDLogIntensity]
  sum_dt[, ProteinIndex := as.integer(ProteinIndex)]

  # Export summary table for this plot.
  fwrite(sum_dt, sub("\\.png$", "_mean_sd_summary.csv", out_png))

  # Use color to identify proteins instead of X-axis labels.
  n_prot <- uniqueN(sum_dt$FeatureID)
  protein_cols <- grDevices::hcl.colors(n_prot, palette = "Dark 3")
  names(protein_cols) <- prot_order$ProteinLabel

  p <- ggplot(sum_dt, aes(x = ProteinIndex, y = MeanLogIntensity, color = ProteinLabel)) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.18, linewidth = 0.35, alpha = 0.90) +
    geom_point(size = 1.9, alpha = 0.95) +
    facet_grid(. ~ SampleLabel) +
    scale_color_manual(values = protein_cols, guide = "none") +
    scale_x_continuous(
      breaks = if (isTRUE(show_protein_names_on_fingerprint_xaxis)) prot_order$ProteinIndex else NULL,
      labels = if (isTRUE(show_protein_names_on_fingerprint_xaxis)) prot_order$ProteinLabel else NULL,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      title = paste0("Well-quantified protein fingerprints by sample: ", search_name),
      subtitle = "Each point is the replicate mean for one protein; error bars show ±SD. Proteins are ordered high-to-low by global average intensity.",
      x = "Proteins ordered by average intensity",
      y = if (wq_intensity_log_base == 10) "log10 intensity, mean ± SD" else "log2 intensity, mean ± SD"
    ) +
    pub_theme(base_size = 15) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text = element_text(face = "bold", color = "black", size = 16),
      legend.position = "none",
      panel.spacing.x = unit(0.22, "lines")
    )

  n_samples <- uniqueN(sum_dt$SampleLabel)
  plot_width <- max(15, min(48, n_samples * max(2.2, n_prot * 0.045)))
  plot_height <- 12

  save_pub(p, out_png, width = plot_width, height = plot_height)
  invisible(p)
}

plot_all_wq_protein_intensities <- function(protein_mats, wq_inter_protein, sample_meta, out_dir) {
  plot_dir <- file.path(out_dir, "Well_Quantified_Protein_Sample_Fingerprint_Plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  wq_int_long <- make_wq_protein_intensity_long(protein_mats, wq_inter_protein, sample_meta)
  fwrite(wq_int_long, file.path(plot_dir, "well_quantified_protein_sample_fingerprint_long.csv"))

  if (nrow(wq_int_long) == 0) return(invisible(NULL))

  for (s in unique(wq_int_long$Search)) {
    plot_wq_protein_intensities_one_search(
      wq_int_long,
      s,
      file.path(plot_dir, paste0("WQ_protein_sample_fingerprint_meanSD_", sanitize_filename(s), ".png"))
    )
  }

  invisible(wq_int_long)
}


# ==============================
# MAIN
# ==============================
search_paths <- discover_search_folders(base_dir, pg_file, pr_file)
search_paths <- prompt_search_short_names(search_paths)

cat("\nReading search folders...\n")
protein_mats <- list()
peptide_mats <- list()
protein_anns <- list()
peptide_anns <- list()

for (s in names(search_paths)) {
  cat("  Reading:", s, "\n")
  pg_path <- file.path(search_paths[[s]], pg_file)
  pr_path <- file.path(search_paths[[s]], pr_file)

  pg <- read_protein_matrix(pg_path)
  pr <- read_peptide_matrix(pr_path)

  # Clean sample column names consistently.
  colnames(pg$matrix) <- make.unique(vapply(colnames(pg$matrix), get_basename_noext, character(1)), sep = "__dup")
  colnames(pr$matrix) <- make.unique(vapply(colnames(pr$matrix), get_basename_noext, character(1)), sep = "__dup")

  protein_mats[[s]] <- pg$matrix
  peptide_mats[[s]] <- pr$matrix
  protein_anns[[s]] <- pg$annotation
  peptide_anns[[s]] <- pr$annotation
}

all_samples <- sort(unique(c(unlist(lapply(protein_mats, colnames)), unlist(lapply(peptide_mats, colnames)))))
sample_meta <- build_sample_meta(all_samples)
validate_sample_meta(sample_meta)
fwrite(sample_meta, file.path(out_dir, "sample_metadata_parsed_all_columns.csv"))

# Keep only successfully parsed replicate sample columns for CV analysis.
sample_meta_all <- copy(sample_meta)
sample_meta <- sample_meta[ParseOK %in% TRUE & !is.na(Day) & !is.na(SampleID) & !is.na(Replicate)]
if (nrow(sample_meta) == 0) {
  stop("No valid Day_SampleID_Replicate sample columns remain after parsing. Check sample_name_parse_failures.csv.")
}
all_samples <- sample_meta$Sample
fwrite(sample_meta, file.path(out_dir, "sample_metadata_parsed_used_for_CV.csv"))

cat("\nLoading FASTA for full-protein PTM site mapping...\n")
fasta_path_used <- pick_fasta_file()
fasta_db <- read_fasta_sequences(fasta_path_used)
if (!is.na(fasta_path_used)) {
  cat("FASTA loaded: ", fasta_path_used, "\n", sep = "")
  cat("Protein sequences indexed: ", length(fasta_db$seq_by_accession), "\n", sep = "")
} else {
  cat("No FASTA selected. PTM sites can only use DIA-NN protein-position columns if present; otherwise they will be marked UNMAPPED.\n")
}


# Align samples within levels, using only parsed replicate sample columns.
protein_mats <- lapply(protein_mats, align_matrix_samples, all_samples = all_samples)
peptide_mats <- lapply(peptide_mats, align_matrix_samples, all_samples = all_samples)

# Combined matrix across search strategies.
protein_mats[[ALL_SEARCHES_LABEL]] <- build_all_searches_matrix(protein_mats)
peptide_mats[[ALL_SEARCHES_LABEL]] <- build_all_searches_matrix(peptide_mats)

searches_with_all <- names(protein_mats)
searches_no_all <- setdiff(searches_with_all, ALL_SEARCHES_LABEL)

protein_ann <- combine_annotations(protein_anns)
peptide_ann <- combine_annotations(peptide_anns)
all_ann <- rbindlist(list(protein_ann, peptide_ann), fill = TRUE)

fwrite(protein_ann, file.path(out_dir, "protein_annotation_table.csv"))
fwrite(peptide_ann, file.path(out_dir, "peptide_annotation_table.csv"))

cat("\nCalculating CVs...\n")
intra_list <- list()
inter_list <- list()
idx <- 1

for (s in searches_with_all) {
  cat("  CV:", s, "Protein\n")
  intra_list[[idx]] <- calc_intra_day_cv(protein_mats[[s]], s, "Protein", sample_meta); idx <- idx + 1
  inter_list[[idx]] <- calc_inter_day_cv(protein_mats[[s]], s, "Protein", sample_meta); idx <- idx + 1

  cat("  CV:", s, "Peptide\n")
  intra_list[[idx]] <- calc_intra_day_cv(peptide_mats[[s]], s, "Peptide", sample_meta); idx <- idx + 1
  inter_list[[idx]] <- calc_inter_day_cv(peptide_mats[[s]], s, "Peptide", sample_meta); idx <- idx + 1
}

intra_cv_all <- rbindlist(intra_list, fill = TRUE)
inter_cv_all <- rbindlist(inter_list, fill = TRUE)

intra_cv_all <- add_annotations(intra_cv_all, all_ann)
inter_cv_all <- add_annotations(inter_cv_all, all_ann)

fwrite(intra_cv_all, file.path(out_dir, "feature_level_intra_day_CVs_all_searches.csv"))
fwrite(inter_cv_all, file.path(out_dir, "feature_level_inter_day_CVs_all_searches.csv"))

cat("\nPlotting CV violin distributions...\n")
plot_all_cv_violins(intra_cv_all, inter_cv_all, searches_with_all, out_dir)

cat("\nIdentifying well-quantified features from inter-day CV...\n")

# Proteins: well quantified = detected across all expected days/replicate sets and inter-day CV < 20%.
wq_inter_protein <- identify_well_quantified_interday(
  inter_cv_all[Level == "Protein"],
  threshold = inter_cv_threshold,
  sample_meta = sample_meta
)
wq_inter_protein <- add_annotations(wq_inter_protein, all_ann)
wq_inter_protein <- standardize_protein_gene_columns(wq_inter_protein)
fwrite(wq_inter_protein, file.path(out_dir, "well_quantified_inter_day_PROTEINS_CVlt20_completeDays_sampleCoverage.csv"))

cat("\nPlotting sample-level protein fingerprints for well-quantified proteins...\n")
plot_all_wq_protein_intensities(
  protein_mats = protein_mats,
  wq_inter_protein = wq_inter_protein,
  sample_meta = sample_meta,
  out_dir = out_dir
)



cat("\nPlotting well-quantified protein count heatmap...\n")
plot_wq_count_heatmap(
  wq_inter_protein,
  file.path(out_dir, "well_quantified_inter_day_protein_counts_by_search_sample.png")
)



# PTMs / modified peptides: well quantified = detected across all expected days/replicate sets and inter-day CV < 40%.
wq_inter_peptide_ptm <- identify_well_quantified_interday(
  inter_cv_all[Level == "Peptide"],
  threshold = ptm_inter_cv_threshold,
  sample_meta = sample_meta
)
wq_inter_peptide_ptm <- add_annotations(wq_inter_peptide_ptm, all_ann)
wq_inter_peptide_ptm <- standardize_protein_gene_columns(wq_inter_peptide_ptm)

cat("\nCreating PTM tables for well-quantified modified peptides...\n")
wq_ptm <- make_ptm_table(wq_inter_peptide_ptm)
wq_ptm <- filter_reportable_ptms(wq_ptm)
fwrite(wq_ptm, file.path(out_dir, "well_quantified_inter_day_PTM_peptides_CVlt40_filtered_with_protein_gene_sampleCoverage.csv"))

ptm_all_sample_sites <- make_ptm_all_sample_plot_table(wq_ptm, sample_meta)
ptm_all_sample_sites <- add_ptm_search_uniqueness(ptm_all_sample_sites)
ptm_all_sample_sites_removed_from_plot <- ptm_all_sample_sites[
  is_peptide_like_label(GeneLabel) |
    grepl("UniMod:(2|23|30)([^0-9]|$)", ModSiteLabel, perl = TRUE) |
    grepl("^pep[0-9]|UNMAPPED", ModSiteLabel, ignore.case = TRUE)
]
fwrite(ptm_all_sample_sites_removed_from_plot, file.path(out_dir, "PTM_sites_removed_from_plot_peptide_like_or_excluded_unimod.csv"))

fwrite(ptm_all_sample_sites, file.path(out_dir, "well_quantified_PTM_sites_present_all_8_samples_by_gene.csv"))
plot_ptm_sites_all_samples_by_protein(
  ptm_all_sample_sites,
  file.path(out_dir, "well_quantified_PTM_sites_present_all_8_samples_by_gene_per_search.png")
)

if ("PTM_Sites" %in% names(wq_ptm)) {
  unmapped_ptms <- wq_ptm[grepl("UNMAPPED", PTM_Sites, fixed = TRUE)]
  fwrite(unmapped_ptms, file.path(out_dir, "well_quantified_unmapped_PTM_sites_need_FASTA_or_position_columns.csv"))
}

if (nrow(wq_ptm) > 0) {
  ptm_summary <- wq_ptm[, list(
    WellQuantifiedPTMPeptides = uniqueN(FeatureID),
    Proteins = uniqueN(Protein, na.rm = TRUE),
    Genes = uniqueN(Gene, na.rm = TRUE),
    SamplesDetected = uniqueN(SampleID),
    SamplesDetectedList = paste(natural_sort_key(unique(paste0("S", SampleID))), collapse = ";"),
    MedianInterDayCV = median(CV, na.rm = TRUE)
  ), by = c("Search", "PTM")]
  fwrite(ptm_summary, file.path(out_dir, "well_quantified_PTM_summary_by_search_CVlt40_sampleCoverage.csv"))
}

cat("
Skipping old well_quantified_inter_day_protein_counts_by_search_sample heatmap; protein tile coverage plot replaces it.
")

cat("\nBuilding protein-centric sample coverage outputs...\n")
protein_sample_coverage <- make_protein_sample_coverage(wq_inter_protein, sample_meta)
protein_sample_coverage <- standardize_protein_gene_columns(protein_sample_coverage)
protein_sample_coverage <- force_genes_column(protein_sample_coverage)
protein_sample_coverage[, GeneLabelDiagnostic := get_gene_display_name(.SD), .SDcols = names(protein_sample_coverage)]
protein_sample_coverage <- force_genes_column(protein_sample_coverage)
diag_cols <- intersect(c("Search", "FeatureID", "Genes", "Gene", "Protein", "GeneLabelDiagnostic"), names(protein_sample_coverage))
fwrite(
  protein_sample_coverage[, ..diag_cols],
  file.path(out_dir, "gene_label_diagnostics.csv")
)
protein_sample_coverage <- attach_ptms_to_protein_coverage(protein_sample_coverage, wq_ptm)
fwrite(protein_sample_coverage, file.path(out_dir, "protein_sample_coverage_interDayCVlt20_completeReplicates_withPTMs.csv"))

protein_sample_matrix <- make_protein_sample_matrix(wq_inter_protein, sample_meta)
fwrite(protein_sample_matrix, file.path(out_dir, "protein_sample_pass_matrix_interDayCVlt20_completeReplicates.csv"))

cat("
Skipping protein_samples_passing_distribution_by_search plot.
")
plot_top_protein_coverage_bar(
  protein_sample_coverage,
  file.path(out_dir, "reproducible_proteins_by_sample_count_withPTMs.png"),
  top_n = Inf,
  rows_per_page = 80
)

cat("
Skipping top_proteins_sample_pass_fail_heatmap plot.
")
cat("\nBuilding PTM sample coverage outputs...\n")
ptm_sample_coverage <- make_ptm_sample_coverage(wq_ptm, sample_meta)
fwrite(ptm_sample_coverage, file.path(out_dir, "ptm_sample_coverage_interDayCVlt40_filtered_completeReplicates.csv"))

plot_ptm_coverage_bar(
  ptm_sample_coverage,
  file.path(out_dir, "all_reproducible_PTMs_gene_tiles_by_sample_count_filtered.png"),
  top_n = Inf,
  rows_per_page = 80
)

# Compact summary tables.
summary_intra <- intra_cv_all[, list(
  FeaturesWithCV = sum(!is.na(CV)),
  MedianCV = median(CV, na.rm = TRUE),
  MeanCV = mean(CV, na.rm = TRUE),
  FeaturesCompleteReplicates = sum(CompleteReplicates %in% TRUE, na.rm = TRUE),
  FeaturesLowCVComplete = sum((!is.na(CV)) & CV <= intra_cv_threshold & CompleteReplicates %in% TRUE, na.rm = TRUE)
), by = c("Search", "Level", "SampleID", "Day")]
fwrite(summary_intra, file.path(out_dir, "summary_intra_day_CV_by_search_level_sample_day.csv"))

summary_inter <- inter_cv_all[, list(
  FeaturesWithCV = sum(!is.na(CV)),
  MedianCV = median(CV, na.rm = TRUE),
  MeanCV = mean(CV, na.rm = TRUE),
  FeaturesCompleteDays = sum(CompleteDays %in% TRUE, na.rm = TRUE),
  FeaturesLowCVComplete = sum((!is.na(CV)) & CV <= inter_cv_threshold & CompleteDays %in% TRUE, na.rm = TRUE)
), by = c("Search", "Level", "SampleID")]
fwrite(summary_inter, file.path(out_dir, "summary_inter_day_CV_by_search_level_sample.csv"))

cat("\nKey outputs:\n")
cat("  - CV violin plots: ", file.path(out_dir, "CV_Violin_Plots"), "\n", sep = "")
cat("  - Intra-day CV table: feature_level_intra_day_CVs_all_searches.csv\n")
cat("  - Inter-day CV table: feature_level_inter_day_CVs_all_searches.csv\n")
cat("  - Protein sample coverage: protein_sample_coverage_interDayCVlt20_completeReplicates_withPTMs.csv\n")
cat("  - Well-quantified peptides: well_quantified_inter_day_PTM_peptides_CVlt40_filtered_with_protein_gene_sampleCoverage.csv\n")
cat("  - PTM sample coverage: ptm_sample_coverage_interDayCVlt40_filtered_completeReplicates.csv\n")
cat("\nDone.\n")

