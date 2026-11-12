# utils_ids.R — protein-group canonicalization, run/sample/preparation parsing.
# Conventions re-derived (fresh) from the actual matrices; validated in-pipeline.

suppressPackageStartupMessages({ library(data.table) })

# ---- protein-group accession canonicalization --------------------------------
# Trim whitespace, standardize separator to ';', sort accessions, PRESERVE isoform
# suffixes (e.g. -2) and the original string. Exact canonical set = primary join key.
canonicalize_pg <- function(pg) {
  vapply(pg, function(x) {
    if (is.na(x) || !nzchar(x)) return(NA_character_)
    parts <- unlist(strsplit(x, "[;,]"))
    parts <- trimws(parts)
    parts <- parts[nzchar(parts)]
    paste(sort(unique(parts)), collapse = ";")
  }, character(1), USE.NAMES = FALSE)
}
leading_accession <- function(pg) sub("[;,].*$", "", trimws(pg))
accession_set <- function(pg) sort(unique(trimws(unlist(strsplit(pg, "[;,]")))))
is_single_accession <- function(pg) !grepl("[;,]", pg)

# ---- run-column preparation parsing ------------------------------------------
# Discovery columns look like (backslash Windows paths):
#   ...\20min_Mix_DIA\251117_DIA_Plasma_Mix_BRH000000_20min.raw   -> Mix20
#   ...\30min_Mix_DIA\251117_DIA_Plasma_Mix_BRH000000.raw         -> Mix30
#   ...\DIA_Trypsin\251117_DIA_Plasma_Try_BRH000000.raw           -> Trypsin2 (Nov batch = mixed)
#   ...\DIA_Trypsin\250807_KreimerSimion_EC_BOB_Plasma_016_BRH000000.raw -> Trypsin1 (Aug batch)
#   ...\DIA_Krakatoa\250807_DIA_Plasma_Kra_BRH000000.raw          -> Krakatoa1
folder_token <- function(col) {
  p <- gsub("\\\\", "/", col)
  parts <- strsplit(p, "/", fixed = TRUE)[[1]]
  if (length(parts) < 2) return(NA_character_)
  parts[length(parts) - 1]
}
acquisition_date <- function(col) {
  bn <- basename(gsub("\\\\", "/", col))
  m <- regmatches(bn, regexpr("^[0-9]{6}", bn))
  if (length(m) == 0) NA_character_ else m
}

# Map (folder token, date) -> canonical preparation label using config.
parse_preparation <- function(cols, cfg) {
  ft <- vapply(cols, folder_token, character(1))
  dt <- vapply(cols, acquisition_date, character(1))
  tok <- cfg$run_mapping$prep_tokens
  bd  <- cfg$run_mapping$batch_dates
  out <- rep(NA_character_, length(cols))
  out[ft == tok$Mix20] <- "Mix20"
  out[ft == tok$Mix30] <- "Mix30"
  is_try <- ft == tok$pure_trypsin
  is_kra <- ft == tok$pure_krakatoa
  bdv <- unlist(bd)                                 # named vector: date -> 1/2
  batch_idx <- bdv[as.character(dt)]                # aligned to each run's date
  out[is_try] <- paste0("Trypsin",  batch_idx[is_try])
  out[is_kra] <- paste0("Krakatoa", batch_idx[is_kra])
  out
}

# ---- biological sample id extraction (multi-scheme) --------------------------
# Extract the BRH######## or HMN########[A|B] token regardless of naming scheme;
# A/B suffixes are KEPT DISTINCT (they are distinct biological samples).
parse_sample_id <- function(cols) {
  bn <- basename(gsub("\\\\", "/", cols)); bn <- sub("\\.raw$", "", bn, ignore.case = TRUE)
  sid <- rep(NA_character_, length(bn))
  # primary: explicit BRH/HMN token anywhere in the name
  m <- regmatches(bn, regexpr("(BRH[0-9]+|HMN[0-9]+[AB]?)", bn))
  has <- vapply(regexpr("(BRH[0-9]+|HMN[0-9]+[AB]?)", bn), function(x) x > 0, logical(1))
  sid[has] <- m
  sid
}

# ---- reproducibility run parsing: Day_SampleID_Replicate ---------------------
# e.g. "1_10_1" -> day=1, sample=10, rep=1 (day==plate; confounded).
parse_repro_run <- function(cols) {
  bn <- basename(gsub("\\\\", "/", cols)); bn <- sub("\\.raw$", "", bn, ignore.case = TRUE)
  toks <- lapply(strsplit(bn, "[_\\.-]"), function(z) z[grepl("^[0-9]+$", z)])
  last3 <- lapply(toks, function(z) if (length(z) >= 3) tail(z, 3) else rep(NA_character_, 3))
  m <- do.call(rbind, last3)
  data.table(original_column = cols, day = m[,1], sample = m[,2], replicate = m[,3])
}
