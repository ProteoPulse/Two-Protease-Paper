# --------------------------------------------------------------------------------------------------
# Two-Protease-Paper | Combine DIA-NN spectral libraries with RT / ion-mobility alignment
# Manuscript: Methods, 'Data Analysis' (DDA and 1 m/z DIA libraries merged per protease).
# Author: Simion Kreimer.
# Run: interactively, combine_diann_libraries() (Windows file pickers), or pass paths directly:
#      combine_diann_libraries(files = c(ref.tsv, other.tsv), output_file = 'combined.tsv').
# Packages: data.table, stringr, tools.
# Released as used for the manuscript: only file paths, line endings and this header were changed;
# the analysis code is untouched. Paths come from environment variables, see ../../README.md.
# --------------------------------------------------------------------------------------------------

# =========================
# Combine DIA-NN Libraries (TSV) with RT/IM alignment
# ONE-BY-ONE file import (keeps prompting until Cancel)
# FIXED: IM_adj column always exists (NA if IM alignment skipped)
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(tools)
})

# ---------------------------
# File chooser helpers (one-by-one)
# ---------------------------
choose_files_one_by_one <- function(prompt_first = "Select a DIA-NN library TSV (Cancel to finish)") {
  cat(prompt_first, "\n")
  files <- character()
  
  repeat {
    f <- tryCatch({
      if (.Platform$OS.type == "windows") {
        choose.files()
      } else {
        file.choose()
      }
    }, error = function(e) "")
    
    if (length(f) == 0 || f == "" || is.na(f)) break
    
    files <- c(files, f)
    cat("Added:", f, "\n")
    cat("Select another file (or Cancel to finish)...\n")
  }
  
  if (length(files) < 2) stop("Please select at least 2 libraries to combine.")
  unique(files)
}

choose_save_path <- function(default_name = "combined_DIANN.tsv", prompt = "Choose output file path") {
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
# Required columns checker
# ---------------------------
require_cols <- function(dt, cols, label = "library") {
  missing <- setdiff(cols, names(dt))
  if (length(missing)) {
    stop(label, " is missing required columns: ", paste(missing, collapse = ", "))
  }
}

# ---------------------------
# Summarize precursors for alignment
# Use ModifiedPeptideSequence + PrecursorCharge as key
# Take median RT/IM across transitions for each precursor
# ---------------------------
summarize_precursors <- function(dt) {
  dt[, .(
    RT = median(NormalizedRetentionTime, na.rm = TRUE),
    IM = if (all(is.na(PrecursorIonMobility))) NA_real_ else median(PrecursorIonMobility, na.rm = TRUE),
    nTransitions = .N
  ), by = .(ModifiedPeptideSequence, PrecursorCharge)]
}

# ---------------------------
# Fit alignment model
# x = library RT/IM, y = reference RT/IM
# Uses LOESS when enough overlap; otherwise linear regression.
# Returns a list with type + predict() function.
# ---------------------------
fit_alignment <- function(x, y, loess_min = 50, span = 0.75) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 5) return(list(type = "none", n = n, predict = function(v) v))
  
  if (n >= loess_min && length(unique(x)) >= 10) {
    mod <- try(loess(y ~ x, span = span, control = loess.control(surface = "direct")), silent = TRUE)
    if (!inherits(mod, "try-error")) {
      return(list(
        type = sprintf("loess(span=%.2f)", span),
        n = n,
        predict = function(v) as.numeric(predict(mod, newdata = data.frame(x = v)))
      ))
    }
  }
  
  mod <- lm(y ~ x)
  co <- coef(mod)
  list(
    type = sprintf("linear(y=%.6f + %.6f*x)", co[1], co[2]),
    n = n,
    predict = function(v) as.numeric(co[1] + co[2] * v)
  )
}

# ---------------------------
# Metrics helper
# ---------------------------
alignment_metrics <- function(before_x, after_x, y) {
  ok1 <- is.finite(before_x) & is.finite(y)
  ok2 <- is.finite(after_x) & is.finite(y)
  
  rmse <- function(a, b) sqrt(mean((a - b)^2))
  mae  <- function(a, b) median(abs(a - b))
  
  list(
    rmse_before = if (sum(ok1) >= 3) rmse(before_x[ok1], y[ok1]) else NA_real_,
    rmse_after  = if (sum(ok2) >= 3) rmse(after_x[ok2],  y[ok2]) else NA_real_,
    mae_before  = if (sum(ok1) >= 3) mae(before_x[ok1],  y[ok1]) else NA_real_,
    mae_after   = if (sum(ok2) >= 3) mae(after_x[ok2],   y[ok2]) else NA_real_,
    n_used      = sum(ok1 & ok2)
  )
}

# ---------------------------
# Deduplicate combined library transitions
# Keep row with highest LibraryIntensity per transition key
# ---------------------------
deduplicate_transitions <- function(dt) {
  key_cols <- c(
    "ModifiedPeptideSequence", "PrecursorCharge",
    "PrecursorMz", "ProductMz",
    "FragmentType", "FragmentCharge", "FragmentSeriesNumber", "FragmentLossType"
  )
  if (!"FragmentLossType" %in% names(dt)) dt[, FragmentLossType := "noloss"]
  require_cols(dt, key_cols, label = "combined library")
  
  setorderv(dt, c(key_cols, "LibraryIntensity"), c(rep(1, length(key_cols)), -1))
  dt[, .SD[1], by = key_cols]
}

# ---------------------------
# Main combiner
# ---------------------------
combine_diann_libraries <- function(files = NULL, output_file = NULL) {
  if (is.null(files)) files <- choose_files_one_by_one()
  if (length(files) < 2) stop("Select at least 2 libraries to combine.")
  if (is.null(output_file)) output_file <- choose_save_path("combined_DIANN.tsv", "Choose output file")
  
  cat("\nLibraries selected (in order; #1 is reference):\n")
  for (i in seq_along(files)) cat(sprintf("  %d) %s\n", i, files[i]))
  
  libs <- lapply(files, fread)
  
  required <- c(
    "PrecursorMz","ProductMz","Annotation","PeptideSequence","ModifiedPeptideSequence",
    "PrecursorCharge","LibraryIntensity","NormalizedRetentionTime",
    "FragmentType","FragmentCharge","FragmentSeriesNumber","FragmentLossType"
  )
  
  # Ensure IM column exists everywhere (optional)
  libs <- lapply(libs, function(dt) {
    if (!"PrecursorIonMobility" %in% names(dt)) dt[, PrecursorIonMobility := NA_real_]
    dt
  })
  
  for (i in seq_along(libs)) {
    require_cols(libs[[i]], required, label = paste0("library ", i, " (", basename(files[i]), ")"))
  }
  
  # Reference is first selected
  ref_idx <- 1
  ref_name <- basename(files[ref_idx])
  cat("\nReference library:", ref_name, "\n")
  
  ref <- libs[[ref_idx]]
  ref_prec <- summarize_precursors(ref)
  setkey(ref_prec, ModifiedPeptideSequence, PrecursorCharge)
  
  report_list <- list()
  overlap_all <- list()
  
  for (i in seq_along(libs)) {
    if (i == ref_idx) next
    lib <- libs[[i]]
    lib_name <- basename(files[i])
    
    lib_prec <- summarize_precursors(lib)
    setkey(lib_prec, ModifiedPeptideSequence, PrecursorCharge)
    
    ov <- lib_prec[ref_prec, nomatch = 0]
    
    setnames(ov,
             old = c("RT","IM","nTransitions","i.RT","i.IM","i.nTransitions"),
             new = c("RT_lib","IM_lib","nTrans_lib","RT_ref","IM_ref","nTrans_ref"))
    
    n_overlap <- nrow(ov)
    cat("\n--- Aligning:", lib_name, "->", ref_name, "\n")
    cat("Overlapping precursors:", n_overlap, "\n")
    
    # RT model
    rt_model <- fit_alignment(ov$RT_lib, ov$RT_ref, loess_min = 50, span = 0.75)
    ov[, RT_adj := rt_model$predict(RT_lib)]
    rt_metrics <- alignment_metrics(ov$RT_lib, ov$RT_adj, ov$RT_ref)
    
    cat("RT model:", rt_model$type, " (n=", rt_model$n, ")\n", sep = "")
    cat(sprintf("RT RMSE before/after: %.4f -> %.4f | median abs err: %.4f -> %.4f\n",
                rt_metrics$rmse_before, rt_metrics$rmse_after, rt_metrics$mae_before, rt_metrics$mae_after))
    
    # IM model (optional)
    im_used <- FALSE
    im_model <- list(type = "none", n = 0, predict = function(v) v)
    im_metrics <- list(rmse_before = NA_real_, rmse_after = NA_real_, mae_before = NA_real_, mae_after = NA_real_, n_used = 0)
    
    if (!all(is.na(ov$IM_lib)) && !all(is.na(ov$IM_ref))) {
      if (sum(is.finite(ov$IM_lib) & is.finite(ov$IM_ref)) >= 10) {
        im_model <- fit_alignment(ov$IM_lib, ov$IM_ref, loess_min = 50, span = 0.75)
        ov[, IM_adj := im_model$predict(IM_lib)]
        im_metrics <- alignment_metrics(ov$IM_lib, ov$IM_adj, ov$IM_ref)
        im_used <- TRUE
        
        cat("IM model:", im_model$type, " (n=", im_model$n, ")\n", sep = "")
        cat(sprintf("IM RMSE before/after: %.4f -> %.4f | median abs err: %.4f -> %.4f\n",
                    im_metrics$rmse_before, im_metrics$rmse_after, im_metrics$mae_before, im_metrics$mae_after))
      } else {
        cat("IM alignment skipped: not enough overlapping non-NA IM values.\n")
      }
    } else {
      cat("IM alignment skipped: IM missing/NA in one of the libraries.\n")
    }
    
    # FIX: ensure IM_adj exists so reporting never crashes
    if (!im_used) {
      ov[, IM_adj := NA_real_]
    }
    
    # Apply adjustments to all rows in this library
    lib[, NormalizedRetentionTime := rt_model$predict(NormalizedRetentionTime)]
    if (im_used) {
      lib[, PrecursorIonMobility := ifelse(is.na(PrecursorIonMobility), NA_real_, im_model$predict(PrecursorIonMobility))]
    }
    libs[[i]] <- lib
    
    # Per-library report (always includes IM_adj now)
    ov_report <- ov[, .(
      ModifiedPeptideSequence, PrecursorCharge,
      RT_ref, RT_lib, RT_adj,
      IM_ref, IM_lib, IM_adj
    )]
    
    rep_path <- file.path(dirname(output_file),
                          paste0(file_path_sans_ext(basename(files[i])), "_alignment_report.tsv"))
    fwrite(ov_report, rep_path, sep = "\t", quote = FALSE, na = "NA")
    cat("Wrote overlap/alignment report:", rep_path, "\n")
    
    report_list[[lib_name]] <- list(
      library = lib_name,
      reference = ref_name,
      n_overlap = n_overlap,
      rt_model = rt_model$type,
      rt_rmse_before = rt_metrics$rmse_before,
      rt_rmse_after  = rt_metrics$rmse_after,
      rt_mae_before  = rt_metrics$mae_before,
      rt_mae_after   = rt_metrics$mae_after,
      im_used = im_used,
      im_model = if (im_used) im_model$type else "skipped",
      im_rmse_before = im_metrics$rmse_before,
      im_rmse_after  = im_metrics$rmse_after,
      im_mae_before  = im_metrics$mae_before,
      im_mae_after   = im_metrics$mae_after
    )
    
    overlap_all[[lib_name]] <- unique(ov$ModifiedPeptideSequence)
  }
  
  # Combine and deduplicate
  combined <- rbindlist(libs, fill = TRUE)
  cat("\nCombined rows (before dedupe):", nrow(combined), "\n")
  combined_dedup <- deduplicate_transitions(combined)
  cat("Combined rows (after dedupe):", nrow(combined_dedup), "\n")
  
  fwrite(combined_dedup, output_file, sep = "\t", quote = FALSE, na = "NA")
  cat("\n✔ Combined DIA-NN library written to:\n", output_file, "\n")
  
  # Write summary report
  if (length(report_list) > 0) {
    summary_dt <- rbindlist(lapply(report_list, as.data.table), fill = TRUE)
    sum_path <- file.path(dirname(output_file), "combine_alignment_summary.tsv")
    fwrite(summary_dt, sum_path, sep = "\t", quote = FALSE, na = "NA")
    cat("✔ Alignment summary written to:\n", sum_path, "\n")
  }
  
  # Write overlapping peptides list
  ov_path <- file.path(dirname(output_file), "overlapping_modified_peptides.tsv")
  ov_dt <- rbindlist(lapply(names(overlap_all), function(nm) {
    data.table(Library = nm, ModifiedPeptideSequence = overlap_all[[nm]])
  }), fill = TRUE)
  fwrite(ov_dt, ov_path, sep = "\t", quote = FALSE, na = "NA")
  cat("✔ Overlapping ModifiedPeptideSequence list written to:\n", ov_path, "\n")
  
  invisible(output_file)
}

# ---------------------------
# HOW TO RUN
# ---------------------------
# combine_diann_libraries()
#
# You will be prompted to choose libraries ONE AT A TIME.
# Click Cancel when you're done selecting.
# The FIRST selected library is the reference.
#
# Outputs next to the combined file:
#  - *_alignment_report.tsv for each non-reference library
#  - combine_alignment_summary.tsv
#  - overlapping_modified_peptides.tsv
