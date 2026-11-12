# qualify_precursors.R — Stage 1: universal precursor qualification.
#
# Applies ONE identical rule to every precursor so protein roll-up (Stage 2) rests only on
# precursors that are actually reproducible. Metrics are LEARNED on the reproducibility runs
# (replicate structure present) and the thresholds are FROZEN before discovery. On discovery data
# (no replicates) we recompute detection_rate and coherence but carry CV/ICC from the frozen
# reproducibility reference (matched by precursor_id).
#
# Per-precursor metrics (over the run universe of one dataset):
#   detection_rate  fraction of runs with a finite intensity
#   median_cv       median across biological samples of the within-sample CV over replicates (linear)
#   icc             ICC(1) of the precursor's per-run log2 intensity grouped by biological sample
#   coherence       median pairwise Pearson of this precursor to its SIBLINGS (same protein+protease)
#   median_qvalue   median DIA-NN q-value across runs
#   is_proteotypic  proteotypic flag
#   n_siblings      number of sibling precursors on the same protein+protease
#
# Depends on: data.table, feature_matrix.R (feature_matrix, median_pairwise_cor), icc.R (icc_oneway).

suppressPackageStartupMessages({ library(data.table) })

# Median within-group CV on the linear scale. values are log2; groups are biological sample ids.
cv_across_reps <- function(log2_values, groups) {
  lin <- 2^log2_values
  d <- data.table(lin = lin, g = groups)[is.finite(lin)]
  if (!nrow(d)) return(NA_real_)
  per <- d[, .(cv = if (.N >= 2 && mean(lin) > 0) stats::sd(lin) / mean(lin) else NA_real_), by = g]
  stats::median(per$cv, na.rm = TRUE)
}

# Per-precursor qualification metrics for one dataset.
#   prec_long : long precursors. Columns: pg_canonical, protease, precursor_id, stripped_sequence,
#               sample_id (biological), run, log2q, and optionally q_value, proteotypic.
#   have_replicates : TRUE for reproducibility data (compute CV + ICC); FALSE for discovery.
precursor_metrics <- function(prec_long, have_replicates = TRUE) {
  d <- prec_long[is.finite(log2q)]
  n_runs <- d[, uniqueN(run)]
  # base per-precursor summaries
  base <- d[, .(
    pg_canonical = pg_canonical[1], protease = protease[1],
    stripped_sequence = stripped_sequence[1],
    n_runs_obs = uniqueN(run),
    detection_rate = uniqueN(run) / n_runs,
    median_qvalue = if ("q_value" %in% names(d)) stats::median(q_value, na.rm = TRUE) else NA_real_,
    is_proteotypic = if ("proteotypic" %in% names(d)) as.logical(max(proteotypic, na.rm = TRUE)) else NA
  ), by = precursor_id]
  # CV + ICC (replicate-based) — only meaningful when replicates exist
  if (have_replicates) {
    cvicc <- d[, .(
      median_cv = cv_across_reps(log2q, sample_id),
      icc = icc_oneway(log2q, sample_id)
    ), by = precursor_id]
    base <- merge(base, cvicc, by = "precursor_id", all.x = TRUE)
  } else {
    base[, `:=`(median_cv = NA_real_, icc = NA_real_)]
  }
  # sibling coherence: per protein x protease, run x precursor matrix; each precursor's median
  # pairwise Pearson to the others of the same protein+protease.
  coh <- vector("list", 0L)
  keys <- unique(d[, .(pg_canonical, protease)])
  for (i in seq_len(nrow(keys))) {
    pg <- keys$pg_canonical[i]; pro <- keys$protease[i]
    sub <- d[pg_canonical == pg & protease == pro]
    prec_ids <- unique(sub$precursor_id); k <- length(prec_ids)
    if (k < 2) { coh[[length(coh)+1L]] <- data.table(precursor_id = prec_ids, coherence = NA_real_, n_siblings = 0L); next }
    m <- feature_matrix(sub, id_col = "run", feature_col = "precursor_id", value_col = "log2q")
    if (is.null(m)) { coh[[length(coh)+1L]] <- data.table(precursor_id = prec_ids, coherence = NA_real_, n_siblings = k-1L); next }
    cc <- suppressWarnings(stats::cor(m, method = "pearson", use = "pairwise.complete.obs"))
    per <- vapply(colnames(m), function(j) stats::median(cc[j, setdiff(colnames(m), j)], na.rm = TRUE), numeric(1))
    coh[[length(coh)+1L]] <- data.table(precursor_id = colnames(m), coherence = as.numeric(per), n_siblings = ncol(m)-1L)
  }
  cohdt <- rbindlist(coh, fill = TRUE)
  out <- merge(base, cohdt, by = "precursor_id", all.x = TRUE)
  out[is.na(n_siblings), n_siblings := 0L]
  out[]
}

# Apply frozen thresholds. Singletons (no siblings) pass the coherence rule with a flag.
#   thr : list(detection_rate_min, median_cv_max, coherence_min, require_proteotypic, qvalue_max)
#   cv_reference : optional data.table(precursor_id, median_cv) to substitute frozen CV on discovery.
qualify_apply <- function(metrics, thr, cv_reference = NULL) {
  d <- copy(metrics)
  if (!is.null(cv_reference)) {
    d[, median_cv := NULL]
    d <- merge(d, cv_reference[, .(precursor_id, median_cv)], by = "precursor_id", all.x = TRUE)
  }
  dmin <- thr$detection_rate_min %||% 0.50
  cvmax <- thr$median_cv_max %||% 0.30
  cmin <- thr$coherence_min %||% 0.40
  qmax <- thr$qvalue_max %||% 0.01
  req_proteo <- isTRUE(thr$require_proteotypic %||% TRUE)
  d[, pass_detection := is.finite(detection_rate) & detection_rate >= dmin]
  d[, pass_cv := !is.finite(median_cv) | median_cv <= cvmax]          # NA CV (no reps) does not fail
  d[, is_singleton := n_siblings == 0L]
  d[, pass_coherence := is_singleton | (is.finite(coherence) & coherence >= cmin)]
  d[, pass_qvalue := !is.finite(median_qvalue) | median_qvalue <= qmax]
  d[, pass_proteotypic := if (!req_proteo) TRUE else (is.na(is_proteotypic) | is_proteotypic)]
  d[, qualifies := pass_detection & pass_cv & pass_coherence & pass_qvalue & pass_proteotypic]
  d[, qualify_flag := fifelse(is_singleton & qualifies, "pass_singleton",
                       fifelse(qualifies, "pass", "fail"))]
  d[]
}

# Threshold sensitivity sweep: how many precursors / proteins survive across a grid of detection and
# coherence cutoffs (holding the others at their base values). Feeds the Stage-1 decision report.
threshold_sensitivity <- function(metrics, thr_base,
                                  detection_grid = seq(0.3, 0.7, by = 0.1),
                                  coherence_grid = seq(0.2, 0.6, by = 0.1)) {
  rows <- CJ(detection_rate_min = detection_grid, coherence_min = coherence_grid)
  res <- rbindlist(lapply(seq_len(nrow(rows)), function(i) {
    thr <- modifyList(thr_base, list(detection_rate_min = rows$detection_rate_min[i],
                                     coherence_min = rows$coherence_min[i]))
    q <- qualify_apply(metrics, thr)
    data.table(detection_rate_min = rows$detection_rate_min[i],
               coherence_min = rows$coherence_min[i],
               n_precursors_pass = q[qualifies == TRUE, .N],
               n_proteins_pass = q[qualifies == TRUE, uniqueN(pg_canonical)])
  }))
  res[]
}

`%||%` <- function(a, b) if (is.null(a)) b else a
