# utils_compare.R — build paired per-protein profiles for a comparison and compute metrics.
# A "comparison" pairs channel A (matrix + preparation) against channel B (matrix + preparation),
# aligning on biological sample id, and (for cross-protease) matching proteins by canonical group.
# No imputation: only samples finite in BOTH channels are used per protein.

suppressPackageStartupMessages({ library(data.table) })

# Long log2 table: one row per (pg_canonical, sample_id) for a given matrix+preparation.
# `mat_long` is precomputed (see harmonize_to_long); we just subset.
subset_channel <- function(mat_long, preparation) {
  mat_long[prep == preparation, .(pg_canonical, sample_id, log2q)]
}

# Build paired table for a comparison.
#   A_long, B_long: long tables (pg_canonical, sample_id, log2q) already subset to a preparation.
#   join = "same_protein" (A,B from same matrix -> identical pg space) or
#          "cross_protease" (match identical pg_canonical across the two matrices).
build_pairs <- function(A_long, B_long) {
  m <- merge(A_long, B_long, by = c("pg_canonical","sample_id"),
             suffixes = c("_A","_B"), all = FALSE)
  m[is.finite(log2q_A) & is.finite(log2q_B)]
}

# Per-protein metrics for one comparison (Spearman primary + descriptive + optional CIs).
#
# FIX (vs. the original pipeline's known bug): per-channel detection counts n_a / n_b and
# completeness are computed from the FULL channel long tables (true per-channel detection over
# the shared sample universe), NOT from the already-paired table (where n_a == n_b == n_paired
# is degenerate). `n_samples` = size of the shared sample universe used for the comparison; when
# NULL it is inferred as the number of distinct samples appearing in either channel.
compare_metrics <- function(A_long, B_long, min_paired = 15, do_ci = FALSE,
                            n_boot = 2000, seed = 20260729, n_samples = NULL,
                            method = "kendall") {
  pr <- build_pairs(A_long, B_long)
  if (!nrow(pr)) return(data.table())
  setorder(pr, pg_canonical, sample_id)
  # true per-channel detection counts (finite log2q rows) per protein, over the shared universe
  det_a <- A_long[is.finite(log2q), .(n_a = uniqueN(sample_id)), by = pg_canonical]
  det_b <- B_long[is.finite(log2q), .(n_b = uniqueN(sample_id)), by = pg_canonical]
  if (is.null(n_samples))
    n_samples <- uniqueN(unique(c(A_long$sample_id, B_long$sample_id)))
  res <- pr[, {
    a <- log2q_A; b <- log2q_B
    # both rank statistics are always computed; `method` picks the canonical stat_value.
    kt <- stat_pair(a, b, method = "kendall")
    sp <- stat_pair(a, b, method = "spearman")
    prim <- if (identical(method, "spearman")) sp else kt
    pm <- pair_metrics(a, b)   # n_a/n_b here are paired-degenerate; overwritten below
    # CI + Fisher-z envelope track the PRIMARY statistic (bootstrap is method-aware).
    ci <- if (do_ci && prim$n >= 5) bootstrap_rho_ci(a, b, n_reps = n_boot, seed = seed, method = method) else c(low=NA_real_, high=NA_real_)
    fz <- fisher_z_ci(prim$value, prim$n)
    .(stat_value = prim$value, stat_method = method,
      kendall_tau = kt$value, spearman_rho = sp$value,
      pearson_r = pearson_pair(a, b), n_paired = prim$n,
      rho_ci_low = ci[["low"]], rho_ci_high = ci[["high"]],
      fisher_low = fz[["low"]], fisher_high = fz[["high"]],
      distinct_a = pm$distinct_a, distinct_b = pm$distinct_b,
      sd_a = pm$sd_a, sd_b = pm$sd_b, iqr_a = pm$iqr_a, iqr_b = pm$iqr_b,
      max_tie_a = pm$max_tie_a, max_tie_b = pm$max_tie_b,
      mean_a = pm$mean_a, mean_b = pm$mean_b,
      amp_log2_iqr_ratio = if (!is.na(pm$iqr_b) && pm$iqr_b > 0) pm$iqr_a/pm$iqr_b else NA_real_,
      median_abs_dz = median_abs_dz(a, b))
  }, by = pg_canonical]
  # attach TRUE per-channel detection + completeness (fix)
  res <- merge(res, det_a, by = "pg_canonical", all.x = TRUE)
  res <- merge(res, det_b, by = "pg_canonical", all.x = TRUE)
  res[is.na(n_a), n_a := 0L][is.na(n_b), n_b := 0L]
  res[, `:=`(completeness_a = n_a / n_samples,
             completeness_b = n_b / n_samples,
             abs_completeness_diff = abs(n_a - n_b) / n_samples,
             n_samples = n_samples)]
  res[]
}

# Convenience: n paired samples available for a comparison (per protein), used for eligibility.
eligibility_flags <- function(dt, min_paired = 15, min_distinct = 5) {
  dt[, `:=`(
    enough_paired = n_paired >= min_paired,
    enough_distinct = distinct_a >= min_distinct & distinct_b >= min_distinct,
    nonzero_var = is.finite(sd_a) & is.finite(sd_b) & sd_a > 0 & sd_b > 0
  )]
  dt[, concordance_evaluable := enough_paired & enough_distinct & nonzero_var]
  dt[, not_assessable_reasons := {
    r <- character(.N)
    r <- ifelse(!enough_paired, paste0(r, "few_paired;"), r)
    r <- ifelse(!enough_distinct, paste0(r, "too_few_distinct;"), r)
    r <- ifelse(!nonzero_var, paste0(r, "zero_variance;"), r)
    sub(";$", "", r)
  }]
  # low dynamic range is a FLAG, not an exclusion
  dt[, low_dynamic_range := is.finite(sd_a) & is.finite(sd_b) & (pmin(sd_a, sd_b) < 0.20)]
  dt[, n_tier := fifelse(n_paired >= 40, ">=40",
                  fifelse(n_paired >= 30, "30-39",
                  fifelse(n_paired >= 20, "20-29",
                  fifelse(n_paired >= 15, "15-19", "<15"))))]
  dt[]
}
