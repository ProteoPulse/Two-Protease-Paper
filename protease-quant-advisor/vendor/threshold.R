# threshold.R — empirical threshold derivation for concordance classification.
# Generalized from CP02. The permutation null is ALWAYS available (needs only the target
# comparison). Same-protease positive controls are OPTIONAL: supply them only when the input
# provides replicate/benchmark preparations of the same channel. Nothing is auto-locked; the
# module returns a candidate table and a recommended value for human review.
#
# Depends on: utils_stats.R, utils_compare.R (build_pairs, spearman_pair).

suppressPackageStartupMessages({ library(data.table) })

# Permutation null: per-protein sample-label shuffles of channel B against A (target comparison).
# Preserves each protein's paired values, n, ties and marginals; frozen seed. Returns all values.
# `method` selects the statistic (kendall tau-b or spearman rho) so the null matches the gate.
perm_null_rhos <- function(A_long, B_long, n_perm = 2000, min_paired = 15, seed = 20260729,
                           method = "kendall") {
  set.seed(seed)
  pr <- build_pairs(A_long, B_long)
  if (!nrow(pr)) return(numeric())
  out <- vector("list", 0L)
  for (pg in unique(pr$pg_canonical)) {
    d <- pr[pg_canonical == pg]
    if (nrow(d) < min_paired) next
    a <- d$log2q_A; b <- d$log2q_B; n <- length(a)
    rr <- numeric(n_perm)
    for (i in seq_len(n_perm)) rr[i] <- suppressWarnings(cor(a, b[sample.int(n, n)], method = method))
    out[[length(out) + 1L]] <- rr
  }
  if (length(out)) unlist(out) else numeric()
}

null_summary <- function(null_rhos) {
  data.table(quantile = c("q50","q90","q95","q99","max"),
             rho = c(median(null_rhos, na.rm = TRUE),
                     quantile(null_rhos, c(.90,.95,.99), names = FALSE, na.rm = TRUE),
                     max(null_rhos, na.rm = TRUE)))
}

# Build the threshold candidate table + operating characteristics.
#   target_rho        : numeric vector of target (cross-channel) Spearman rhos (evaluable proteins)
#   null_rhos         : permutation-null rho pool
#   positive_rho      : OPTIONAL numeric vector of same-channel positive-control rhos (the weaker /
#                       binding control family if two are available). NULL -> control-based
#                       candidates are omitted and the null q95 is the recommendation.
#   target_rho_valid  : OPTIONAL second target (e.g. cross-gradient) for class-consistency reporting,
#                       aligned by name to target_rho; used only descriptively.
# `null_percentile` (default 0.99) is the operational gate percentile of the permutation null.
threshold_candidates <- function(target_rho, null_rhos, positive_rho = NULL,
                                 null_percentile = 0.99) {
  grid <- seq(-1, 1, by = 0.01)
  has_pos <- !is.null(positive_rho) && length(positive_rho) > 0
  pct_int <- round(null_percentile * 100)
  pct_name <- sprintf("permutation_%dth_percentile", pct_int)
  fc_pct <- round((1 - null_percentile) * 100)
  oc <- rbindlist(lapply(grid, function(c) {
    dt <- data.table(
      threshold = c,
      null_specificity = mean(null_rhos < c, na.rm = TRUE),
      false_concordance_null = mean(null_rhos >= c, na.rm = TRUE),
      target_frac_concordant = mean(target_rho >= c, na.rm = TRUE),
      target_n_concordant = sum(target_rho >= c, na.rm = TRUE),
      target_n_discordant = sum(target_rho < c, na.rm = TRUE))
    if (has_pos) {
      dt[, pos_sensitivity := mean(positive_rho >= c, na.rm = TRUE)]
      dt[, false_discordance_weak_pos := mean(positive_rho < c, na.rm = TRUE)]
      dt[, balanced_accuracy := 0.5 * (mean(positive_rho >= c, na.rm = TRUE) +
                                       mean(null_rhos < c, na.rm = TRUE))]
    }
    dt
  }), fill = TRUE)

  null_pct <- quantile(null_rhos, null_percentile, names = FALSE, na.rm = TRUE)
  cand <- data.table(candidate = pct_name,
                     rho_threshold = null_pct,
                     basis = sprintf("%dth percentile of the permutation null (false-concordance <= %d%%)",
                                     pct_int, fc_pct))
  if (has_pos) {
    pos10 <- quantile(positive_rho, .10, names = FALSE, na.rm = TRUE)
    ba_opt <- oc[which.max(balanced_accuracy), threshold]
    sep_mid <- if (null_pct < pos10) mean(c(null_pct, pos10)) else NA_real_
    cand <- rbind(cand,
      data.table(candidate = "balanced_accuracy_optimum", rho_threshold = ba_opt,
                 basis = "max balanced accuracy (weaker positive control vs null)"),
      data.table(candidate = "positive_control_10th_percentile", rho_threshold = pos10,
                 basis = "10th percentile of the weaker same-channel positive control"),
      data.table(candidate = "separation_midpoint", rho_threshold = sep_mid,
                 basis = sprintf("midpoint(null%d, weakPos10) if controls separate from null", pct_int)))
  }
  cand <- rbind(cand, data.table(candidate = "fixed_0.50_reference", rho_threshold = 0.50,
                                 basis = "fixed reference (not auto-selected)"))
  cand <- merge(cand, oc, by.x = "rho_threshold", by.y = "threshold", all.x = TRUE)
  setorder(cand, -rho_threshold)

  # recommendation: if separated, the separation midpoint; else the null percentile (operational cut).
  rec_name <- if (has_pos && !is.na(cand[candidate=="separation_midpoint", rho_threshold][1]))
    "separation_midpoint" else pct_name
  rec_val <- cand[candidate == rec_name, rho_threshold][1]
  list(candidates = cand, oc = oc, null_pct = null_pct, null_percentile = null_percentile,
       null95 = quantile(null_rhos, .95, names = FALSE, na.rm = TRUE),
       recommended_candidate = rec_name, recommended_threshold = rec_val,
       controls_available = has_pos)
}
