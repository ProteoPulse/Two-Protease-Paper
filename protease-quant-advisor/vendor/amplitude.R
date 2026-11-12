# amplitude.R — AMPLITUDE-agreement axis for the concordance classifier.
#
# The rank axis (Kendall tau-b / Spearman rho) is scale-free: it asks only whether the two
# proteases ORDER the samples the same way. A protein whose HTA channel is dynamically
# compressed (nearly flat across samples) can still earn a modest rank correlation while its
# fold-changes disagree — the "looks concordant on rank but not on intensity" failure mode.
#
# This module adds an amplitude-aware agreement axis on the within-channel MEAN-CENTERED log2
# profiles. Centering each channel separately removes the incomparable absolute-intensity offset
# between proteases (different peptides -> different absolute levels) while PRESERVING the
# log2 fold-change spread, which IS on a common scale. Two metrics:
#
#   1. Lin's concordance correlation coefficient (CCC) on the centered profiles. With per-channel
#      centering the location term vanishes, so CCC = 2*cov / (var_a + var_b): it rewards agreeing
#      fold-changes and PENALIZES dynamic-range mismatch (unlike Spearman/Pearson).
#   2. Deming (errors-in-variables) slope of B on A. Both channels carry measurement error, so an
#      ordinary least-squares slope is attenuated; Deming (here with error-variance ratio lambda=1)
#      is the symmetric estimate. slope ~ 1 => matched amplitude; slope -> 0 => B compressed.
#
# Diagnostic/annotation only — never mutates the underlying quantities.
# Depends on: utils_compare.R (build_pairs).

suppressPackageStartupMessages({ library(data.table) })

# Lin's CCC on two aligned vectors after per-vector mean-centering (location forgiven, scale kept).
lin_ccc_centered <- function(a, b) {
  ok <- is.finite(a) & is.finite(b); a <- a[ok]; b <- b[ok]; n <- length(a)
  if (n < 3) return(NA_real_)
  a <- a - mean(a); b <- b - mean(b)
  va <- stats::var(a); vb <- stats::var(b)
  if (!is.finite(va) || !is.finite(vb) || (va + vb) == 0) return(NA_real_)
  cab <- stats::cov(a, b)
  2 * cab / (va + vb)             # centered => mean-difference term is 0
}

# Deming slope of b ~ a with error-variance ratio lambda = var(err_a)/var(err_b) (default 1).
deming_slope <- function(a, b, lambda = 1) {
  ok <- is.finite(a) & is.finite(b); a <- a[ok]; b <- b[ok]; n <- length(a)
  if (n < 3) return(NA_real_)
  sxx <- stats::var(a); syy <- stats::var(b); sxy <- stats::cov(a, b)
  if (!is.finite(sxy) || sxy == 0) return(NA_real_)
  (syy - lambda * sxx + sqrt((syy - lambda * sxx)^2 + 4 * lambda * sxy^2)) / (2 * sxy)
}

# Per-protein amplitude metrics for a comparison. Returns pg_canonical + ccc + deming_slope.
amplitude_metrics <- function(A_long, B_long, min_paired = 15) {
  pr <- build_pairs(A_long, B_long)
  if (!nrow(pr)) return(data.table(pg_canonical = character(), ccc = numeric(), deming_slope = numeric()))
  pr[, .(ccc = lin_ccc_centered(log2q_A, log2q_B),
         deming_slope = deming_slope(log2q_A, log2q_B),
         n_amp = sum(is.finite(log2q_A) & is.finite(log2q_B))),
     by = pg_canonical][n_amp >= min_paired | n_amp > 0]
}
