# utils_stats.R — transforms, correlation, CIs, permutation null. No imputation anywhere.

suppressPackageStartupMessages({ library(data.table) })

# ---- transform: linear -> log2 with zero/blank/neg/non-finite -> NA ----------
# Auto-detect whether a matrix is already log-transformed (heuristic on positive values).
detect_is_log <- function(x) {
  v <- x[is.finite(x)]; v <- v[v != 0]
  if (!length(v)) return(FALSE)
  # linear MaxLFQ intensities span many orders of magnitude and are large/positive;
  # log2 values are typically within ~[-5, 45] and include small magnitudes.
  stats::quantile(abs(v), 0.99, names = FALSE) < 100 && min(v) > -30
}
to_log2 <- function(x, already_log = NA) {
  x <- as.numeric(x)
  x[!is.finite(x) | x <= 0] <- NA_real_
  if (is.na(already_log)) already_log <- detect_is_log(x)
  if (already_log) x else log2(x)
}

# ---- rank correlation on complete pairs (method = "kendall" | "spearman") ------
# Generalized from spearman_pair. Kendall returns tau-b (tie-corrected) via stats::cor.
# Returns list(value, n); `rho` is kept as an alias for back-compat with callers.
stat_pair <- function(a, b, method = "kendall") {
  ok <- is.finite(a) & is.finite(b)
  n <- sum(ok)
  if (n < 3) return(list(value = NA_real_, rho = NA_real_, n = n))
  v <- suppressWarnings(stats::cor(a[ok], b[ok], method = method))
  list(value = v, rho = v, n = n)
}
# Back-compat thin wrapper (Spearman): existing callers / tests keep working.
spearman_pair <- function(a, b) {
  s <- stat_pair(a, b, method = "spearman")
  list(rho = s$value, n = s$n)
}
pearson_pair <- function(a, b) {
  ok <- is.finite(a) & is.finite(b); n <- sum(ok)
  if (n < 3) return(NA_real_)
  suppressWarnings(stats::cor(a[ok], b[ok], method = "pearson"))
}

# ---- Fisher-z CI envelope for a single rho at given n ------------------------
fisher_z_ci <- function(rho, n, level = 0.95) {
  if (is.na(rho) || n < 4) return(c(low = NA_real_, high = NA_real_))
  z <- atanh(pmin(pmax(rho, -0.999999), 0.999999))
  se <- 1 / sqrt(n - 3)
  zc <- stats::qnorm(1 - (1 - level) / 2)
  c(low = tanh(z - zc * se), high = tanh(z + zc * se))
}

# ---- descriptive eligibility metrics for one paired protein profile ----------
pair_metrics <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  ao <- a[ok]; bo <- b[ok]; n <- sum(ok)
  frac_tie <- function(v) if (length(v)) max(table(v)) / length(v) else NA_real_
  list(
    n_paired = n,
    n_a = sum(is.finite(a)), n_b = sum(is.finite(b)),
    distinct_a = length(unique(ao)), distinct_b = length(unique(bo)),
    sd_a = if (n > 1) stats::sd(ao) else NA_real_, sd_b = if (n > 1) stats::sd(bo) else NA_real_,
    iqr_a = if (n > 1) stats::IQR(ao) else NA_real_, iqr_b = if (n > 1) stats::IQR(bo) else NA_real_,
    max_tie_a = frac_tie(ao), max_tie_b = frac_tie(bo),
    mean_a = if (n) mean(ao) else NA_real_, mean_b = if (n) mean(bo) else NA_real_
  )
}

# ---- bootstrap CI over SAMPLES (jointly resample sample IDs, recompute statistic) ----
# a, b are aligned per-sample vectors for ONE protein. Resampling sample indices
# with replacement respects that protein statistics share the same samples. `method`
# selects the statistic (kendall tau-b or spearman rho) so the CI matches the gate.
bootstrap_rho_ci <- function(a, b, n_reps = 2000, level = 0.95, seed = NULL, method = "kendall") {
  if (!is.null(seed)) set.seed(seed)
  ok <- is.finite(a) & is.finite(b); ao <- a[ok]; bo <- b[ok]; n <- length(ao)
  if (n < 5) return(c(low = NA_real_, high = NA_real_))
  reps <- numeric(n_reps)
  for (i in seq_len(n_reps)) {
    idx <- sample.int(n, n, replace = TRUE)
    reps[i] <- suppressWarnings(stats::cor(ao[idx], bo[idx], method = method))
  }
  a2 <- (1 - level) / 2
  stats::quantile(reps, c(a2, 1 - a2), na.rm = TRUE, names = FALSE) |>
    setNames(c("low", "high"))
}

# ---- amplitude annotation: log2 IQR ratio (from ORIGINAL log2, not z-scored) --
log2_iqr_ratio <- function(a, b) {
  ma <- pair_metrics(a, b)
  if (is.na(ma$iqr_a) || is.na(ma$iqr_b) || ma$iqr_b == 0) return(NA_real_)
  ma$iqr_a / ma$iqr_b
}

# ---- diagnostic standardized discrepancy (z-scored; VISUALIZATION/diag only) --
zprofile <- function(x) {
  ok <- is.finite(x)
  if (sum(ok) < 2 || stats::sd(x[ok]) == 0) return(rep(NA_real_, length(x)))
  m <- mean(x[ok]); s <- stats::sd(x[ok]); z <- (x - m) / s; z[!ok] <- NA_real_; z
}
median_abs_dz <- function(a, b) {
  d <- zprofile(b) - zprofile(a)
  stats::median(abs(d), na.rm = TRUE)
}
