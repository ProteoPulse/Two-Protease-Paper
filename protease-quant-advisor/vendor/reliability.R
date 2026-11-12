# reliability.R — the RELIABILITY GATE for the concordance classifier.
#
# A cross-protease correlation can only be as high as each channel is internally reproducible.
# If HTA measures a protein with a single ragged peptide, the HTA profile is mostly noise and NO
# honest concordance call is possible — such a protein must be labelled "measurement_limited",
# not "discordant". This module estimates, per protein and channel, how reproducibly the protease
# even measures that protein, and turns it into a CEILING on the attainable cross-protease
# correlation.
#
# Estimator: split-half reliability with the Spearman-Brown prophecy step, computed from the core
# peptide matrix (sample x peptide log2). Repeatedly split a protein's core peptides into two
# halves, build a surrogate profile from each half (median of z-scored peptides), correlate the
# two half-profiles across samples, and step the half-length correlation up to the full peptide
# set: rel = 2 r_half / (1 + r_half). Median over random splits. This is a peptide-RESAMPLING
# reliability: it measures whether the protein's own peptides agree on the sample ordering.
#
#   - >= 2 core peptides required; a single peptide cannot self-validate -> rel = NA -> the protein
#     is treated as measurement-limited (exactly the HTA-sparse case).
#   - median pairwise peptide coherence is kept as a cheaper fallback / cross-check.
#
# Ceiling: ceiling = sqrt(max(0,rel_A) * max(0,rel_B)) — the classic attenuation bound on the
# correlation between two noisily-measured variables. Observed cross-protease rho is compared to
# this ceiling in classify_multiaxis.R (the disattenuated concordance score).
#
# Diagnostic only. Depends on: diagnostics_peptide.R (core_matrix, zcols, surrogate_profile,
# median_pairwise_coherence) and utils_stats.R.

suppressPackageStartupMessages({ library(data.table) })

# Spearman correlation of two per-sample vectors on their finite overlap (>= min_paired).
.prof_cor <- function(x, y, min_paired = 8) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < min_paired) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
}

# Split-half + Spearman-Brown reliability for one sample x peptide matrix m.
channel_reliability <- function(m, n_splits = 200, seed = 20260729, min_paired = 8) {
  if (is.null(m) || ncol(m) < 2) return(NA_real_)
  k <- ncol(m); mz <- zcols(m)
  set.seed(seed)
  rr <- numeric(n_splits)
  for (i in seq_len(n_splits)) {
    idx <- sample.int(k, k)                      # random peptide permutation
    h1 <- idx[seq_len(floor(k / 2))]; h2 <- idx[(floor(k / 2) + 1L):k]
    p1 <- surrogate_profile(mz[, h1, drop = FALSE])
    p2 <- surrogate_profile(mz[, h2, drop = FALSE])
    rr[i] <- .prof_cor(p1, p2, min_paired)
  }
  r_half <- stats::median(rr, na.rm = TRUE)
  if (!is.finite(r_half)) return(NA_real_)
  r_half <- max(min(r_half, 0.999), -0.999)
  r_full <- 2 * r_half / (1 + r_half)            # Spearman-Brown: half-length -> full-length
  max(0, min(1, r_full))
}

# Per-protein/channel reliability + attenuation ceiling for a set of protein groups.
#   pepA, pepB : load_peptides_long() output for channels A, B (all shared proteins).
# Returns pg_canonical + rel_A/rel_B, coherence_rel_A/B (fallback), n_core_A/B, ceiling.
reliability_ceiling <- function(pepA, pepB, pgs, core_frac = 0.50, min_obs = 20,
                                n_splits = 200, seed = 20260729, min_paired = 8) {
  rows <- vector("list", length(pgs))
  for (i in seq_along(pgs)) {
    pg <- pgs[i]
    mA <- core_matrix(pepA, pg, core_frac, min_obs)
    mB <- core_matrix(pepB, pg, core_frac, min_obs)
    kA <- if (is.null(mA)) 0L else ncol(mA); kB <- if (is.null(mB)) 0L else ncol(mB)
    relA <- channel_reliability(mA, n_splits, seed, min_paired)
    relB <- channel_reliability(mB, n_splits, seed, min_paired)
    cohA <- if (kA >= 2) median_pairwise_coherence(mA) else NA_real_
    cohB <- if (kB >= 2) median_pairwise_coherence(mB) else NA_real_
    # fall back to coherence if the split-half estimate is unavailable but >=2 peptides exist
    urA <- if (is.finite(relA)) relA else cohA
    urB <- if (is.finite(relB)) relB else cohB
    ceiling <- if (is.finite(urA) && is.finite(urB)) sqrt(max(0, urA) * max(0, urB)) else NA_real_
    rows[[i]] <- data.table(pg_canonical = pg, rel_A = relA, rel_B = relB,
                            coherence_rel_A = cohA, coherence_rel_B = cohB,
                            n_core_A = kA, n_core_B = kB, ceiling = ceiling)
  }
  rbindlist(rows, fill = TRUE)
}
