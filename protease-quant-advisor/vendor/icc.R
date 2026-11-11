# icc.R — one-way random-effects intraclass correlation (ICC(1)) for the reliability lookup.
#
# The reproducibility design measures each biological sample with many technical replicates
# (8 samples x 24 digests x 3 days). For a rolled-up protein estimate, ICC(1) is the fraction of
# total variance that is BETWEEN samples — i.e. how reliably the protease ranks samples through the
# technical noise. This is the true replicate-based reliability the Aug-18 plan calls for, and it
# REPLACES the split-half attenuation surrogate (reliability.R) when replicate data exist: it feeds
# classify_multiaxis() through the same ceiling / rel_A / rel_B columns.
#
# Estimator (one-way random ANOVA, unbalanced-safe):
#   ICC = (MSB - MSW) / (MSB + (k0 - 1) * MSW)
#   k0  = (N - sum(n_i^2)/N) / (m - 1)     # average group size correction for unbalanced designs
# Negative ICC (MSB < MSW: no between-sample signal) is clamped to 0.
#
# Depends on: data.table. (irr::icc is an optional external cross-check, not required.)

suppressPackageStartupMessages({ library(data.table) })

# ICC(1) for a vector of `values` grouped by `groups` (biological sample id).
# Returns NA if fewer than 2 groups or fewer than 3 finite observations.
icc_oneway <- function(values, groups) {
  ok <- is.finite(values) & !is.na(groups)
  values <- values[ok]; groups <- factor(groups[ok])
  N <- length(values); m <- nlevels(groups)
  if (N < 3 || m < 2) return(NA_real_)
  grand <- mean(values)
  n_i <- as.numeric(table(groups))
  grp_mean <- tapply(values, groups, mean)
  SSB <- sum(n_i * (grp_mean - grand)^2)
  SSW <- sum((values - grp_mean[as.character(groups)])^2)
  dfB <- m - 1; dfW <- N - m
  if (dfW <= 0) return(NA_real_)
  MSB <- SSB / dfB; MSW <- SSW / dfW
  k0 <- (N - sum(n_i^2) / N) / dfB
  if (!is.finite(k0) || k0 <= 0) return(NA_real_)
  icc <- (MSB - MSW) / (MSB + (k0 - 1) * MSW)
  if (!is.finite(icc)) return(NA_real_)
  max(0, min(1, icc))
}

# Per protein x protease ICC from a long per-run roll-up table.
#   roll_long  : data.table with pg_canonical, protease, sample_id (biological), value_col (theta)
#   value_col  : the rolled-up per-run estimate column (default "theta_hat")
#   labelA/labelB : protease labels to pivot into R_<labelA> / R_<labelB>
#   floor      : reliabilities below this are flagged (inherently noisy; r_adjusted -> NA later)
# Returns pg_canonical | R_<labelA> | R_<labelB> | ceiling | flag_low_A | flag_low_B | flag_low.
build_icc_lookup <- function(roll_long, value_col = "theta_hat",
                             labelA = "trypsin", labelB = "hta", floor = 0.30) {
  d <- roll_long[is.finite(get(value_col))]
  icc <- d[, .(R = icc_oneway(get(value_col), sample_id)), by = .(pg_canonical, protease)]
  wide <- dcast(icc, pg_canonical ~ protease, value.var = "R")
  ca <- labelA; cb <- labelB
  if (!ca %in% names(wide)) wide[, (ca) := NA_real_]
  if (!cb %in% names(wide)) wide[, (cb) := NA_real_]
  setnames(wide, c(ca, cb), c(paste0("R_", ca), paste0("R_", cb)))
  RA <- paste0("R_", ca); RB <- paste0("R_", cb)
  wide[, ceiling := fifelse(is.finite(get(RA)) & is.finite(get(RB)),
                            sqrt(pmax(0, get(RA)) * pmax(0, get(RB))), NA_real_)]
  wide[, flag_low_A := is.finite(get(RA)) & get(RA) < floor]
  wide[, flag_low_B := is.finite(get(RB)) & get(RB) < floor]
  wide[, flag_low := (is.finite(get(RA)) & get(RA) < floor) | (is.finite(get(RB)) & get(RB) < floor)]
  # rel_A / rel_B aliases so classify_multiaxis() consumes ICC exactly like the split-half ceiling
  wide[, rel_A := get(RA)][, rel_B := get(RB)]
  wide[]
}
