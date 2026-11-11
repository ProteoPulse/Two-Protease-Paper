# detection.R — the DETECTION / COMPLEMENTARITY layer (protease-exclusive proteins).
#
# The concordance engine only judges the SHARED proteins (both proteases quantify them). This layer
# turns the discarded trypsin-only and HTA-only sets into signal and asks the complementary
# question: WHY does each protease see a different proteome — sensitivity or digestion chemistry?
# It also measures per-sample detection agreement and the coverage a second protease adds.
#
# Four pieces (see build_detection_report.R for the driver):
#   1. detection_table / exclusivity_tiers — per-protein detection counts per channel over the
#      shared sample universe; robust vs marginal (near-limit) exclusivity.
#   2. tryptic_peptide_count — in-silico trypsin digest: # detectable-length peptides (the trypsin
#      "sensitivity" covariate; HTA cleaves non-specifically so its detectability is length-driven).
#   3. detection_concordance — per-protein phi / Jaccard agreement of the two channels' 0/1
#      per-sample detection vectors (missingness as signal).
#   4. (coverage accounting is a small aggregation done in the driver).
#
# Depends on: data.table.

suppressPackageStartupMessages({ library(data.table) })

# Per-protein detection over the shared sample universe, from the two channel long tables
# (each already restricted to the shared sample columns; rows are finite-log2 detections).
detection_table <- function(A_long, B_long, n_samples, labelA = "trypsin", labelB = "hta") {
  da <- A_long[, .(n_a = uniqueN(sample_id)), by = pg_canonical]
  db <- B_long[, .(n_b = uniqueN(sample_id)), by = pg_canonical]
  d <- merge(da, db, by = "pg_canonical", all = TRUE)
  d[is.na(n_a), n_a := 0L][is.na(n_b), n_b := 0L]
  d[, `:=`(completeness_a = n_a / n_samples, completeness_b = n_b / n_samples,
           n_samples = n_samples)]
  d[, detection_class := fifelse(n_a > 0 & n_b > 0, "shared",
                          fifelse(n_a > 0 & n_b == 0, paste0(labelA, "_only"),
                                  paste0(labelB, "_only")))]
  d[]
}

# Robust vs marginal exclusivity: an exclusive protein is only "robust" if the detecting channel
# sees it in >= robust_frac of samples (else it is a near-detection-limit / stochastic call).
exclusivity_tiers <- function(det, robust_frac = 0.50, labelA = "trypsin", labelB = "hta") {
  d <- copy(det); ns <- d$n_samples[1]
  d[, exclusivity_tier := "shared"]
  d[detection_class == paste0(labelA, "_only"),
    exclusivity_tier := fifelse(n_a >= robust_frac * ns, "robust_exclusive", "marginal_exclusive")]
  d[detection_class == paste0(labelB, "_only"),
    exclusivity_tier := fifelse(n_b >= robust_frac * ns, "robust_exclusive", "marginal_exclusive")]
  d[]
}

# In-silico trypsin digest: number of fully-tryptic peptides in a detectable length window.
# Cleave C-terminal to K/R except when followed by P. Sequence-derived (channel-independent).
tryptic_peptide_count <- function(seq, min_len = 7, max_len = 30) {
  if (is.na(seq) || !nchar(seq)) return(NA_integer_)
  s <- gsub("([KR])(?!P)", "\\1|", seq, perl = TRUE)
  L <- nchar(strsplit(s, "|", fixed = TRUE)[[1]])
  sum(L >= min_len & L <= max_len)
}

# phi coefficient + Jaccard between two 0/1 detection vectors aligned by sample.
.phi_jaccard <- function(av, bv) {
  n11 <- sum(av & bv); n10 <- sum(av & !bv); n01 <- sum(!av & bv); n00 <- sum(!av & !bv)
  denom <- sqrt(prod(c(n11 + n10, n01 + n00, n11 + n01, n10 + n00)))
  phi <- if (denom > 0) (n11 * n00 - n10 * n01) / denom else NA_real_
  jac <- if ((n11 + n10 + n01) > 0) n11 / (n11 + n10 + n01) else NA_real_
  list(phi = phi, jaccard = jac, n_both = n11, n_either = n11 + n10 + n01)
}

# Per-protein per-sample detection concordance for a set of proteins over the shared sample list.
detection_concordance <- function(A_long, B_long, pgs, samples) {
  detA <- A_long[, .(pg_canonical, sample_id, a = 1L)]
  detB <- B_long[, .(pg_canonical, sample_id, b = 1L)]
  rows <- vector("list", length(pgs))
  for (i in seq_along(pgs)) {
    pg <- pgs[i]
    av <- samples %in% detA[pg_canonical == pg, sample_id]
    bv <- samples %in% detB[pg_canonical == pg, sample_id]
    r <- .phi_jaccard(av, bv)
    rows[[i]] <- data.table(pg_canonical = pg, det_phi = r$phi, det_jaccard = r$jaccard,
                            n_both = r$n_both, n_either = r$n_either)
  }
  rbindlist(rows, fill = TRUE)
}
