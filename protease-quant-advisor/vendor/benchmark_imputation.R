# benchmark_imputation.R — is cross-protease filling actually better than the standard methods?
#
# The only honest way to answer is to hide values we already have and see who gets them back. Cells
# observed in BOTH channels are masked from the recipient, imputed, and compared to the withheld truth.
#
# TWO MASKING REGIMES, both reported:
#   MCAR  a random mask_frac of co-observed recipient cells (repeated mcar_repeats times). Every
#         method looks good under MCAR — the missing values are ordinary ones — so this is the
#         flattering regime, and on its own it would overstate any method.
#   MNAR  the lowest-abundance mask_frac of each protein's co-observed cells: the left-censored
#         analogue of why proteomics values are actually missing. This is the regime that matches
#         the problem, and the one that decides the question.
#
# BASELINES: MinDet and MinProb (the classic downshifted MNAR draws), protein-mean, sample-median, and
# protein-kNN (the standard MAR method). Cross-protease has to beat these, not a straw man.
#
# NO LEAKAGE: every offset, mean, median and neighbour set is recomputed from the UNMASKED cells only.
#
# NEGATIVE CONTROL: the donor's sample labels are shuffled and the whole cross-protease procedure is
# rerun. If the method is really using the donor's per-sample signal, its error must collapse toward
# the baselines when that signal is destroyed. If it does NOT collapse, the method is exploiting
# protein-level structure (an average, essentially) and the "cross-protease" framing is a fiction.
#
# Depends on: data.table.

suppressPackageStartupMessages({ library(data.table) })

# ---- masking -----------------------------------------------------------------------------------
# X : long table pg_canonical, sample_id, theta_rec, theta_don (both finite = maskable)
mask_mcar <- function(X, frac, seed) {
  set.seed(seed)
  cand <- which(is.finite(X$theta_rec) & is.finite(X$theta_don))
  if (!length(cand)) return(integer(0))
  sample(cand, max(1L, floor(frac * length(cand))))
}
mask_mnar <- function(X, frac) {
  X2 <- copy(X)[, .I0 := .I]
  cand <- X2[is.finite(theta_rec) & is.finite(theta_don)]
  if (!nrow(cand)) return(integer(0))
  cand[, r := frank(theta_rec, ties.method = "first"), by = pg_canonical]
  cand[, n := .N, by = pg_canonical]
  cand[r <= pmax(1L, floor(frac * n)), .I0]
}

# ---- imputers ------------------------------------------------------------------------------------
# Each returns a numeric vector aligned to `target` rows of X. `train` marks the cells an imputer may
# learn from (observed AND not masked).
.imp_cross_protease <- function(X, train, target, min_paired) {
  tr <- X[train == TRUE & is.finite(theta_rec) & is.finite(theta_don)]
  off <- tr[, .(offset = stats::median(theta_rec - theta_don), n = .N), by = pg_canonical]
  off <- off[n >= min_paired]
  m <- merge(X[target, .(pg_canonical, sample_id, theta_don)], off, by = "pg_canonical", all.x = TRUE,
             sort = FALSE)
  m$theta_don + m$offset
}
.imp_protein_mean <- function(X, train, target) {
  pm <- X[train == TRUE & is.finite(theta_rec), .(mu = mean(theta_rec)), by = pg_canonical]
  merge(X[target, .(pg_canonical)], pm, by = "pg_canonical", all.x = TRUE, sort = FALSE)$mu
}
.imp_sample_median <- function(X, train, target) {
  sm <- X[train == TRUE & is.finite(theta_rec), .(md = stats::median(theta_rec)), by = sample_id]
  merge(X[target, .(sample_id)], sm, by = "sample_id", all.x = TRUE, sort = FALSE)$md
}
.imp_mindet <- function(X, train, target, q = 0.01) {
  sm <- X[train == TRUE & is.finite(theta_rec),
          .(v = stats::quantile(theta_rec, probs = q, names = FALSE)), by = sample_id]
  merge(X[target, .(sample_id)], sm, by = "sample_id", all.x = TRUE, sort = FALSE)$v
}
.imp_minprob <- function(X, train, target, shift_sd, width_sd, seed) {
  set.seed(seed)
  sm <- X[train == TRUE & is.finite(theta_rec),
          .(mu = mean(theta_rec), sdv = stats::sd(theta_rec)), by = sample_id]
  m <- merge(X[target, .(sample_id)], sm, by = "sample_id", all.x = TRUE, sort = FALSE)
  m[, sdv := fifelse(is.finite(sdv) & sdv > 0, sdv, 0.1)]
  stats::rnorm(nrow(m), mean = m$mu - shift_sd * m$sdv, sd = width_sd * m$sdv)
}
# protein-kNN on the recipient channel alone: neighbours by correlation over co-observed samples,
# donor value shifted by the protein-mean difference (the standard MAR construction).
.imp_protein_knn <- function(X, train, target, k) {
  w <- dcast(X[train == TRUE & is.finite(theta_rec)], pg_canonical ~ sample_id, value.var = "theta_rec")
  pgs <- w$pg_canonical; M <- as.matrix(w[, -1, drop = FALSE]); rownames(M) <- pgs
  if (nrow(M) < 3L) return(rep(NA_real_, length(target)))
  C <- suppressWarnings(stats::cor(t(M), use = "pairwise.complete.obs"))
  mu <- rowMeans(M, na.rm = TRUE)
  tg <- X[target, .(pg_canonical, sample_id)]
  vapply(seq_len(nrow(tg)), function(i) {
    pg <- tg$pg_canonical[i]; s <- as.character(tg$sample_id[i])
    if (!(pg %in% rownames(C)) || !(s %in% colnames(M))) return(NA_real_)
    sim <- C[pg, ]; sim[pg] <- NA_real_
    col <- M[, s]
    ok <- is.finite(sim) & is.finite(col)
    if (!any(ok)) return(NA_real_)
    idx <- order(-sim[ok])[seq_len(min(k, sum(ok)))]
    nb <- names(sim[ok])[idx]
    mean(col[nb] + (mu[pg] - mu[nb]), na.rm = TRUE)
  }, numeric(1))
}

# ---- one masking experiment ----------------------------------------------------------------------
.evaluate <- function(truth, pred) {
  ok <- is.finite(truth) & is.finite(pred)
  if (!any(ok)) return(data.table(n_cells = 0L, coverage = 0, median_abs_error = NA_real_,
                                  rmse = NA_real_, median_signed_error = NA_real_, spearman = NA_real_))
  e <- pred[ok] - truth[ok]
  data.table(n_cells = sum(ok), coverage = sum(ok) / length(truth),
             median_abs_error = stats::median(abs(e)), rmse = sqrt(mean(e^2)),
             median_signed_error = stats::median(e),
             spearman = if (sum(ok) > 3 && stats::sd(truth[ok]) > 0 && stats::sd(pred[ok]) > 0)
               suppressWarnings(stats::cor(truth[ok], pred[ok], method = "spearman")) else NA_real_)
}

run_one_mask <- function(X, mask_idx, cfg_b, min_paired, seed, shuffle_donor = FALSE) {
  Y <- copy(X)
  if (shuffle_donor) {                    # negative control: destroy the donor's sample-level signal
    set.seed(seed + 7L)
    smap <- data.table(sample_id = unique(Y$sample_id))
    smap[, shuffled := sample(sample_id)]
    Yd <- Y[, .(pg_canonical, sample_id, theta_don)]
    Yd <- merge(Yd, smap, by = "sample_id")[, .(pg_canonical, sample_id = shuffled, theta_don)]
    Y[, theta_don := NULL]
    Y <- merge(Y, Yd, by = c("pg_canonical","sample_id"), all.x = TRUE)
  }
  truth <- Y$theta_rec[mask_idx]
  Y[, train := TRUE]; Y[mask_idx, train := FALSE]
  Yh <- copy(Y); Yh[mask_idx, theta_rec := NA_real_]
  res <- list()
  res$cross_protease <- .imp_cross_protease(Yh, Yh$train, mask_idx, min_paired)
  bl <- unlist(cfg_b$baselines)
  if ("protein_mean"  %in% bl) res$protein_mean  <- .imp_protein_mean(Yh, Yh$train, mask_idx)
  if ("sample_median" %in% bl) res$sample_median <- .imp_sample_median(Yh, Yh$train, mask_idx)
  if ("mindet"        %in% bl) res$mindet        <- .imp_mindet(Yh, Yh$train, mask_idx)
  if ("minprob"       %in% bl) res$minprob       <- .imp_minprob(Yh, Yh$train, mask_idx,
                                                     cfg_b$minprob_shift_sd, cfg_b$minprob_width_sd, seed)
  if ("protein_knn"   %in% bl) res$protein_knn   <- .imp_protein_knn(Yh, Yh$train, mask_idx, 10L)
  rbindlist(lapply(names(res), function(nm) cbind(data.table(method = nm), .evaluate(truth, res[[nm]]))))
}

# ---- full benchmark --------------------------------------------------------------------------------
run_benchmark <- function(X, cfg_b, min_paired = 8L, seed = 20260820) {
  out <- list()
  if ("mcar" %in% unlist(cfg_b$regimes)) {
    reps <- as.integer(cfg_b$mcar_repeats)
    for (r in seq_len(reps)) {
      mi <- mask_mcar(X, cfg_b$mask_frac, seed + r)
      if (!length(mi)) next
      out[[length(out)+1L]] <- cbind(regime = "mcar", rep = r, run_one_mask(X, mi, cfg_b, min_paired, seed + r))
      if (isTRUE(cfg_b$negative_control) && r <= 5L)
        out[[length(out)+1L]] <- cbind(regime = "mcar_shuffled_donor", rep = r,
                                       run_one_mask(X, mi, cfg_b, min_paired, seed + r, shuffle_donor = TRUE)[method == "cross_protease"])
    }
  }
  if ("mnar" %in% unlist(cfg_b$regimes)) {
    mi <- mask_mnar(X, cfg_b$mask_frac)
    if (length(mi)) {
      out[[length(out)+1L]] <- cbind(regime = "mnar", rep = 1L, run_one_mask(X, mi, cfg_b, min_paired, seed))
      if (isTRUE(cfg_b$negative_control))
        out[[length(out)+1L]] <- cbind(regime = "mnar_shuffled_donor", rep = 1L,
                                       run_one_mask(X, mi, cfg_b, min_paired, seed, shuffle_donor = TRUE)[method == "cross_protease"])
    }
  }
  raw <- rbindlist(out, fill = TRUE)
  summ <- raw[, .(n_reps = .N, n_cells = as.numeric(stats::median(as.numeric(n_cells))), coverage = mean(coverage),
                  median_abs_error = mean(median_abs_error, na.rm = TRUE),
                  mae_sd = stats::sd(median_abs_error, na.rm = TRUE),
                  rmse = mean(rmse, na.rm = TRUE),
                  median_signed_error = mean(median_signed_error, na.rm = TRUE),
                  spearman = mean(spearman, na.rm = TRUE)), by = .(regime, method)]
  setorder(summ, regime, median_abs_error)
  list(raw = raw, summary = summ)
}

# ---- post-imputation QC ------------------------------------------------------------------------------
# Does filling change each protein's spread? A method that compresses dynamic range (everything pulled
# to a mean) or spikes it (wild donor translations) shows up here even when the error metrics look fine.
imputation_qc <- function(mat) {
  before <- mat[value_source == "measured", .(iqr_before = stats::IQR(value, na.rm = TRUE),
                                              sd_before = stats::sd(value, na.rm = TRUE),
                                              n_before = .N), by = pg_canonical]
  after <- mat[, .(iqr_after = stats::IQR(value, na.rm = TRUE),
                   sd_after = stats::sd(value, na.rm = TRUE), n_after = .N), by = pg_canonical]
  d <- merge(before, after, by = "pg_canonical")
  d[, `:=`(delta_iqr = iqr_after - iqr_before, delta_sd = sd_after - sd_before,
           n_imputed = n_after - n_before)]
  d[]
}

# Technical-CV inflation (reproducibility only): CV across replicate runs before vs after imputation.
# Imputation happens at SAMPLE grain, so an imputed cell has no within-sample replicate spread; this
# reports the protein-level CV over samples instead, and states the count it is based on.
cv_inflation <- function(mat) {
  f <- function(sub) { lin <- 2^sub$value; if (length(lin) < 2 || mean(lin) <= 0) NA_real_ else stats::sd(lin)/mean(lin) }
  b <- mat[value_source == "measured", .(cv_before = f(.SD)), by = pg_canonical, .SDcols = "value"]
  a <- mat[, .(cv_after = f(.SD), n_imputed = sum(value_source != "measured")), by = pg_canonical, .SDcols = c("value","value_source")]
  merge(b, a, by = "pg_canonical")
}
