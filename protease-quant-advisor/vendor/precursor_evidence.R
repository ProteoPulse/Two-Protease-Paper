# precursor_evidence.R — Part A core: per protein x protease evidence on two axes.
#
# WHY THIS EXISTS. The Aug-19 framework reconciled the two proteases into ONE number per shared
# protein (pooled peptides or IVW). Both channels come from the same plasma injection searched against
# two libraries, so that reconciliation does not combine independent measurements — and the result was
# silently trypsin (pooled sat closer to trypsin in 140/150 Mix30 proteins). Phase 3 therefore stops
# reconciling and starts CHOOSING: measure each channel's precursor-level evidence and recommend the
# channel whose roll-up is better supported.
#
# Two axes, each a family of metrics (see the recommendation weights in config.yml):
#
#   COMPLETENESS — is the protein actually measured everywhere we need it?
#     n_prec_present / n_prec_qual / n_peptides_qual / n_features
#     cell_fill_rate          observed (precursor x unit) cells / (n_prec_qual * n_units)
#     protein_detection_rate  fraction of SAMPLES with >= 1 qualifying precursor
#     frac_samples_ge2_prec   fraction of samples where the protein can self-validate (>= 2 precursors)
#     median_prec_detection   median per-precursor detection rate
#     peptide_set_jaccard     median pairwise Jaccard of per-unit observed peptide sets. Prior work
#                             (HTA_Peptide_Heterogeneity_20260731) put trypsin ~0.87 vs HTA ~0.49:
#                             HTA measures a DIFFERENT peptide set in each run, so its roll-up is not
#                             the same estimator sample to sample. This metric is that effect.
#     seq_coverage_frac / median_unit_coverage_frac   FASTA-gated residue coverage (NA without FASTA)
#
#   COHERENCE — do the peptides that support the protein tell the same story?
#     median_sibling_coherence / frac_prec_coherent   per-precursor median Pearson to its siblings
#     splithalf_reliability   split-half + Spearman-Brown over the protein's own feature matrix
#     lopo_rank_stability     leave-one-peptide-out: re-polish with each feature dropped, Spearman to
#                             the full estimate across samples, MINIMISED over dropped features. Low =>
#                             the protein's number hangs on one peptide.
#     lopo_max_delta          max median |change| in the estimate under LOPO (offset-centred, log2)
#     median_rollup_se        median per-sample median-polish SE
#     frozen_technical_cv     median frozen reproducibility CV over the protein's qualifying precursors
#     icc_repro               frozen replicate ICC (reproducibility), where the protein has one
#
# GRAIN. Presence metrics (cell fill, peptide-set stability, per-precursor detection) are computed at
# the observation unit = RUN; profile metrics (split-half, LOPO, SE) at the SAMPLE grain, replicate
# runs collapsed by median. In discovery run == sample, so the two coincide; in reproducibility they
# differ and each metric uses the grain that makes it mean what its name says.
#
# COHERENCE SCOPE. Coherence is measured over ALL present precursors of the protein, not only the
# qualifying ones: frozen qualification already gates on coherence, so within the qualifying subset
# the metric is near-degenerate and cannot discriminate between channels. The qualifying-subset value
# is emitted alongside as *_qual for transparency.
#
# PER-PROTEASE CORRECTNESS. The frozen Stage-1 computed per-precursor metrics grouped by
# precursor_id ALONE, so a precursor id seen in both searches had its two channels' metrics merged.
# That is inherited as-is for qualification (frozen means frozen), but every metric here is computed
# per (protein x protease), and channel_independence.R reports how many ids actually collide.
#
# Depends on: data.table; vendored feature_matrix.R, rollup_tmp.R (rollup_protein_tmp),
# reliability.R + diagnostics_peptide.R (channel_reliability, zcols, surrogate_profile), seq_coverage.R.

suppressPackageStartupMessages({ library(data.table) })

# ---- sibling coherence, computed per (protein x protease) ------------------------------------
# Same estimator as qualify_precursors.R::precursor_metrics (median pairwise Pearson of a precursor
# to its siblings across units) but grouped by protein AND protease, and done group-wise so it costs
# one pass instead of one full-table scan per protein.
sibling_coherence_by_channel <- function(prec, unit_col = "run") {
  d <- prec[is.finite(log2q)]
  d[, {
    ids <- unique(precursor_id)
    if (length(ids) < 2L) {
      list(precursor_id = ids, coherence = NA_real_, n_siblings = 0L)
    } else {
      m <- feature_matrix(.SD, id_col = unit_col, feature_col = "precursor_id", value_col = "log2q")
      if (is.null(m) || ncol(m) < 2L) {
        list(precursor_id = ids, coherence = NA_real_, n_siblings = max(0L, length(ids) - 1L))
      } else {
        cc <- suppressWarnings(stats::cor(m, method = "pearson", use = "pairwise.complete.obs"))
        cn <- colnames(m)
        list(precursor_id = cn,
             coherence = vapply(cn, function(j) stats::median(cc[j, setdiff(cn, j)], na.rm = TRUE), numeric(1)),
             n_siblings = rep(ncol(m) - 1L, ncol(m)))
      }
    }
  }, by = .(pg_canonical, protease), .SDcols = c("precursor_id", unit_col, "log2q")]
}

# ---- peptide-set stability -------------------------------------------------------------------
# Median pairwise Jaccard of the observed peptide SETS across units. 1 = every unit sees the same
# peptides; low = the estimator's inputs change from unit to unit.
.peptide_set_jaccard <- function(sub, unit_col, min_units = 3L) {
  tab <- unique(sub[, .(u = as.character(get(unit_col)), p = stripped_sequence)])
  if (uniqueN(tab$u) < min_units || uniqueN(tab$p) < 2L) return(NA_real_)
  M <- table(tab$p, tab$u) > 0                      # peptides x units
  storage.mode(M) <- "numeric"
  inter <- crossprod(M)                             # units x units, |A n B|
  sz <- diag(inter)
  un <- outer(sz, sz, "+") - inter                  # |A u B|
  J <- inter / un
  stats::median(J[upper.tri(J)], na.rm = TRUE)
}

# ---- leave-one-peptide-out roll-up stability --------------------------------------------------
# msf : (feature x sample) log2 matrix, the exact input rollup_protein_tmp() consumes.
# Drops one feature at a time, re-polishes, and compares the per-sample estimate to the full one.
# Rank comparison is Spearman across samples (minimised over drops: the WEAKEST link is the risk);
# level comparison is the median |difference| after removing the constant offset a dropped feature
# induces in the median-polish overall term.
lopo_stability <- function(msf, max_features = 30L) {
  na <- list(lopo_rank_stability = NA_real_, lopo_max_delta = NA_real_, lopo_n_evaluated = 0L)
  if (is.null(msf) || nrow(msf) < 2L || ncol(msf) < 3L) return(na)
  full <- rollup_protein_tmp(msf)
  if (is.null(full)) return(na)
  th_full <- full$theta_hat
  n_obs_feat <- apply(msf, 1, function(v) sum(is.finite(v)))
  drop_idx <- order(-n_obs_feat)[seq_len(min(nrow(msf), max_features))]
  rr <- numeric(0); dd <- numeric(0)
  for (i in drop_idx) {
    r <- rollup_protein_tmp(msf[-i, , drop = FALSE])
    if (is.null(r)) next
    th <- r$theta_hat
    ok <- is.finite(th_full) & is.finite(th)
    if (sum(ok) >= 3L && stats::sd(th_full[ok]) > 0 && stats::sd(th[ok]) > 0)
      rr <- c(rr, suppressWarnings(stats::cor(th_full[ok], th[ok], method = "spearman")))
    if (sum(ok) >= 1L) {
      dlt <- th_full[ok] - th[ok]
      dd <- c(dd, stats::median(abs(dlt - stats::median(dlt))))
    }
  }
  list(lopo_rank_stability = if (length(rr)) min(rr, na.rm = TRUE) else NA_real_,
       lopo_max_delta      = if (length(dd)) max(dd, na.rm = TRUE) else NA_real_,
       lopo_n_evaluated    = length(drop_idx))
}

# ---- main -------------------------------------------------------------------------------------
# prec        : ALL ingested precursors, both proteases (pre-qualification)
# qual_ids    : precursor_ids that pass the FROZEN qualification for this condition
# collapsed   : collapse_peptide_hierarchy()$collapsed (frozen decisions applied)
# roll        : rollup_all() output, run grain
# cv_ref      : frozen CV reference (precursor_id, median_cv)
# icc_lookup  : frozen ICC lookup (pg_canonical, R_<labelA>, R_<labelB>)
# fa_index    : fasta_index() or NULL
# pgs         : proteins to evaluate (default: present in both channels)
# Returns one row per (pg_canonical, protease).
precursor_evidence <- function(prec, qual_ids, collapsed, roll, cv_ref = NULL, icc_lookup = NULL,
                               fa_index = NULL, pgs = NULL, labelA = "trypsin", labelB = "hta",
                               splithalf_reps = 200L, splithalf_min_paired = 6L,
                               lopo_max_features = 30L, jaccard_min_units = 3L, seed = 20260820) {
  units   <- sort(unique(prec$run))
  samples <- sort(unique(prec$sample_id))
  n_units <- length(units); n_samples <- length(samples)
  if (is.null(pgs)) {
    byp <- prec[, .(n_pro = uniqueN(protease)), by = pg_canonical]
    pgs <- byp[n_pro == 2L, pg_canonical]
  }
  P <- prec[pg_canonical %in% pgs & is.finite(log2q)]
  P[, is_qual := precursor_id %in% qual_ids]

  # --- coherence over ALL present precursors, per channel (and over the qualifying subset) ------
  coh_all  <- sibling_coherence_by_channel(P, unit_col = "run")
  coh_qual <- sibling_coherence_by_channel(P[is_qual == TRUE], unit_col = "run")
  setnames(coh_qual, c("coherence", "n_siblings"), c("coherence_qual", "n_siblings_qual"))

  # --- per-precursor detection within this condition -------------------------------------------
  pdet <- P[, .(prec_detection = uniqueN(run) / n_units), by = .(pg_canonical, protease, precursor_id)]

  # --- completeness (vectorised over all protein x protease) ------------------------------------
  cov_cols <- c("pg_canonical", "protease")
  comp <- P[, {
    q <- .SD[is_qual == TRUE]
    nq <- uniqueN(q$precursor_id)
    det_samp <- if (nrow(q)) q[, uniqueN(sample_id)] else 0L
    ge2 <- if (nrow(q)) q[, .(k = uniqueN(precursor_id)), by = sample_id][k >= 2L, .N] else 0L
    .(n_prec_present = uniqueN(precursor_id),
      n_prec_qual    = nq,
      n_peptides_present = uniqueN(stripped_sequence),
      n_peptides_qual = if (nrow(q)) uniqueN(q$stripped_sequence) else 0L,
      n_runs_obs = uniqueN(run),
      cell_fill_rate = if (nq > 0) nrow(unique(q[, .(precursor_id, run)])) / (nq * n_units) else NA_real_,
      protein_detection_rate = det_samp / n_samples,
      frac_samples_ge2_prec  = ge2 / n_samples,
      peptide_set_jaccard    = if (nrow(q)) .peptide_set_jaccard(q, "run", jaccard_min_units) else NA_real_)
  }, by = cov_cols, .SDcols = c("precursor_id","stripped_sequence","run","sample_id","is_qual","log2q")]

  med_pdet <- P[is_qual == TRUE, .(precursor_id = unique(precursor_id)), by = .(pg_canonical, protease)]
  med_pdet <- merge(med_pdet, pdet, by = c("pg_canonical","protease","precursor_id"))
  med_pdet <- med_pdet[, .(median_prec_detection = stats::median(prec_detection)), by = .(pg_canonical, protease)]

  coh_sum <- merge(coh_all, unique(P[, .(pg_canonical, protease, precursor_id, is_qual)]),
                   by = c("pg_canonical","protease","precursor_id"), all.x = TRUE)
  coh_sum <- merge(coh_sum, coh_qual, by = c("pg_canonical","protease","precursor_id"), all.x = TRUE)
  coh_agg <- coh_sum[, .(
    median_sibling_coherence = stats::median(coherence, na.rm = TRUE),
    frac_prec_coherent = mean(is.finite(coherence) & coherence >= 0.30, na.rm = TRUE),
    n_prec_with_siblings = sum(n_siblings > 0L),
    median_sibling_coherence_qual = if (any(is_qual)) stats::median(coherence_qual[is_qual], na.rm = TRUE) else NA_real_
  ), by = .(pg_canonical, protease)]

  # --- profile-grain metrics (split-half, LOPO, SE) --------------------------------------------
  set.seed(seed)
  C <- collapsed[pg_canonical %in% pgs]
  keys <- unique(C[, .(pg_canonical, protease)])
  prof <- vector("list", nrow(keys))
  for (i in seq_len(nrow(keys))) {
    pg <- keys$pg_canonical[i]; pro <- keys$protease[i]
    sub <- C[pg_canonical == pg & protease == pro]
    m_samp <- feature_matrix(sub, id_col = "sample_id", feature_col = "feature_id", value_col = "log2q")
    nf <- if (is.null(m_samp)) 0L else ncol(m_samp)
    sh <- if (nf >= 2L) channel_reliability(m_samp, n_splits = splithalf_reps, seed = seed,
                                            min_paired = splithalf_min_paired) else NA_real_
    lp <- lopo_stability(if (is.null(m_samp)) NULL else t(m_samp), max_features = lopo_max_features)
    prof[[i]] <- data.table(pg_canonical = pg, protease = pro, n_features = nf,
                            splithalf_reliability = sh,
                            lopo_rank_stability = lp$lopo_rank_stability,
                            lopo_max_delta = lp$lopo_max_delta,
                            lopo_n_evaluated = lp$lopo_n_evaluated)
  }
  prof <- rbindlist(prof, fill = TRUE)

  se_agg <- roll[pg_canonical %in% pgs & is.finite(se_hat),
                 .(median_rollup_se = stats::median(se_hat)), by = .(pg_canonical, protease)]
  est_agg <- unique(roll[pg_canonical %in% pgs, .(pg_canonical, protease, validity_flag, estimator)])

  # --- frozen technical CV + replicate ICC ------------------------------------------------------
  cvagg <- NULL
  if (!is.null(cv_ref) && nrow(cv_ref)) {
    q <- unique(P[is_qual == TRUE, .(pg_canonical, protease, precursor_id)])
    q <- merge(q, cv_ref[, .(precursor_id, median_cv)], by = "precursor_id", all.x = TRUE)
    cvagg <- q[, .(frozen_technical_cv = stats::median(median_cv, na.rm = TRUE),
                   n_prec_with_frozen_cv = sum(is.finite(median_cv))), by = .(pg_canonical, protease)]
  }
  iccagg <- NULL
  if (!is.null(icc_lookup) && nrow(icc_lookup)) {
    ra <- paste0("R_", labelA); rb <- paste0("R_", labelB)
    if (all(c(ra, rb) %in% names(icc_lookup)))
      iccagg <- rbind(icc_lookup[, .(pg_canonical, protease = labelA, icc_repro = get(ra))],
                      icc_lookup[, .(pg_canonical, protease = labelB, icc_repro = get(rb))])
  }

  # --- sequence coverage ------------------------------------------------------------------------
  covagg <- NULL
  if (!is.null(fa_index)) {
    qsub <- P[is_qual == TRUE]
    kk <- unique(qsub[, .(pg_canonical, protease)])
    rows <- vector("list", nrow(kk))
    for (i in seq_len(nrow(kk))) {
      pg <- kk$pg_canonical[i]; pro <- kk$protease[i]
      s <- qsub[pg_canonical == pg & protease == pro]
      byunit <- split(s$stripped_sequence, s$run)
      cv <- protein_coverage(pg, s$stripped_sequence, byunit, fa_index)
      rows[[i]] <- data.table(pg_canonical = pg, protease = pro,
                              seq_coverage_frac = cv$seq_coverage_frac,
                              median_unit_coverage_frac = cv$median_unit_coverage_frac,
                              protein_len = cv$protein_len, n_peptides_mapped = cv$n_peptides_mapped)
    }
    covagg <- rbindlist(rows, fill = TRUE)
  }

  out <- comp
  for (tb in list(med_pdet, coh_agg, prof, se_agg, est_agg, cvagg, iccagg, covagg))
    if (!is.null(tb) && nrow(tb)) out <- merge(out, tb, by = c("pg_canonical","protease"), all.x = TRUE)
  out[is.na(n_features), n_features := 0L]
  out[, `:=`(n_samples = n_samples, n_units = n_units)]
  setorder(out, pg_canonical, protease)
  out[]
}
