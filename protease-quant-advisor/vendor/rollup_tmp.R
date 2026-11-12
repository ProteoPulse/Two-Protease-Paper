# rollup_tmp.R — Stage 2: within-protease protein roll-up from QUALIFYING precursors.
#
# This is the headline new capability of the Aug-18 plan: instead of consuming DIA-NN's MaxLFQ pg
# quantity as-is, we build our own per-sample protein estimate from the precursors we have
# qualified ourselves (Stage 1). Two steps:
#
#   (a) collapse_peptide_hierarchy(): group qualifying precursors by stripped sequence. Charge/mod
#       variants of one peptide that agree (median pairwise Pearson >= corr_min across samples) are
#       collapsed to a single representative (the best-detected, lowest-CV variant); variants that
#       DISAGREE are kept as separate features (they may report different proteoforms). Decisions are
#       LEARNED on reproducibility and re-APPLIED verbatim to discovery (no re-learning downstream).
#   (b) rollup_protein_tmp(): Tukey median polish (stats::medpolish) on the log2 (feature x sample)
#       matrix. The per-sample protein estimate is overall + column effect; residuals give a robust
#       per-sample SE. Median polish is resistant to the missing/ragged features HTA produces; when a
#       matrix is too sparse for a stable polish it falls back to a column-median summary (recorded).
#
# Output is a long roll-up table that rollup_to_channel_long() splits by protease into the exact
# A_long / B_long contract the rest of the engine speaks, plus SE and validity side-tables.
#
# Depends on: data.table, feature_matrix.R (feature_matrix, median_pairwise_cor), utils_ids.R.

suppressPackageStartupMessages({ library(data.table) })

# ---- (a) peptide-hierarchy collapsing --------------------------------------------------------
# prec_long : qualifying precursors, long. Columns: pg_canonical, protease, stripped_sequence,
#             precursor_id, sample_id, log2q. Optional per-precursor metrics for "best" selection:
#             detection_rate, median_cv (joined by precursor_id) — else best = most-observed.
# decisions : NULL (learn) or a prior data.table (apply). When applying, reuse the kept representative
#             per (protein,protease,stripped); if it is absent here, fall back to the next-best
#             present variant (deterministic; recorded as fallback=TRUE).
# Returns list(collapsed, decisions):
#   collapsed : pg_canonical, protease, feature_id, sample_id, log2q   (feature_id = stripped when
#               collapsed, else precursor_id when kept separate)
#   decisions : pg_canonical, protease, stripped_sequence, action(collapse|separate|singleton),
#               kept_precursor_id, variant_corr, n_variants, fallback
collapse_peptide_hierarchy <- function(prec_long, decisions = NULL, corr_min = 0.70) {
  d <- prec_long[is.finite(log2q)]
  if (!"run" %in% names(d)) d[, run := sample_id]            # discovery: one run per sample
  # per-precursor selection stats (highest detection, then lowest CV) if not provided
  if (!"detection_rate" %in% names(d)) {
    det <- d[, .(detection_rate = uniqueN(run)), by = precursor_id]
    d <- merge(d, det, by = "precursor_id", all.x = TRUE)
  }
  if (!"median_cv" %in% names(d)) d[, median_cv := NA_real_]
  best_variant <- function(sub) {
    st <- unique(sub[, .(precursor_id, detection_rate, median_cv)])
    setorder(st, -detection_rate, median_cv, precursor_id)  # NA CV sorts last within a detection tie
    st$precursor_id[1]
  }
  keys <- unique(d[, .(pg_canonical, protease, stripped_sequence)])
  dec <- vector("list", nrow(keys)); coll <- vector("list", nrow(keys))
  learn <- is.null(decisions)
  for (i in seq_len(nrow(keys))) {
    pg <- keys$pg_canonical[i]; pro <- keys$protease[i]; sq <- keys$stripped_sequence[i]
    sub <- d[pg_canonical == pg & protease == pro & stripped_sequence == sq]
    vars <- unique(sub$precursor_id); nv <- length(vars)
    fallback <- FALSE
    if (nv == 1L) {
      action <- "singleton"; kept <- vars[1]; vc <- NA_real_
    } else {
      m <- feature_matrix(sub, id_col = "run", feature_col = "precursor_id", value_col = "log2q")
      vc <- median_pairwise_cor(m, method = "pearson")
      prior <- if (!learn) decisions[pg_canonical == pg & protease == pro & stripped_sequence == sq] else NULL
      if (!is.null(prior) && nrow(prior)) {
        action <- prior$action[1]; kept <- prior$kept_precursor_id[1]
        if (action == "collapse" && !(kept %in% vars)) { kept <- best_variant(sub); fallback <- TRUE }
      } else {
        action <- if (is.finite(vc) && vc >= corr_min) "collapse" else "separate"
        kept <- if (action == "collapse") best_variant(sub) else NA_character_
      }
    }
    dec[[i]] <- data.table(pg_canonical = pg, protease = pro, stripped_sequence = sq,
                           action = action, kept_precursor_id = kept, variant_corr = vc,
                           n_variants = nv, fallback = fallback)
    if (action == "separate") {
      coll[[i]] <- sub[, .(pg_canonical, protease, feature_id = precursor_id, sample_id, run, log2q)]
    } else {
      rep_id <- if (action == "singleton") vars[1] else kept
      coll[[i]] <- sub[precursor_id == rep_id,
                       .(pg_canonical, protease, feature_id = sq, sample_id, run, log2q)]
    }
  }
  list(collapsed = rbindlist(coll, fill = TRUE), decisions = rbindlist(dec, fill = TRUE))
}

# ---- (b) Tukey median polish roll-up for ONE protein x protease ------------------------------
# m : (feature x sample) log2 matrix (rows=features, cols=samples). Returns a data.table with one
# row per sample: sample_id, theta_hat, se_hat, n_obs, plus attr "estimator" (medpolish|median).
rollup_protein_tmp <- function(m, min_features_for_tmp = 2L, max_iter = 20L) {
  if (is.null(m) || !nrow(m) || !ncol(m)) return(NULL)
  samples <- colnames(m)
  n_obs <- apply(m, 2, function(v) sum(is.finite(v)))
  estimator <- "medpolish"
  if (nrow(m) < min_features_for_tmp) {
    # single feature (or none): the estimate IS that feature; no residual spread within a column
    theta <- apply(m, 2, function(v) { v <- v[is.finite(v)]; if (!length(v)) NA_real_ else stats::median(v) })
    se <- rep(NA_real_, length(samples)); estimator <- "median"
  } else {
    mp <- tryCatch(stats::medpolish(m, eps = 1e-4, maxiter = max_iter, trace.iter = FALSE, na.rm = TRUE),
                   error = function(e) NULL)
    if (is.null(mp)) {
      theta <- apply(m, 2, function(v) { v <- v[is.finite(v)]; if (!length(v)) NA_real_ else stats::median(v) })
      se <- rep(NA_real_, length(samples)); estimator <- "median"
    } else {
      theta <- mp$overall + mp$col                    # per-sample protein estimate (relative log2)
      # robust per-sample SE from the column's residuals: MAD-scaled, shrunk by sqrt(#features seen)
      res <- mp$residuals
      se <- vapply(seq_along(samples), function(j) {
        r <- res[, j]; r <- r[is.finite(r)]
        if (length(r) < 1) return(NA_real_)
        (stats::median(abs(r)) * 1.4826) / sqrt(max(1, length(r)))
      }, numeric(1))
      theta[n_obs == 0] <- NA_real_
    }
  }
  out <- data.table(sample_id = samples, theta_hat = as.numeric(theta),
                    se_hat = as.numeric(se), n_obs = as.integer(n_obs))
  setattr(out, "estimator", estimator)
  out
}

# ---- validity flag per protein x protease ----------------------------------------------------
validity_flag <- function(n_features, n_samples_detected, n_samples_total, detection_min = 0.50) {
  if (n_features == 0L) return("absent")
  det <- if (n_samples_total > 0) n_samples_detected / n_samples_total else 0
  if (det < detection_min) return("insufficient")
  if (n_features == 1L) return("single_peptide")
  "valid"
}

# ---- roll up ALL proteins x proteases --------------------------------------------------------
# collapsed : output of collapse_peptide_hierarchy()$collapsed (carries run + sample_id)
# Roll-up is done at the finest UNIT (run) so the replicate structure survives for ICC (Stage 3A)
# and per-run performance (Stage 4). Detection-based validity is judged at the SAMPLE level.
# Returns long at run grain: pg_canonical, protease, run, sample_id, theta_hat, se_hat, n_obs,
#   validity_flag, n_features, estimator.
rollup_all <- function(collapsed, n_samples_total = NULL, detection_min = 0.50,
                       min_features_for_tmp = 2L, unit_col = "run") {
  if (!unit_col %in% names(collapsed)) collapsed[, (unit_col) := sample_id]
  if (is.null(n_samples_total))
    n_samples_total <- collapsed[, uniqueN(sample_id)]
  unit_map <- unique(collapsed[, .(unit = get(unit_col), sample_id)])
  keys <- unique(collapsed[, .(pg_canonical, protease)])
  out <- vector("list", nrow(keys))
  for (i in seq_len(nrow(keys))) {
    pg <- keys$pg_canonical[i]; pro <- keys$protease[i]
    sub <- collapsed[pg_canonical == pg & protease == pro]
    n_features <- sub[, uniqueN(feature_id)]
    m <- feature_matrix(sub, id_col = unit_col, feature_col = "feature_id", value_col = "log2q")
    mt <- if (!is.null(m)) t(m) else NULL          # feature_matrix gives unit x feature; TMP wants feature x unit
    rt <- rollup_protein_tmp(mt, min_features_for_tmp = min_features_for_tmp)
    if (is.null(rt)) next
    est <- attr(rt, "estimator")
    setnames(rt, "sample_id", "unit")              # rollup_protein_tmp labels columns as "sample_id"
    rt <- merge(rt, unit_map, by = "unit", all.x = TRUE)
    n_samp_det <- rt[is.finite(theta_hat), uniqueN(sample_id)]
    vflag <- validity_flag(n_features, n_samp_det, n_samples_total, detection_min)
    rt[, `:=`(pg_canonical = pg, protease = pro, validity_flag = vflag,
              n_features = n_features, estimator = est)]
    out[[i]] <- rt
  }
  res <- rbindlist(out, fill = TRUE)
  setnames(res, "unit", "run")
  setcolorder(res, c("pg_canonical","protease","run","sample_id","theta_hat","se_hat","n_obs",
                     "validity_flag","n_features","estimator"))
  res[]
}

# ---- ADAPTER: roll-up long -> A_long / B_long (+ SE + validity side-tables) -------------------
# This is the seam that lets everything downstream (compare_metrics, classify_multiaxis, ...) run
# unchanged on custom-TMP estimates. A_long/B_long match harmonize_to_long()'s schema exactly.
#   roll   : rollup_all() output
#   annA/annB : optional pg_annotation() tables to recover protein_group_original (else pg_canonical)
# Returns a list mimicking harmonize_two(): A_long, B_long, annA, annB, shared_samples, n_samples,
#   shared_pg, disposition, labelA, labelB, PLUS se_side and validity (custom-TMP extras).
# Roll-up is at run grain; concordance is per biological sample, so run estimates are aggregated to
# a per-sample estimate (median over replicate runs — identity for discovery's one-run-per-sample).
rollup_to_channel_long <- function(roll, labelA = "trypsin", labelB = "hta",
                                   annA = NULL, annB = NULL) {
  mk_long <- function(pro, ann) {
    d <- roll[protease == pro & is.finite(theta_hat) & validity_flag != "absent",
              .(log2q = stats::median(theta_hat)), by = .(pg_canonical, sample_id)]
    if (!is.null(ann)) {
      d <- merge(d, unique(ann[, .(pg_canonical, protein_group_original)]),
                 by = "pg_canonical", all.x = TRUE)
      d[is.na(protein_group_original), protein_group_original := pg_canonical]
    } else d[, protein_group_original := pg_canonical]
    d[, .(pg_canonical, protein_group_original, sample_id, log2q)]
  }
  A_long <- mk_long(labelA, annA); B_long <- mk_long(labelB, annB)
  shared_samples <- sort(unique(c(A_long$sample_id, B_long$sample_id)))
  a_det <- unique(A_long$pg_canonical); b_det <- unique(B_long$pg_canonical)
  shared_pg <- intersect(a_det, b_det)
  disposition <- data.table(
    class = c(paste0(labelA, "_only"), "shared", paste0(labelB, "_only")),
    n = c(length(setdiff(a_det, b_det)), length(shared_pg), length(setdiff(b_det, a_det))))
  # per-sample SE: combine replicate-run SEs (mean of variances / n_runs), fallback median se_hat
  se_side <- roll[is.finite(theta_hat), {
    v <- se_hat[is.finite(se_hat)]
    .(se_hat = if (length(v)) sqrt(mean(v^2)) / sqrt(uniqueN(run)) else NA_real_,
      n_runs = uniqueN(run))
  }, by = .(pg_canonical, protease, sample_id)]
  validity <- unique(roll[, .(pg_canonical, protease, validity_flag, n_features, estimator)])
  list(A_long = A_long, B_long = B_long, annA = annA, annB = annB,
       shared_samples = shared_samples, n_samples = length(shared_samples),
       shared_pg = shared_pg, disposition = disposition, labelA = labelA, labelB = labelB,
       se_side = se_side, validity = validity)
}
