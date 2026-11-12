# reproducibility.R — Stage 6: how many proteins are WELL QUANTIFIED, per protease
# and combined, before and after cross-protease imputation.
#
# This is the stage the app was missing. It brings the Sep-9 manuscript framing
# (corrected inter-day CV, an explicit completeness definition, imputation tiers,
# and the union decomposition) inside the tool, and it PICKS THE RIGHT METRIC for
# the data it is given instead of assuming one:
#
#   * replicate + day structure present  ->  PRECISION axis. Inter-day CV, but
#     computed the defensible way: between-day SD in log2 from ALL points via a
#     variance-component estimator (injection noise removed), not from 3 day-means.
#     The published 3-day-mean CV is also reported, for continuity with the published
#     panel. Completeness = detected on every expected day.
#   * replicates but no day map          ->  PRECISION axis, coarse. Total
#     between-replicate log2 SD across the replicates of a sample.
#   * no replicates (discovery cohort)   ->  DETECTION axis. Cohort completeness:
#     the fraction of runs in which the protein is quantified. Inter-day CV is
#     undefined here and is NOT reported.
#
# In every case: Combined = UNION (a protein counts if EITHER protease quantifies
# it well). The two proteases are never averaged. Counts are reported per protease
# and as the union, measured-only and after imputation, with the decomposition
#   best single protease -> union (+X) -> + imputation (+Y)
# so the reader sees that the second protease, not the fill, carries the benefit.
#
# Depends on: data.table. Consumes pg_run (run grain: pg_canonical, protease, run,
# sample_id, log2q) and the imputation ledger (imp$matrix, if Stage 5 ran).

suppressPackageStartupMessages({ library(data.table) })

# ---- metric: precision (inter-day CV, corrected) ----------------------------
# moment-estimator variance components on log2 values, per protein x sample:
#   sigma2_day = var(day means) - MS_within / n_per_day   (floored at 0)
# reported as a between-day CV so it sits on the same axis as the published number.
.corrected_interday <- function(D) {   # D: pg_canonical, sample_id, day, rep, log2q (finite)
  dm <- D[, .(day_mean = mean(log2q), n = .N,
              ss = sum((log2q - mean(log2q))^2)), by = .(pg_canonical, sample_id, day)]
  agg <- dm[, {
    nd <- .N
    within_df <- sum(n) - nd
    ms_within <- if (within_df > 0) sum(ss) / within_df else NA_real_
    npd <- mean(n)
    s2_between_raw <- if (nd >= 2) stats::var(day_mean) else NA_real_
    s2_day <- if (is.finite(s2_between_raw) && is.finite(ms_within))
      max(s2_between_raw - ms_within / npd, 0) else s2_between_raw
    .(MeasuredDays = nd, sd_day_log2 = sqrt(s2_day),
      daymean_cv = if (nd >= 2) 100 * stats::sd(day_mean) / mean(day_mean) else NA_real_,
      # the published CV is on LINEAR day-means; recompute on linear scale to match exactly
      grand = mean(day_mean))
  }, by = .(pg_canonical, sample_id)]
  agg[]
}

# The published number is a LINEAR-scale CV of the 3 day-means. Compute it
# directly on linear intensities so the continuity check is exact.
.daymean_interday <- function(Dlin) {   # Dlin: pg_canonical, sample_id, day, intensity (finite,>0)
  dm <- Dlin[, .(day_mean = mean(intensity)), by = .(pg_canonical, sample_id, day)]
  dm[, .(CV = if (.N < 2) NA_real_ else 100 * stats::sd(day_mean) / mean(day_mean),
         MeasuredDays = .N), by = .(pg_canonical, sample_id)]
}

# ---- the stage --------------------------------------------------------------
# pg_run  : pg_canonical, protease, run, sample_id, log2q  (measured, both channels)
# imp_mat : imputed cells from Stage 5 (pg_canonical, sample_id, recipient, imputed_value)
#           or NULL. Imputation in this app fills the RECOMMENDED channel only; the
#           per-protease "after" count credits a protease only where it received a fill.
# cfg_rep : list(run_day_regex, cv_threshold, completeness_min, completeness_grid)
reproducibility_stage <- function(pg_run, imp_mat, labelA, labelB, cfg_rep, have_reps) {
  stopifnot(all(c("pg_canonical","protease","run","sample_id","log2q") %in% names(pg_run)))
  cv_thr  <- cfg_rep$cv_threshold %||% 20
  comp_min <- cfg_rep$completeness_min %||% 0.50
  day_rx  <- cfg_rep$run_day_regex %||% ""
  X <- copy(pg_run)[is.finite(log2q)]

  # try to recover a day label from the run name
  has_day <- FALSE
  if (nzchar(day_rx)) {
    d <- sub(day_rx, "\\1", X$run)
    if (!all(d == X$run)) { X[, day := d]; has_day <- TRUE }
  }

  mode <- if (have_reps && has_day) "interday_cv" else if (have_reps) "replicate_sd" else "completeness"

  # ---- per (protein, sample, protease) quality + completeness ----------------
  n_runs <- X[, .(n_run_total = uniqueN(run)), by = .(sample_id, protease)]
  if (mode == "interday_cv") {
    Dlin <- X[, .(day = day, intensity = 2^log2q, log2q = log2q, rep = run), by = .(pg_canonical, sample_id, protease)]
    exp_days <- X[, .(ExpectedDays = uniqueN(day)), by = .(sample_id, protease)]
    q <- rbindlist(lapply(split(X, X$protease), function(sub) {
      cc <- .corrected_interday(sub[, .(pg_canonical, sample_id, day, log2q)])
      ss <- .daymean_interday(sub[, .(pg_canonical, sample_id, day, intensity = 2^log2q)])
      m <- merge(cc, ss[, .(pg_canonical, sample_id, daymean_cv_linear = CV)], by = c("pg_canonical","sample_id"))
      m[, protease := sub$protease[1]][]
    }))
    q <- merge(q, exp_days, by = c("sample_id","protease"), all.x = TRUE)
    q[, complete := MeasuredDays == ExpectedDays]
    # corrected between-day CV (%) from the log2 SD, on the same footing as the published one
    q[, corrected_cv := 100 * (2^sd_day_log2 - 1)]
    q[, well := complete %in% TRUE & MeasuredDays >= 2 & is.finite(corrected_cv) & corrected_cv <= cv_thr]
    q[, well_daymean := complete %in% TRUE & MeasuredDays >= 2 & is.finite(daymean_cv_linear) & daymean_cv_linear <= cv_thr]
  } else if (mode == "replicate_sd") {
    q <- X[, .(sd_rep_log2 = if (.N >= 2) stats::sd(log2q) else NA_real_, n = .N),
           by = .(pg_canonical, sample_id, protease)]
    q <- merge(q, n_runs, by = c("sample_id","protease"))
    q[, corrected_cv := 100 * (2^sd_rep_log2 - 1)]
    q[, complete := n == n_run_total]
    q[, well := complete %in% TRUE & is.finite(corrected_cv) & corrected_cv <= cv_thr]
  } else {  # completeness
    q <- X[, .(n = uniqueN(run)), by = .(pg_canonical, sample_id, protease)]
    q <- merge(q, n_runs, by = c("sample_id","protease"))
    q[, completeness := n / n_run_total]
    q[, well := completeness >= comp_min]
  }

  # ---- imputation "after" (app fills the recommended channel only) -----------
  after_extra <- data.table(pg_canonical = character(), sample_id = character(), protease = character())
  if (!is.null(imp_mat) && nrow(imp_mat)) {
    # an imputed cell adds an observation to the recipient protease of that protein
    after_extra <- unique(imp_mat[, .(pg_canonical, sample_id,
                                      protease = if ("recipient" %in% names(imp_mat)) recipient else labelA)])
  }

  # ---- well-quantified counts, per protease + UNION, at protein x sample grain -
  # The manuscript's Combined row is a protein x sample union: a (protein, sample)
  # counts if EITHER protease quantifies it well there. Report the protein x sample
  # totals (comparable to Panel A's per-sample counts), the per-sample vector, and
  # the distinct-protein rollup (how many proteins clear the bar in >=1 sample).
  ps <- function(flag) {
    w <- q[get(flag) == TRUE, .(pg_canonical, sample_id, protease)]
    A <- w[protease == labelA, .(pg_canonical, sample_id)]
    B <- w[protease == labelB, .(pg_canonical, sample_id)]
    U <- unique(rbind(A, B))
    per_samp <- U[, .N, by = sample_id]
    list(A_ps = nrow(A), B_ps = nrow(B), union_ps = nrow(U),
         best_single_ps = max(nrow(A), nrow(B)),
         A_prot = uniqueN(A$pg_canonical), B_prot = uniqueN(B$pg_canonical),
         union_prot = uniqueN(U$pg_canonical),
         per_sample = per_samp)
  }
  main_flag <- if (mode == "interday_cv") "well" else "well"   # corrected CV is the headline metric
  m <- ps(main_flag)
  counts <- data.table(
    metric = mode, cv_definition = if (mode == "interday_cv") "corrected_interday" else mode,
    A_label = labelA, B_label = labelB,
    A_wellq_ps = m$A_ps, B_wellq_ps = m$B_ps, union_wellq_ps = m$union_ps,
    best_single_ps = m$best_single_ps, union_gain_over_single_ps = m$union_ps - m$best_single_ps,
    A_wellq_proteins = m$A_prot, B_wellq_proteins = m$B_prot, union_wellq_proteins = m$union_prot,
    n_imputed_cells = if (is.null(imp_mat)) 0L else nrow(imp_mat))
  # continuity: the published day-mean CV, when available
  if (mode == "interday_cv") {
    s <- ps("well_daymean")
    counts[, `:=`(daymean_A_ps = s$A_ps, daymean_B_ps = s$B_ps, daymean_union_ps = s$union_ps)]
  }
  list(mode = mode, cv_threshold = cv_thr, completeness_min = comp_min,
       per_protein = q, per_sample = m$per_sample, counts = counts)
}

# short human-readable summary line for the methods doc / stage report
reproducibility_headline <- function(rs, labelA, labelB) {
  c <- rs$counts
  metric_txt <- switch(rs$mode,
    interday_cv  = sprintf("corrected inter-day CV <= %g%% and complete days", rs$cv_threshold),
    replicate_sd = sprintf("between-replicate CV <= %g%% and all replicates present", rs$cv_threshold),
    completeness = sprintf(">= %g%% run completeness", 100 * rs$completeness_min))
  sprintf(paste0("Well quantified (%s), protein x sample: %s %d | %s %d | union %d (+%d over the best ",
                 "single protease). The two proteases are never averaged; the union is the second ",
                 "protease's contribution, and imputation is a refinement on top of it."),
          metric_txt, labelA, c$A_wellq_ps, labelB, c$B_wellq_ps, c$union_wellq_ps,
          c$union_gain_over_single_ps)
}
