# validate_recommendation.R — does the recommendation predict anything, or is it just a score?
#
# A weighted score always produces a winner. The question is whether that winner is the channel whose
# roll-up is genuinely the better measurement. The reproducibility data can answer it directly: it has
# technical replicates, so each protein x protease has a MEASURED technical CV — an outcome the score
# is supposed to anticipate.
#
# The test only means something if the predictor cannot see the answer. So the score is rebuilt with
# the technical_cv and replicate_icc criteria REMOVED (config recommend.holdout_drop) and asked to
# predict the lower-CV channel from completeness/coherence/stability alone. Reported against the 50%
# coin flip with an exact binomial p. Ties in measured CV are excluded (unpredictable by construction),
# and the count of excluded proteins is reported so the denominator is visible.
#
# Failing this test is informative, not fatal — it would say the score's usable criteria do not track
# measurement quality, and the weights need rework BEFORE the discovery runs are worth doing.
#
# Depends on: data.table.

suppressPackageStartupMessages({ library(data.table) })

# Measured within-sample technical CV (linear scale) of the rolled-up estimate, per protein x
# protease: CV across replicate runs within a biological sample, median over samples. NA without
# replicates.
measured_technical_cv <- function(roll) {
  d <- roll[is.finite(theta_hat)]
  per <- d[, {
    lin <- 2^theta_hat
    .(cv = if (.N >= 2L && mean(lin) > 0) stats::sd(lin) / mean(lin) else NA_real_)
  }, by = .(pg_canonical, protease, sample_id)]
  per[, .(technical_cv_measured = stats::median(cv, na.rm = TRUE),
          n_samples_with_reps = sum(is.finite(cv))), by = .(pg_canonical, protease)]
}

# Held-out test: does the CV/ICC-free score pick the lower-measured-CV channel?
holdout_cv_test <- function(rec_holdout, cvm, labelA = "trypsin", labelB = "hta") {
  cw <- dcast(cvm, pg_canonical ~ protease, value.var = "technical_cv_measured")
  if (!all(c(labelA, labelB) %in% names(cw)))
    return(list(summary = data.table(n_testable = 0L), detail = data.table()))
  setnames(cw, c(labelA, labelB), c("cv_a", "cv_b"))
  d <- merge(rec_holdout[, .(pg_canonical, recommended_channel, recommendation_margin,
                             recommendation_confidence, n_usable_criteria)],
             cw[, .(pg_canonical, cv_a, cv_b)], by = "pg_canonical")
  d[, both_cv := is.finite(cv_a) & is.finite(cv_b)]
  d[, truth := fifelse(!both_cv, NA_character_,
                fifelse(cv_a < cv_b, labelA, fifelse(cv_b < cv_a, labelB, NA_character_)))]
  t <- d[!is.na(truth)]
  hits <- t[recommended_channel == truth, .N]; n <- nrow(t)
  bt <- if (n > 0) stats::binom.test(hits, n, 0.5) else NULL
  by_conf <- t[, .(n = .N, hit_rate = mean(recommended_channel == truth)),
               by = recommendation_confidence]
  summ <- data.table(
    n_shared_scored = nrow(d), n_excluded_no_cv = d[both_cv == FALSE, .N],
    n_excluded_cv_tie = d[both_cv == TRUE & is.na(truth), .N], n_testable = n,
    n_hits = hits, hit_rate = if (n) hits / n else NA_real_,
    binom_p = if (!is.null(bt)) bt$p.value else NA_real_,
    ci_low = if (!is.null(bt)) bt$conf.int[1] else NA_real_,
    ci_high = if (!is.null(bt)) bt$conf.int[2] else NA_real_)
  list(summary = summ, by_confidence = by_conf, detail = d)
}

# Do two conditions recommend the same channel for the same protein?
cross_condition_stability <- function(rec1, rec2, name1 = "mix30", name2 = "mix20") {
  a <- rec1[, .(pg_canonical, rec_1 = recommended_channel, conf_1 = recommendation_confidence,
                margin_1 = recommendation_margin)]
  b <- rec2[, .(pg_canonical, rec_2 = recommended_channel, conf_2 = recommendation_confidence,
                margin_2 = recommendation_margin)]
  d <- merge(a, b, by = "pg_canonical")
  setnames(d, c("rec_1","conf_1","margin_1","rec_2","conf_2","margin_2"),
           c(paste0("rec_", name1), paste0("conf_", name1), paste0("margin_", name1),
             paste0("rec_", name2), paste0("conf_", name2), paste0("margin_", name2)))
  d[, agree := get(paste0("rec_", name1)) == get(paste0("rec_", name2))]
  d[]
}
