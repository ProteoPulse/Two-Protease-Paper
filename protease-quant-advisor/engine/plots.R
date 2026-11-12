# plots.R — every figure the app and the report show. Each returns a ggplot so the GUI can render it
# live and the HTML report can embed the identical object.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
TRY_COL <- "#2C6FB3"; HTA_COL <- "#E08214"; OK_COL <- "#1B7837"; BAD_COL <- "#D6604D"; GREY <- "#8A94A2"

theme_adv <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(colour = "grey90"),
          plot.title = element_text(face = "bold", size = base + 1),
          plot.subtitle = element_text(colour = "grey35", size = base - 1),
          strip.text = element_text(face = "bold"), legend.position = "bottom")
}
.chan_scale <- function(res, aes_fn = scale_fill_manual)
  aes_fn(values = setNames(c(TRY_COL, HTA_COL), c(res$cfg$label_a, res$cfg$label_b)), name = NULL)

# 1. What the two searches see, and how much of it needs a decision at all.
p_detection <- function(res) {
  d <- copy(res$detection)
  d[, category := factor(category, levels = rev(category))]
  ggplot(d, aes(n_proteins, category, fill = grepl("shared", category))) +
    geom_col(width = .62) +
    geom_text(aes(label = n_proteins), hjust = -.15, size = 3.4) +
    scale_fill_manual(values = c(`TRUE` = OK_COL, `FALSE` = GREY), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, .16))) +
    labs(title = "Only the shared proteins need adjudicating",
         subtitle = "Proteins seen by one search are used as they are; there is no competing value to choose between.",
         x = "proteins", y = NULL) + theme_adv()
}

# 2. Stage 1 — does the qualification rule hold quality without throwing the data away?
p_qual_sweep <- function(res) {
  s <- res$qual_sweep; if (is.null(s) || !nrow(s)) return(NULL)
  q <- res$cfg$qualify
  ggplot(s, aes(factor(detection_rate_min), factor(coherence_min), fill = n_precursors_pass)) +
    geom_tile(colour = "white", linewidth = .6) +
    geom_text(aes(label = scales::comma(n_precursors_pass)), size = 3, colour = "white") +
    annotate("point", x = factor(q$detection_rate_min), y = factor(q$coherence_min),
             shape = 21, size = 13, stroke = 1.4, colour = "black", fill = NA) +
    scale_fill_gradient(low = "#9FB6CC", high = "#1F4E79", guide = "none") +
    labs(title = "Stage 1 — qualification sensitivity",
         subtitle = "Qualifying precursors across the threshold grid. The ring marks the setting in use.",
         x = "minimum detection rate", y = "minimum sibling coherence") + theme_adv()
}

# 3. Stage 2 — the evidence that separates the two channels.
p_evidence <- function(res) {
  keep <- c(completeness = "protein_detection_rate", depth = "n_features",
            coherence = "median_sibling_coherence", `peptide-set stability` = "peptide_set_jaccard",
            `self-consistency` = "splithalf_reliability")
  e <- res$evidence
  rows <- rbindlist(lapply(names(keep), function(nm) {
    cn <- keep[[nm]]; if (!cn %in% names(e)) return(NULL)
    data.table(metric = nm, protease = e$protease, value = e[[cn]])
  }), fill = TRUE)
  rows <- rows[is.finite(value)]
  rows[metric == "depth", value := log10(pmax(value, 1))]
  rows[, metric := factor(metric, levels = names(keep))]
  ggplot(rows, aes(value, metric, fill = protease)) +
    ggplot2::geom_violin(colour = NA, alpha = .65, scale = "width", position = position_dodge(.75)) +
    stat_summary(fun = median, geom = "point", position = position_dodge(.75), size = 1.9) +
    .chan_scale(res) +
    labs(title = "Stage 2 — what separates the two channels",
         subtitle = "Per protein, per protease, over the shared set. Depth is log10(features).",
         x = "value", y = NULL) + theme_adv()
}

# 4. Stage 4 — where the agreement threshold comes from.
p_concordance <- function(res) {
  obs <- res$concordance[is.finite(stat_value), .(stat_value)]
  nul <- data.table(stat_value = as.numeric(res$null_rhos))
  d <- rbindlist(list(cbind(obs, src = "observed proteins"), cbind(nul, src = "permutation null")))
  ggplot(d, aes(stat_value, fill = src)) +
    geom_histogram(aes(y = after_stat(density)), bins = 46, colour = "white",
                   linewidth = .18, alpha = .72, position = "identity") +
    geom_vline(xintercept = res$threshold, linetype = "dashed", linewidth = .7) +
    annotate("text", x = res$threshold, y = Inf, vjust = 1.6, hjust = -.08, size = 3.2,
             label = sprintf("threshold %.3f", res$threshold)) +
    scale_fill_manual(values = c(`observed proteins` = TRY_COL, `permutation null` = GREY), name = NULL) +
    labs(title = "Stage 4 — do the two channels agree about this protein?",
         subtitle = paste0("The threshold is the ", 100*as.numeric(res$cfg$concordance$null_percentile),
                           "th percentile of a permutation null, not a chosen number."),
         x = paste0(res$cfg$concordance$method, " statistic"), y = "density") + theme_adv()
}

# 5. Stage 3 — the calls, and how confident they are.
p_recommendation <- function(res) {
  r <- res$rec[!is.na(recommended_channel)]
  d <- r[, .N, by = .(recommended_channel, recommendation_confidence)]
  d[, recommendation_confidence := factor(recommendation_confidence, c("high","moderate","low"))]
  ggplot(d, aes(N, recommended_channel, fill = recommendation_confidence)) +
    geom_col(width = .6) +
    geom_text(aes(label = N), position = position_stack(vjust = .5), size = 3.2, colour = "white") +
    scale_fill_manual(values = c(high = "#1F4E79", moderate = "#5E8CB8", low = "#B7C7D6"),
                      name = "confidence") +
    labs(title = "Stage 3 — which protease to trust, per shared protein",
         subtitle = "Every shared protein gets a call. Thin evidence shows up as low confidence, never as an abstention.",
         x = "proteins", y = NULL) + theme_adv()
}

# 6. Which calls actually change the answer.
p_consequential <- function(res) {
  r <- copy(res$rec)
  r[, grp := fifelse(is.na(consequential), "not classified",
              fifelse(consequential, "channels disagree — the choice changes the answer",
                      "channels agree — the choice is near-cosmetic"))]
  d <- r[, .N, by = grp][order(-N)]
  d[, grp := factor(grp, levels = rev(grp))]
  ggplot(d, aes(N, grp, fill = grepl("disagree", grp))) +
    geom_col(width = .58) + geom_text(aes(label = N), hjust = -.15, size = 3.4) +
    scale_fill_manual(values = c(`TRUE` = BAD_COL, `FALSE` = GREY), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, .17))) +
    labs(title = "How many recommendations actually matter?",
         subtitle = "Where the two channels already agree, picking one over the other barely moves the number.",
         x = "shared proteins", y = NULL) + theme_adv()
}

# 7. Stage 5 — what the gap filling did, and what it refused.
p_imputation <- function(res) {
  im <- res$imputation; if (is.null(im)) return(NULL)
  d <- im$refusal_breakdown
  if (!nrow(d) && !nrow(im$ledger)) return(NULL)
  d <- rbindlist(list(
    if (nrow(im$ledger)) data.table(refusal_reason = "FILLED from the other protease", N = nrow(im$ledger)),
    d), fill = TRUE)
  d[, refusal_reason := factor(refusal_reason, levels = rev(refusal_reason))]
  ggplot(d, aes(N, refusal_reason, fill = grepl("^FILLED", refusal_reason))) +
    geom_col(width = .6) + geom_text(aes(label = N), hjust = -.15, size = 3.3) +
    scale_fill_manual(values = c(`TRUE` = OK_COL, `FALSE` = GREY), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, .19))) +
    labs(title = "Stage 5 — every gap, filled or refused, with the reason",
         subtitle = "On well-measured data most gaps are refused, and that is the intended behaviour: a donor that disagrees is not a surrogate.",
         x = "cells", y = NULL) + theme_adv()
}

# 8. The probe's diagnostic — proteins whose value rests on one peptide.
p_probe <- function(res) {
  pf <- res$probe_flags; if (is.null(pf) || !nrow(pf)) return(NULL)
  ggplot(pf[is.finite(lopo_rank_stability)],
         aes(pmax(n_features, 1), lopo_rank_stability, colour = fragile_rollup)) +
    geom_point(alpha = .55, size = 1.5) + scale_x_log10() +
    facet_wrap(~ protease) +
    scale_colour_manual(values = c(`FALSE` = GREY, `TRUE` = BAD_COL),
                        labels = c("ordinary", "value rests on one peptide"), name = NULL) +
    labs(title = "Diagnostic — how much does the protein depend on a single peptide?",
         subtitle = "Leave-one-peptide-out stability. Reported per protein; deliberately NOT counted as evidence in the score.",
         x = "qualifying features (log scale)", y = "LOPO rank stability") + theme_adv()
}

# ---------------------------------------------------------------------------------------------
# Stage 1 — representative precursor QC. The sensitivity heatmap above shows that the thresholds
# are not knife-edge; these two show what the thresholds are actually doing to real precursors.
# ---------------------------------------------------------------------------------------------

# 9. Where each qualification rule cuts, drawn on the distribution it cuts.
p_qc_metrics <- function(res) {
  m <- res$metrics; if (is.null(m) || !nrow(m)) return(NULL)
  q <- res$cfg$qualify
  spec <- list(
    list(col = "detection_rate", lab = "detection rate\n(fraction of runs seen in)", thr = as.numeric(q$detection_rate_min), side = "above"),
    list(col = "coherence",      lab = "sibling coherence\n(agreement with the protein's other precursors)", thr = as.numeric(q$coherence_min), side = "above"),
    list(col = "median_cv",      lab = "median CV\n(within-sample replicate scatter)", thr = as.numeric(q$median_cv_max), side = "below"))
  rows <- rbindlist(lapply(spec, function(s) {
    if (!s$col %in% names(m)) return(NULL)
    v <- m[[s$col]]; ok <- is.finite(v)
    if (!any(ok)) return(NULL)
    data.table(metric = s$lab, protease = m$protease[ok], value = v[ok], thr = s$thr, side = s$side)
  }), fill = TRUE)
  if (is.null(rows) || !nrow(rows)) return(NULL)
  rows[, metric := factor(metric, levels = unique(vapply(spec, `[[`, "", "lab")))]
  thr <- unique(rows[, .(metric, thr, side)])
  npass <- res$qual[qualifies == TRUE, .N]; ntot <- nrow(res$qual)
  ggplot(rows, aes(value, fill = protease)) +
    geom_histogram(bins = 44, colour = "white", linewidth = .15, alpha = .7, position = "identity") +
    geom_vline(data = thr, aes(xintercept = thr), linetype = "dashed", linewidth = .65) +
    geom_text(data = thr, aes(x = thr, y = Inf, label = sprintf("keep %s %.2f", ifelse(side == "above", "≥", "≤"), thr)),
              inherit.aes = FALSE, vjust = 1.5, hjust = -0.06, size = 2.9, colour = "grey20") +
    facet_wrap(~ metric, scales = "free", ncol = 3) +
    .chan_scale(res) +
    labs(title = "Stage 1 — what the qualification rules actually cut",
         subtitle = sprintf("Every precursor in both searches. %s of %s precursors (%.1f%%) qualify to carry evidence into Stage 2.",
                            format(npass, big.mark = ","), format(ntot, big.mark = ","), 100 * npass / max(ntot, 1)),
         x = NULL, y = "precursors") + theme_adv()
}

# 10. What a qualifying precursor looks like next to a failing one, on the same protein.
p_qc_traces <- function(res) {
  prec <- res$prec; qual <- res$qual
  if (is.null(prec) || !nrow(prec) || is.null(qual) || !nrow(qual)) return(NULL)
  qq <- unique(qual[, .(precursor_id, protease, qualifies, coherence, detection_rate)])
  # A protein that carries BOTH kinds of precursor makes the contrast honest: same protein, same
  # samples, same load — so any difference in the traces is the precursor's own behaviour.
  cand <- qual[, .(n_pass = sum(qualifies), n_fail = sum(!qualifies)), by = .(pg_canonical, protease)]
  cand <- cand[n_pass >= 2 & n_fail >= 2 & (n_pass + n_fail) <= 14][order(-pmin(n_pass, n_fail))]
  if (!nrow(cand)) return(NULL)
  pick <- cand[1]
  d <- merge(prec[pg_canonical == pick$pg_canonical & protease == pick$protease,
                  .(precursor_id, sample_id, log2q)],
             qq[, .(precursor_id, protease, qualifies, coherence)][protease == pick$protease],
             by = "precursor_id")
  if (!nrow(d)) return(NULL)
  ord <- d[, .(m = stats::median(log2q, na.rm = TRUE)), by = sample_id][order(m)]
  d[, sample_id := factor(sample_id, levels = ord$sample_id)]
  d[, status := fifelse(qualifies, "qualifies — tracks its siblings", "rejected — does not track its siblings")]
  gene <- unique(prec[pg_canonical == pick$pg_canonical, pg_canonical])[1]
  ggplot(d, aes(as.integer(sample_id), log2q, group = precursor_id, colour = status)) +
    geom_line(alpha = .8, linewidth = .5) + geom_point(size = .8, alpha = .7) +
    scale_colour_manual(values = c(OK_COL, BAD_COL) |> setNames(c("qualifies — tracks its siblings",
                                                                 "rejected — does not track its siblings")), name = NULL) +
    labs(title = sprintf("Stage 1 — a representative protein (%s, %s search)", gene, pick$protease),
         subtitle = paste0("Every precursor this search assigns to the protein, across samples ordered by the protein's level. ",
                           "Qualification keeps the precursors that move together; the rejected ones are what roll-up noise looks like."),
         x = "sample (ordered by protein level)", y = "log2 intensity") + theme_adv()
}

# ---------------------------------------------------------------------------------------------
# Stage 5 — what imputation did to each value it wrote. Dumbbells, ordered by intensity.
# ---------------------------------------------------------------------------------------------

# 11. Cell level: the donor's own value, and the value we actually write after translating it.
p_impute_dumbbell <- function(res) {
  im <- res$imputation; if (is.null(im) || is.null(im$ledger) || !nrow(im$ledger)) return(NULL)
  d <- copy(im$ledger)
  d[, key := paste(pg_canonical, sample_id)]
  setorder(d, imputed_value)
  d[, key := factor(key, levels = key)]
  d[, dir := fifelse(offset_used >= 0, "translated up", "translated down")]
  med <- stats::median(d$offset_used, na.rm = TRUE)
  ggplot(d, aes(y = key)) +
    geom_segment(aes(x = theta_don, xend = imputed_value, yend = key, colour = dir), linewidth = .55, alpha = .85) +
    geom_point(aes(x = theta_don), colour = GREY, size = 1.5) +
    geom_point(aes(x = imputed_value, colour = dir), size = 1.7) +
    scale_colour_manual(values = c(`translated up` = TRY_COL, `translated down` = HTA_COL), name = NULL) +
    labs(title = "Stage 5 — every imputed cell, before and after translation",
         subtitle = sprintf(paste0("Grey = the donor protease's own value. Coloured = the value written into the matrix, after adding that ",
                                   "protein's offset. %d cells; median offset %+.2f log2. Ordered by the written value."),
                            nrow(d), med),
         x = "log2 intensity", y = "imputed cells (ordered by intensity)") +
    theme_adv() + theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())
}

# 12. Protein level: do the filled values land inside the protein's measured range, or below it?
# This is the MNAR question, and it is the one a reviewer will ask.
p_impute_protein_dumbbell <- function(res) {
  im <- res$imputation; if (is.null(im) || is.null(im$ledger) || !nrow(im$ledger)) return(NULL)
  fm <- res$final_matrix
  meas <- fm[value_source == "measured" & pg_canonical %in% unique(im$ledger$pg_canonical),
             .(lo = stats::quantile(value, .05, na.rm = TRUE),
               mid = stats::median(value, na.rm = TRUE),
               hi = stats::quantile(value, .95, na.rm = TRUE),
               n_meas = .N), by = pg_canonical]
  if (!nrow(meas)) return(NULL)
  impd <- im$ledger[, .(n_imp = .N, imp_med = stats::median(imputed_value, na.rm = TRUE),
                        imp_lo = min(imputed_value, na.rm = TRUE),
                        imp_hi = max(imputed_value, na.rm = TRUE)), by = pg_canonical]
  d <- merge(meas, impd, by = "pg_canonical")
  setorder(d, mid)
  d[, pg := factor(pg_canonical, levels = pg_canonical)]
  below <- d[imp_med < mid, .N]
  ggplot(d, aes(y = pg)) +
    geom_segment(aes(x = lo, xend = hi, yend = pg), colour = GREY, linewidth = 1.6, alpha = .55) +
    geom_point(aes(x = mid), colour = "grey25", size = 1.6) +
    geom_segment(aes(x = imp_lo, xend = imp_hi, yend = pg), colour = OK_COL, linewidth = .5, alpha = .8) +
    geom_point(aes(x = imp_med), colour = OK_COL, size = 1.8) +
    labs(title = "Stage 5 — where the filled values land inside each protein",
         subtitle = sprintf(paste0("Grey bar = the protein's measured 5–95%% range in the recommended channel, dark dot = its median. ",
                                   "Green = the imputed values. %d of %d proteins are filled below their own median, which is what ",
                                   "missing-not-at-random looks like."), below, nrow(d)),
         x = "log2 intensity", y = "imputed proteins (ordered by measured median)") +
    theme_adv() + theme(axis.text.y = element_text(size = 5.5), panel.grid.major.y = element_blank())
}
