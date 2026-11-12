# recommend_rollup.R — Part A decision: for each shared protein, trypsin or HTA.
#
# Extends the .criterion_winner head-to-head pattern of Concordance_Framework/engine/recommend.R with
# WEIGHTS (not all evidence is equally decisive) and with a confidence tier in place of an abstention.
# The user's constraint is explicit: "neither cannot be an option" — every shared protein gets a call.
# Where the evidence is thin that shows up as recommendation_confidence = "low", never as a refusal.
#
#   score_channel = (weight won by that channel) / (weight of all USABLE criteria)
# A criterion is usable only when the metric is finite in BOTH channels; unusable criteria leave the
# numerator AND the denominator, so a protein judged on 4 criteria is scored on the same 0-1 scale as
# one judged on 8 (and n_usable_criteria is reported so the reader can discount it). Exact ties award
# no weight to either side but stay in the denominator: a protein where the channels tie everywhere
# scores 0 vs 0 and lands in "low" confidence, which is exactly right.
#
# Two flags keep the output honest:
#   consequential   TRUE when the frozen concordance class is `discordant` — the channels disagree, so
#                   the choice changes the answer. FALSE when concordant: the call is near-cosmetic.
#                   Lets us report how many recommendations MATTER instead of implying all of them do.
#   flag_low_independence  carried from channel_independence.R: the two channels rest on largely the
#                   same peptides, so there is little to choose between.
#
# `drop_criteria` powers the held-out validation on the reproducibility data: rebuild the score with
# the CV and ICC criteria removed, then ask whether it still predicts the channel with the lower
# MEASURED technical CV. Dropping them is what makes that a test rather than a tautology.
#
# Depends on: data.table.

suppressPackageStartupMessages({ library(data.table) })

# criterion definition: metric column, direction, and an optional fallback metric used when the
# primary is unavailable in EITHER channel (both sides must use the same metric to be comparable).
CRITERIA <- list(
  completeness     = list(metric = "protein_detection_rate",   higher = TRUE),
  depth            = list(metric = "n_features",               higher = TRUE),
  coherence        = list(metric = "median_sibling_coherence", higher = TRUE),
  self_consistency = list(metric = "splithalf_reliability",    higher = TRUE,
                          fallback = "frac_prec_coherent"),
  rollup_stability = list(metric = "lopo_rank_stability",      higher = TRUE),
  precision        = list(metric = "median_rollup_se",         higher = FALSE),
  technical_cv     = list(metric = "frozen_technical_cv",      higher = FALSE),
  replicate_icc    = list(metric = "icc_repro",                higher = TRUE)
)

# Wide (one row per protein) view of the evidence table, channel columns suffixed _a / _b.
evidence_wide <- function(evidence, labelA = "trypsin", labelB = "hta") {
  a <- evidence[protease == labelA]; b <- evidence[protease == labelB]
  drop <- c("protease")
  a <- a[, setdiff(names(a), drop), with = FALSE]; b <- b[, setdiff(names(b), drop), with = FALSE]
  merge(a, b, by = "pg_canonical", suffixes = c("_a", "_b"), all = FALSE)
}

.pick_metric <- function(w, cr) {
  ca <- paste0(cr$metric, "_a"); cb <- paste0(cr$metric, "_b")
  usable <- ca %in% names(w) & cb %in% names(w)
  if (!usable) return(NULL)
  va <- w[[ca]]; vb <- w[[cb]]
  if (!is.null(cr$fallback)) {
    fa <- paste0(cr$fallback, "_a"); fb <- paste0(cr$fallback, "_b")
    if (fa %in% names(w) && fb %in% names(w)) {
      need <- !(is.finite(va) & is.finite(vb))
      va[need] <- w[[fa]][need]; vb[need] <- w[[fb]][need]
    }
  }
  list(va = va, vb = vb)
}

# Weighted head-to-head over shared proteins.
#   evidence      : precursor_evidence() output (long, pg x protease)
#   weights       : named list, one entry per criterion in CRITERIA
#   drop_criteria : criteria to exclude entirely (held-out validation)
#   class_tbl     : optional data.table(pg_canonical, concordance_class_final) — frozen, read-only
#   indep         : optional channel_independence() output
recommend_rollup <- function(evidence, weights, labelA = "trypsin", labelB = "hta",
                             drop_criteria = character(), class_tbl = NULL, indep = NULL,
                             high_margin = 0.34, high_min_criteria = 4L,
                             moderate_margin = 0.17, moderate_min_criteria = 3L,
                             selection_mode = c("completeness", "multi_criteria"),
                             completeness_high_gap = 0.20, completeness_moderate_gap = 0.05) {
  selection_mode <- match.arg(selection_mode)
  w <- evidence_wide(evidence, labelA, labelB)
  if (!nrow(w)) return(data.table())
  crits <- CRITERIA[setdiff(names(CRITERIA), drop_criteria)]
  n <- nrow(w)
  won_a <- numeric(n); won_b <- numeric(n); usable_w <- numeric(n); n_usable <- integer(n)
  best_a <- rep(NA_character_, n); best_b <- rep(NA_character_, n)
  best_aw <- numeric(n); best_bw <- numeric(n)
  for (nm in names(crits)) {
    cr <- crits[[nm]]; wt <- as.numeric(weights[[nm]] %||% 1)
    mv <- .pick_metric(w, cr); if (is.null(mv)) next
    va <- mv$va; vb <- mv$vb
    ok <- is.finite(va) & is.finite(vb)
    a_wins <- ok & (if (cr$higher) va > vb else va < vb)
    b_wins <- ok & (if (cr$higher) vb > va else vb < va)
    w[, (paste0("win_", nm)) := fifelse(!ok, "na", fifelse(a_wins, "A", fifelse(b_wins, "B", "tie")))]
    w[, (paste0("val_", nm, "_a")) := va][, (paste0("val_", nm, "_b")) := vb]
    won_a <- won_a + wt * a_wins; won_b <- won_b + wt * b_wins
    usable_w <- usable_w + wt * ok; n_usable <- n_usable + as.integer(ok)
    upd_a <- a_wins & wt > best_aw; best_a[upd_a] <- nm; best_aw[upd_a] <- wt
    upd_b <- b_wins & wt > best_bw; best_b[upd_b] <- nm; best_bw[upd_b] <- wt
  }
  w[, `:=`(score_a = fifelse(usable_w > 0, won_a / usable_w, NA_real_),
           score_b = fifelse(usable_w > 0, won_b / usable_w, NA_real_),
           usable_weight = usable_w, n_usable_criteria = n_usable)]
  w[!is.finite(score_a), score_a := 0][!is.finite(score_b), score_b := 0]
  w[, recommendation_margin := abs(score_a - score_b)]

  # deterministic tie-break: detection rate -> feature count -> trypsin
  w[, recommended_channel := fifelse(score_a > score_b, labelA,
                             fifelse(score_b > score_a, labelB, NA_character_))]
  w[, tie_break_used := is.na(recommended_channel)]
  w[is.na(recommended_channel), recommended_channel :=
      fifelse(protein_detection_rate_a > protein_detection_rate_b, labelA,
      fifelse(protein_detection_rate_b > protein_detection_rate_a, labelB,
      fifelse(n_features_a > n_features_b, labelA,
      fifelse(n_features_b > n_features_a, labelB, labelA))))]   # final fallback: trypsin

  # confidence tier — this is where an abstention WOULD have gone; it is reported, not acted on
  w[, winner_valid := fifelse(recommended_channel == labelA,
                              validity_flag_a == "valid", validity_flag_b == "valid")]
  w[is.na(winner_valid), winner_valid := FALSE]
  w[, recommendation_confidence := fifelse(
      recommendation_margin >= high_margin & n_usable_criteria >= high_min_criteria & winner_valid, "high",
      fifelse(recommendation_margin >= moderate_margin & n_usable_criteria >= moderate_min_criteria,
              "moderate", "low"))]

  # rationale: the highest-weight criterion each channel won, with its numbers
  fmt <- function(x) ifelse(is.finite(x), formatC(x, digits = 3, format = "g"), "NA")
  get_val <- function(crit, side) {
    out <- rep(NA_real_, nrow(w))
    for (nm in unique(crit[!is.na(crit)])) {
      cc <- paste0("val_", nm, "_", side)
      if (cc %in% names(w)) { idx <- which(crit == nm); out[idx] <- w[[cc]][idx] }
    }
    out
  }
  wa <- get_val(best_a, "a"); wb <- get_val(best_a, "b")
  la <- get_val(best_b, "a"); lb <- get_val(best_b, "b")
  w[, rationale := {
    winner <- recommended_channel
    top_w <- fifelse(winner == labelA, best_a, best_b)
    tw_a  <- fifelse(winner == labelA, wa, la); tw_b <- fifelse(winner == labelA, wb, lb)
    opp   <- fifelse(winner == labelA, best_b, best_a)
    op_a  <- fifelse(winner == labelA, la, wa); op_b <- fifelse(winner == labelA, lb, wb)
    paste0(winner, " wins on ", fifelse(is.na(top_w), "no criterion", top_w),
           " (", labelA, " ", fmt(tw_a), " vs ", labelB, " ", fmt(tw_b), ")",
           fifelse(is.na(opp), "; unopposed",
                   paste0("; ", fifelse(winner == labelA, labelB, labelA), " leads on ", opp,
                          " (", labelA, " ", fmt(op_a), " vs ", labelB, " ", fmt(op_b), ")")),
           "; score ", fmt(score_a), " vs ", fmt(score_b), " over ", n_usable_criteria, " criteria")
  }]

  # -------------------------------------------------------------------------------------------------
  # SELECTION MODE. Everything above computed the WEIGHTED VOTE, which the manuscript describes as the
  # general case the app offers. The manuscript's actual method — used throughout the paper and its
  # figures — is simpler: choose the channel with the higher COMPLETENESS (mean precursor detection
  # rate across samples), held constant per protein. We keep the vote's decision as `vote_*` for
  # reporting, and in `completeness` mode override the primary decision with the completeness winner.
  # In `multi_criteria` mode the vote IS the primary decision (backward-compatible).
  w[, `:=`(vote_channel = recommended_channel,
           vote_confidence = recommendation_confidence,
           vote_margin = recommendation_margin)]
  if (selection_mode == "completeness") {
    dra <- w$protein_detection_rate_a; drb <- w$protein_detection_rate_b
    # completeness winner, with the same deterministic tie-break cascade the vote used
    w[, recommended_channel :=
        fifelse(dra > drb, labelA,
        fifelse(drb > dra, labelB,
        fifelse(n_features_a > n_features_b, labelA,
        fifelse(n_features_b > n_features_a, labelB, labelA))))]   # final fallback: trypsin
    w[, tie_break_used := !(is.finite(dra) & is.finite(drb) & dra != drb)]
    w[, recommendation_margin := abs(dra - drb)]
    # confidence is tiered on the size of the completeness gap (an app-specific tier: the manuscript
    # defines high/moderate/low only for the weighted vote). The winning channel must also be `valid`.
    w[, winner_valid := fifelse(recommended_channel == labelA,
                                validity_flag_a == "valid", validity_flag_b == "valid")]
    w[is.na(winner_valid), winner_valid := FALSE]
    w[, recommendation_confidence := fifelse(
        is.finite(recommendation_margin) & recommendation_margin >= completeness_high_gap & winner_valid, "high",
        fifelse(is.finite(recommendation_margin) & recommendation_margin >= completeness_moderate_gap,
                "moderate", "low"))]
  }
  w[, selection_mode := selection_mode]
  w[, selection_disagrees_with_vote := recommended_channel != vote_channel]

  if (!is.null(class_tbl) && nrow(class_tbl)) {
    cls <- unique(class_tbl[, .(pg_canonical, concordance_class_final)])
    w <- merge(w, cls, by = "pg_canonical", all.x = TRUE)
  } else w[, concordance_class_final := NA_character_]
  w[, consequential := !is.na(concordance_class_final) & concordance_class_final == "discordant"]
  if (!is.null(indep) && nrow(indep))
    w <- merge(w, indep[, .(pg_canonical, shared_peptide_frac, complementarity, det_phi,
                            flag_low_independence)], by = "pg_canonical", all.x = TRUE)
  setorder(w, -recommendation_margin)
  w[]
}

`%||%` <- function(a, b) if (is.null(a)) b else a
