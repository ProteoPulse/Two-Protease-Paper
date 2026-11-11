# impute_pg.R — Stage 5: fill gaps in the recommended channel from the other protease's PG value.
#
# The earlier implementation filled gaps in OUR OWN roll-up. This one fills gaps in the vendor PG
# matrix, because that is what the framework now delivers. The statistics are otherwise the Aug-20
# design with the Aug-24 benchmark's conclusions applied:
#
#   * default offset model is SAMPLE-kNN, not the per-protein median. The offset is computed from the
#     k co-observed samples nearest the target sample in donor-profile space, which catches offsets
#     that drift with sample composition. It beat the per-protein median in both conditions on paired
#     per-cell tests (p = 3.5e-07 Mix30, 1.0e-08 Mix20) — by a small margin (~0.04 log2).
#   * BORROWED offsets (the old tier 3, taking an offset from similar proteins when a protein has too
#     few co-observed samples) are OFF by default. Benchmarked at 1.10-1.12 log2 overall and 4.53 on
#     gap-bearing proteins, they are worse than leaving the value missing.
#   * no slope is offered. A per-protein slope is the best model under MCAR and among the worst under
#     MNAR: it is fitted on the observed range and asked to extrapolate into the censored tail.
#
# THE HARD BLOCK STAYS. If the two channels disagree across samples (concordance class `discordant`)
# the donor is not a valid surrogate, and translating it would manufacture a trend the recipient never
# showed. Those gaps stay missing. This is a decision, not a knob.
#
# Every filled cell lands in a ledger with its donor value, offset, model, n, MAD, tier and the gates
# it passed. Every REFUSED gap lands in the same ledger with the reason. On real data the refusals are
# the more informative half.

suppressPackageStartupMessages({ library(data.table) })

# Lin's concordance correlation coefficient on paired vectors.
.ccc <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); if (sum(ok) < 3L) return(NA_real_)
  x <- x[ok]; y <- y[ok]
  vx <- stats::var(x); vy <- stats::var(y)
  if (!is.finite(vx) || !is.finite(vy) || (vx + vy) == 0) return(NA_real_)
  2 * stats::cov(x, y) / (vx + vy + (mean(x) - mean(y))^2)
}

# Sample-sample similarity in DONOR space, over proteins the donor measures nearly everywhere.
# A sparse basis would make the distance a missingness artefact, so the basis is deliberately dense.
.sample_similarity <- function(pg_channel, samples, min_detection = 0.90) {
  d <- pg_channel[is.finite(log2q)]
  ns <- length(samples)
  det <- d[, .(det = uniqueN(sample_id) / ns), by = pg_canonical]
  keep <- det[det >= min_detection, pg_canonical]
  if (length(keep) < 10L) keep <- det[order(-det)][seq_len(min(50L, .N)), pg_canonical]
  m <- dcast(d[pg_canonical %in% keep], pg_canonical ~ sample_id, value.var = "log2q")
  M <- as.matrix(m[, -1, drop = FALSE])
  suppressWarnings(stats::cor(M, use = "pairwise.complete.obs", method = "spearman"))
}

impute_from_pg <- function(pg, rec, cls, cfg_imp, labelA, labelB, samples) {
  k <- as.integer(cfg_imp$knn_k); min_paired <- as.integer(cfg_imp$min_paired_offset)
  G <- cfg_imp$gates; E <- cfg_imp$entry
  block <- unlist(cfg_imp$block_classes)
  model <- cfg_imp$offset_model %||% "sample_knn"

  # full (protein x sample) grid for the shared proteins, both channels side by side
  grid <- CJ(pg_canonical = rec$pg_canonical, sample_id = samples, unique = TRUE)
  w <- merge(grid, dcast(pg[pg_canonical %in% rec$pg_canonical],
                         pg_canonical + sample_id ~ protease, value.var = "log2q"),
             by = c("pg_canonical","sample_id"), all.x = TRUE)
  for (cc in c(labelA, labelB)) if (!cc %in% names(w)) w[, (cc) := NA_real_]

  w <- merge(w, rec[, .(pg_canonical, recipient = recommended_channel)], by = "pg_canonical")
  w[, donor := fifelse(recipient == labelA, labelB, labelA)]
  w[, theta_rec := fifelse(recipient == labelA, get(labelA), get(labelB))]
  w[, theta_don := fifelse(donor     == labelA, get(labelA), get(labelB))]

  ## ---- per-protein state: completeness, offset, agreement -------------------------------------
  ns <- length(samples)
  pp <- w[, .(completeness_recipient = sum(is.finite(theta_rec)) / ns,
              completeness_donor     = sum(is.finite(theta_don)) / ns,
              n_gaps                 = sum(!is.finite(theta_rec)),
              n_gaps_with_donor      = sum(!is.finite(theta_rec) & is.finite(theta_don)),
              n_paired               = sum(is.finite(theta_rec) & is.finite(theta_don)),
              offset                 = stats::median(theta_rec - theta_don, na.rm = TRUE),
              offset_mad             = stats::mad(theta_rec - theta_don, na.rm = TRUE),
              ccc                    = .ccc(theta_rec, theta_don)),
          by = .(pg_canonical, recipient, donor)]
  ct <- if (!is.null(cls) && nrow(cls))
    unique(cls[, .(pg_canonical, concordance_class_final,
                   kendall_tau = if ("stat_value" %in% names(cls)) stat_value else NA_real_)])
    else data.table(pg_canonical = pp$pg_canonical, concordance_class_final = NA_character_,
                    kendall_tau = NA_real_)
  pp <- merge(pp, ct, by = "pg_canonical", all.x = TRUE)

  ## ---- gates, each recorded separately so a refusal can always be explained --------------------
  pp[, gate_has_gap        := n_gaps > 0]
  pp[, gate_donor_covers   := completeness_donor >= as.numeric(E$min_completeness_donor)]
  pp[, gate_not_complete   := completeness_recipient <= as.numeric(E$max_completeness_recipient)]
  pp[, gate_class          := is.na(concordance_class_final) | !(concordance_class_final %in% block)]
  pp[, gate_paired         := n_paired >= min_paired]
  pp[, gate_agree          := (is.finite(ccc) & ccc >= G$ccc_min) |
                              (is.finite(kendall_tau) & kendall_tau >= G$tau_min)]
  pp[, gate_offset_stable  := is.finite(offset_mad) & offset_mad <= G$offset_mad_max]
  pp[, eligible := gate_has_gap & gate_donor_covers & gate_not_complete & gate_class &
                   gate_paired & gate_agree & gate_offset_stable]
  pp[, refusal_reason := fifelse(!gate_has_gap, "no gaps to fill",
                          fifelse(!gate_class, "channels disagree (class discordant)",
                          fifelse(!gate_donor_covers, "donor does not cover the protein",
                          fifelse(!gate_paired, "too few co-observed samples",
                          fifelse(!gate_agree, "channels do not track each other (ccc/tau gate)",
                          fifelse(!gate_offset_stable, "offset is not stable across samples",
                                  NA_character_))))))]

  ## ---- offsets ---------------------------------------------------------------------------------
  elig <- pp[eligible == TRUE, pg_canonical]
  gaps <- w[!is.finite(theta_rec) & is.finite(theta_don) & pg_canonical %in% elig]
  ledger <- data.table()
  if (nrow(gaps)) {
    if (identical(model, "sample_knn")) {
      # One similarity matrix per donor channel: a gap whose donor is trypsin finds its neighbours in
      # trypsin space, a gap whose donor is HTA in HTA space.
      S_by_donor <- lapply(setNames(nm = unique(gaps$donor)), function(dn)
        tryCatch(.sample_similarity(pg[protease == dn], samples), error = function(e) NULL))
      obs <- w[is.finite(theta_rec) & is.finite(theta_don),
               .(pg_canonical, sample_id = as.character(sample_id), off_s = theta_rec - theta_don)]
      by_pg <- split(obs, obs$pg_canonical)
      gaps[, offset_used := vapply(seq_len(.N), function(i) {
        o <- by_pg[[pg_canonical[i]]]
        if (is.null(o) || !nrow(o)) return(NA_real_)
        S <- S_by_donor[[donor[i]]]
        s <- as.character(sample_id[i])
        if (is.null(S) || !(s %in% rownames(S))) return(stats::median(o$off_s))
        sim <- S[s, o$sample_id]; ok <- is.finite(sim)
        if (!any(ok)) return(stats::median(o$off_s))
        idx <- order(-sim[ok])[seq_len(min(k, sum(ok)))]
        stats::median(o$off_s[which(ok)][idx])
      }, numeric(1))]
      gaps[, offset_model_used := "sample_knn"]
    } else {
      gaps <- merge(gaps, pp[, .(pg_canonical, offset)], by = "pg_canonical", all.x = TRUE)
      gaps[, `:=`(offset_used = offset, offset_model_used = "protein_median")]
    }
    gaps[, imputed_value := theta_don + offset_used]
    gaps <- gaps[is.finite(imputed_value)]
    ledger <- merge(gaps[, .(pg_canonical, sample_id, recipient, donor, theta_don,
                             offset_used, offset_model_used, imputed_value)],
                    pp[, .(pg_canonical, n_paired, offset_mad, ccc, kendall_tau,
                           concordance_class_final, completeness_recipient, completeness_donor)],
                    by = "pg_canonical")
    ledger[, tier := fifelse(!is.na(concordance_class_final) & concordance_class_final == "concordant",
                             1L, 2L)]
  }

  ## ---- refused gaps, recorded with the same care as the filled ones ---------------------------
  refused <- w[!is.finite(theta_rec) & !(pg_canonical %in% elig),
               .(pg_canonical, sample_id, recipient, donor, theta_don)]
  refused <- merge(refused, pp[, .(pg_canonical, refusal_reason, ccc, kendall_tau, n_paired,
                                   completeness_recipient, completeness_donor,
                                   concordance_class_final)],
                   by = "pg_canonical", all.x = TRUE)
  refused[is.na(refusal_reason), refusal_reason := "donor also missing in this sample"]
  refused[!is.finite(theta_don), refusal_reason := "donor also missing in this sample"]

  # LABELLING — this is the one place the decision grain could be misread, so it is explicit.
  # The framework picks ONE channel per protein and holds it across every sample. An imputed cell does
  # NOT switch channel: the donor observation is translated onto the RECIPIENT's scale by that
  # protein's offset, so the number still lives in the recommended channel's frame. Labelling the cell
  # with the donor made a protein look as though its channel changed sample-by-sample, which is
  # precisely the thing this framework never does. `source_channel` therefore stays the recommended
  # channel for every cell of the protein; the donor that supplied the evidence is recorded separately
  # in `donor_channel`, alongside the untranslated `donor_value` and the `offset_used`.
  mat <- if (nrow(ledger)) ledger[, .(pg_canonical, sample_id, value = imputed_value,
                                      source_channel = recipient, donor_channel = donor,
                                      donor_value = theta_don, offset_used,
                                      value_source = "imputed_cross_protease")]
         else data.table()
  summary <- data.table(
    n_shared = nrow(pp),
    n_cells_recipient = nrow(w),
    n_measured = w[is.finite(theta_rec), .N],
    pct_complete = round(100 * w[is.finite(theta_rec), .N] / nrow(w), 2),
    n_gaps = w[!is.finite(theta_rec), .N],
    n_gaps_with_donor = w[!is.finite(theta_rec) & is.finite(theta_don), .N],
    n_imputed = nrow(ledger),
    n_proteins_imputed = if (nrow(ledger)) ledger[, uniqueN(pg_canonical)] else 0L,
    n_refused = nrow(refused))
  list(per_protein = pp, ledger = ledger, refused = refused, matrix = mat, summary = summary,
       refusal_breakdown = refused[, .N, by = refusal_reason][order(-N)])
}
