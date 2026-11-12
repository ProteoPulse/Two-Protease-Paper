# pipeline.R — the whole analysis, as one function. The Shiny app and the CLI both call this, so the
# GUI and the command line cannot drift apart.
#
# SHAPE OF THE FRAMEWORK
#   PR matrix ─┐
#              ├─► Stage 1 qualify precursors  ─► Stage 2 evidence per protein x protease
#   PG matrix ─┘                                        │
#        │                                              ▼
#        │                                   Stage 3 recommend: trypsin or HTA
#        ├──────► Stage 4 concordance ───────────────► (gates + `consequential` flag)
#        │                                              │
#        └──────────────────────────────────────────────┴─► the recommended PG value per protein
#                                                             │
#                                                    Stage 5 optional imputation
#
# WHAT WE DO NOT DO: we never emit our own protein quant. The deliverable is a decision about which
# existing PG-matrix value to trust. The median polish that appears in Stage 2 is a QC probe whose
# output is reported as a diagnostic and, by default, kept out of the score entirely.

suppressPackageStartupMessages({ library(data.table) })

.PROGRESS <- function(msg, frac) {
  f <- getOption("advisor.progress", NULL)
  if (is.function(f)) f(msg, frac) else cat(sprintf("[%3.0f%%] %s\n", 100 * frac, msg))
}

run_advisor <- function(pr_a, pr_b, pg_a, pg_b, cfg, prec = NULL, pg = NULL, pg_run = NULL) {
  LA <- cfg$label_a; LB <- cfg$label_b
  set.seed(as.integer(cfg$seed))
  ing <- cfg$ingest
  notes <- character()

  ## ---- Stage 0: ingest --------------------------------------------------------------------------
  if (is.null(prec) || is.null(pg)) {
    .PROGRESS("Reading precursor matrices", 0.03)
    prec <- rbind(
      load_pr_matrix(pr_a, LA, ing$run_include, ing$run_exclude, ing$run_sample_regex, as.integer(ing$min_runs_obs)),
      load_pr_matrix(pr_b, LB, ing$run_include, ing$run_exclude, ing$run_sample_regex, as.integer(ing$min_runs_obs)))
    .PROGRESS("Reading protein-group matrices", 0.16)
    pg_run <- rbind(
      load_pg_matrix(pg_a, LA, ing$run_include, ing$run_exclude, ing$run_sample_regex),
      load_pg_matrix(pg_b, LB, ing$run_include, ing$run_exclude, ing$run_sample_regex))
    pg <- pg_to_sample_grain(pg_run)
  }
  notes <- c(notes, paste0(
    "Input is a pivoted pr_matrix, which DIA-NN has already filtered at 1% q-value and from which ",
    "per-run q-values are absent. The qvalue rule is therefore not applicable; detection, CV and ",
    "coherence rules are unaffected."))

  samples_shared <- intersect(pg[protease == LA, unique(sample_id)], pg[protease == LB, unique(sample_id)])
  if (!length(samples_shared)) stop("the two matrices share no sample columns after run filtering")
  pg <- pg[sample_id %in% samples_shared]
  if (!is.null(pg_run)) pg_run <- pg_run[sample_id %in% samples_shared]
  prec <- prec[sample_id %in% samples_shared]
  have_reps <- prec[, uniqueN(run)] > prec[, uniqueN(sample_id)]

  ## proteins each channel sees at PROTEIN level, and the shared set we adjudicate
  seen <- pg[, .(n_ch = uniqueN(protease)), by = pg_canonical]
  shared_pgs <- seen[n_ch == 2L, pg_canonical]
  only_a <- setdiff(pg[protease == LA, unique(pg_canonical)], shared_pgs)
  only_b <- setdiff(pg[protease == LB, unique(pg_canonical)], shared_pgs)
  detection <- data.table(
    category = c("shared (adjudicated)", paste0(LA, "-only"), paste0(LB, "-only")),
    n_proteins = c(length(shared_pgs), length(only_a), length(only_b)),
    action = c("recommend a channel", "use it as-is; nothing to decide", "use it as-is; nothing to decide"))
  .PROGRESS(sprintf("%d shared proteins to adjudicate", length(shared_pgs)), 0.28)

  ## ---- Stage 1: qualification ------------------------------------------------------------------
  .PROGRESS("Stage 1 — qualifying precursors", 0.32)
  q <- cfg$qualify
  thr <- list(detection_rate_min = q$detection_rate_min, median_cv_max = q$median_cv_max,
              coherence_min = q$coherence_min, require_proteotypic = isTRUE(q$require_proteotypic),
              qvalue_max = NULL)
  mets <- precursor_metrics(prec, have_replicates = have_reps)
  qual <- qualify_apply(mets, thr)
  qual_ids <- qual[qualifies == TRUE, unique(precursor_id)]
  sweep <- tryCatch(threshold_sensitivity(mets, thr,
                      detection_grid = unlist(q$sweep$detection_grid),
                      coherence_grid = unlist(q$sweep$coherence_grid)),
                    error = function(e) NULL)

  ## ---- Stage 2: evidence (shared proteins only — that is where the ambiguity is) ----------------
  .PROGRESS("Stage 2 — collapsing peptide variants", 0.40)
  prec_shared <- prec[pg_canonical %in% shared_pgs]
  qual_shared <- prec_shared[precursor_id %in% qual_ids]
  ch <- collapse_peptide_hierarchy(qual_shared, decisions = NULL, corr_min = 0.70)

  roll <- data.table(pg_canonical = character(), protease = character(), run = character(),
                     sample_id = character(), theta_hat = numeric(), se_hat = numeric(),
                     n_obs = integer(), validity_flag = character(), n_features = integer(),
                     estimator = character())
  if (isTRUE(cfg$probe$enabled)) {
    .PROGRESS("Stage 2 — internal QC probe (median polish; diagnostic only)", 0.50)
    roll <- rollup_all(ch$collapsed, detection_min = 0.50, min_features_for_tmp = 2L)
  }

  .PROGRESS("Stage 2 — measuring evidence per protein and protease", 0.56)
  ev <- cfg$evidence
  fa <- if (nzchar(cfg$fasta %||% "") && file.exists(cfg$fasta)) fasta_index(cfg$fasta) else NULL
  evidence <- precursor_evidence(
    prec_shared, qual_ids, ch$collapsed, roll, cv_ref = NULL, icc_lookup = NULL, fa_index = fa,
    pgs = shared_pgs, labelA = LA, labelB = LB,
    splithalf_reps = as.integer(ev$splithalf_reps),
    splithalf_min_paired = as.integer(ev$splithalf_min_paired),
    lopo_max_features = as.integer(ev$lopo_max_features),
    jaccard_min_units = as.integer(ev$jaccard_min_units), seed = as.integer(cfg$seed))

  ## Replicate-based criteria, when the data has replicates. Both are measured on the PG value we are
  ## recommending — NOT on the internal probe — so they stay available even with the probe switched
  ## off, and so they describe the quantity the user will actually use.
  if (have_reps && !is.null(pg_run) && nrow(pg_run)) {
    cvm <- pg_technical_cv(pg_run[pg_canonical %in% shared_pgs])
    evidence <- merge(evidence, cvm[, .(pg_canonical, protease, frozen_technical_cv = technical_cv_measured)],
                      by = c("pg_canonical","protease"), all.x = TRUE)
    iccl <- tryCatch(build_icc_lookup(pg_run[pg_canonical %in% shared_pgs], value_col = "log2q",
                                      labelA = LA, labelB = LB), error = function(e) NULL)
    if (!is.null(iccl) && nrow(iccl)) {
      ra <- paste0("R_", LA); rb <- paste0("R_", LB)
      il <- rbind(iccl[, .(pg_canonical, protease = LA, icc_repro = get(ra))],
                  iccl[, .(pg_canonical, protease = LB, icc_repro = get(rb))])
      evidence <- merge(evidence, il, by = c("pg_canonical","protease"), all.x = TRUE)
    }
  }

  ## the probe's own diagnostic flag, kept separate from the score
  probe_flags <- NULL
  if (isTRUE(cfg$probe$enabled) && "lopo_rank_stability" %in% names(evidence)) {
    cut <- evidence[is.finite(lopo_rank_stability),
                    .(cut = stats::quantile(lopo_rank_stability, cfg$probe$fragile_quantile, na.rm = TRUE)),
                    by = protease]
    probe_flags <- merge(evidence[, .(pg_canonical, protease, n_features, lopo_rank_stability,
                                      lopo_max_delta, median_rollup_se)], cut, by = "protease")
    probe_flags[, fragile_rollup := is.finite(lopo_rank_stability) & n_features >= 3L &
                                    lopo_rank_stability < cut]
  }

  ## ---- Stage 2b: channel independence -----------------------------------------------------------
  .PROGRESS("Stage 2 — checking the channels are independent evidence", 0.66)
  A_long_pre <- pg[protease == LA, .(pg_canonical, sample_id, log2q)]
  B_long_pre <- pg[protease == LB, .(pg_canonical, sample_id, log2q)]
  indep <- tryCatch(channel_independence(prec_shared, qual_ids, A_long_pre, B_long_pre, shared_pgs,
                                         samples_shared, labelA = LA, labelB = LB),
                    error = function(e) NULL)

  ## ---- Stage 4: concordance from the PG values --------------------------------------------------
  .PROGRESS("Stage 4 — do the two channels agree? (permutation null)", 0.72)
  cc <- cfg$concordance
  A_long <- A_long_pre; B_long <- B_long_pre

  # min_paired has to be REACHABLE. A default tuned for a discovery set (15) is impossible on a
  # reproducibility set with 8 biological samples, and the failure mode is silent and severe: every
  # protein falls out as `not_assessable`, the permutation null comes back empty, the threshold
  # becomes NA, and the run still reports "0 discordant" — which reads as "the channels never
  # disagree" when it actually means "we never tested". Clamp it, and say so in the notes.
  min_paired <- as.integer(cc$min_paired)
  n_s <- length(samples_shared)
  if (min_paired > n_s) {
    adj <- max(3L, as.integer(floor(0.75 * n_s)))
    msg <- sprintf(paste0("concordance.min_paired was %d, but only %d samples are shared, so no ",
      "protein could reach it. Clamped to %d. An agreement statistic on %d paired samples is weak: ",
      "treat the concordant/discordant split in this run as indicative, not decisive."),
      min_paired, n_s, adj, adj)
    notes <- c(notes, msg)
    .PROGRESS(sprintf("min_paired %d unreachable with %d samples -> clamped to %d", min_paired, n_s, adj), 0.71)
    min_paired <- adj
  }
  cmet <- compare_metrics(A_long, B_long, min_paired = min_paired, method = cc$method)
  # compare_metrics returns the statistics; eligibility_flags decides which proteins are even
  # ASSESSABLE (enough paired samples, enough distinct values, non-zero variance). A protein that
  # fails these is not "discordant" — it is unmeasurable, and classify_concordance requires the flag.
  cmet <- eligibility_flags(cmet, min_paired = min_paired, min_distinct = 5L)
  nullr <- perm_null_rhos(A_long, B_long, n_perm = as.integer(cc$n_perm),
                          min_paired = min_paired, seed = as.integer(cfg$seed),
                          method = cc$method)
  # A degenerate null must stop the run, not be papered over with an NA threshold.
  if (!length(nullr) || !any(is.finite(nullr)))
    stop(sprintf(paste0("the permutation null is empty: no protein reached min_paired = %d across ",
      "%d shared samples. Lower concordance.min_paired, or check the run/sample filters."),
      min_paired, n_s))
  thr_val <- cc$threshold
  if (is.null(thr_val) || is.na(thr_val))
    thr_val <- unname(stats::quantile(nullr, as.numeric(cc$null_percentile), na.rm = TRUE))
  if (!is.finite(thr_val))
    stop("the agreement threshold is not finite; the permutation null could not be summarised.")
  cfg$concordance$min_paired_used <- min_paired
  cls <- classify_concordance(cmet, threshold = thr_val)
  setnames(cls, "concordance_class", "concordance_class_final", skip_absent = TRUE)
  not_assessable <- cmet[concordance_evaluable == FALSE,
                         .(pg_canonical, concordance_class_final = "not_assessable",
                           not_assessable_reasons, n_paired, stat_value)]
  cls <- rbindlist(list(cls, not_assessable), fill = TRUE)

  ## ---- Stage 3: the recommendation ---------------------------------------------------------------
  .PROGRESS("Stage 3 — recommending a channel per protein", 0.80)
  drop <- if (isTRUE(cfg$probe$in_score)) character() else c("rollup_stability", "precision")
  conf <- cfg$recommend$confidence
  sel_mode <- cfg$recommend$selection_mode %||% "completeness"
  rec <- recommend_rollup(evidence, cfg$recommend$weights, LA, LB, drop_criteria = drop,
                          class_tbl = cls, indep = indep,
                          high_margin = conf$high_margin, high_min_criteria = conf$high_min_criteria,
                          moderate_margin = conf$moderate_margin,
                          moderate_min_criteria = conf$moderate_min_criteria,
                          selection_mode = sel_mode,
                          completeness_high_gap = as.numeric(cfg$recommend$completeness_high_gap %||% 0.20),
                          completeness_moderate_gap = as.numeric(cfg$recommend$completeness_moderate_gap %||% 0.05))
  notes <- c(notes, if (sel_mode == "completeness")
    sprintf(paste0("Channel selection = COMPLETENESS (the manuscript method): each protein takes the ",
      "channel with the higher mean precursor detection rate, held across all samples. The weighted ",
      "vote is still reported as vote_channel; it would have changed %d of %d calls."),
      sum(rec$selection_disagrees_with_vote, na.rm = TRUE), nrow(rec))
    else "Channel selection = MULTI_CRITERIA (the weighted-vote generalization the app offers).")
  if (!is.null(probe_flags)) {
    pf <- dcast(probe_flags, pg_canonical ~ protease, value.var = "fragile_rollup")
    setnames(pf, setdiff(names(pf), "pg_canonical"), paste0("fragile_rollup_", setdiff(names(pf), "pg_canonical")))
    rec <- merge(rec, pf, by = "pg_canonical", all.x = TRUE)
  }

  ## ---- the deliverable: which PG value to use ----------------------------------------------------
  .PROGRESS("Assembling the recommended abundance matrix", 0.86)
  chosen <- merge(pg, rec[, .(pg_canonical, recommended_channel, recommendation_confidence,
                              recommendation_margin, consequential)],
                  by = "pg_canonical")[protease == recommended_channel]
  chosen[, `:=`(value_source = "measured", source_channel = protease)]
  advised <- chosen[, .(pg_canonical, gene, sample_id, value = log2q,
                        recommended_channel, source_channel,
                        donor_channel = NA_character_, donor_value = NA_real_,
                        offset_used = NA_real_,
                        recommendation_confidence, recommendation_margin, consequential,
                        value_source)]
  ## proteins only one channel saw: pass through untouched, clearly labelled
  passthru <- pg[pg_canonical %in% c(only_a, only_b),
                 .(pg_canonical, gene, sample_id, value = log2q,
                   recommended_channel = protease, source_channel = protease,
                   donor_channel = NA_character_, donor_value = NA_real_, offset_used = NA_real_,
                   recommendation_confidence = "not applicable (single channel)",
                   recommendation_margin = NA_real_, consequential = FALSE, value_source = "measured")]
  final_matrix <- rbindlist(list(advised, passthru), use.names = TRUE)

  ## ---- Stage 5: optional imputation --------------------------------------------------------------
  imp <- NULL
  if (isTRUE(cfg$impute$enabled)) {
    .PROGRESS("Stage 5 — cross-protease gap filling", 0.92)
    imp <- impute_from_pg(pg, rec, cls, cfg$impute, LA, LB, samples_shared)
    if (!is.null(imp$matrix) && nrow(imp$matrix)) {
      # An imputed cell belongs to a protein that already has a decision, so it inherits that
      # protein's recommendation verbatim. Leaving these fields blank (as an earlier version did)
      # made the imputed rows look like a separate, unadjudicated class of value.
      im <- merge(copy(imp$matrix),
                  rec[, .(pg_canonical, recommended_channel, recommendation_confidence,
                          recommendation_margin, consequential)],
                  by = "pg_canonical", all.x = TRUE)
      gmap <- unique(pg[, .(pg_canonical, gene)])[!duplicated(pg_canonical)]
      im <- merge(im, gmap, by = "pg_canonical", all.x = TRUE)
      final_matrix <- rbindlist(list(
        final_matrix[!paste(pg_canonical, sample_id) %in%
                       paste(im$pg_canonical, im$sample_id)],
        im), use.names = TRUE, fill = TRUE)
    }
  }

  ## ---- Stage 6: reproducibility / well-quantified counts ----------------------------------------
  # Brings the manuscript's well-quantified framing into the app and PICKS THE METRIC for the data:
  # inter-day CV (corrected) when replicates + a day map exist, else cohort completeness. Combined is
  # always the UNION; the proteases are never averaged. Never fatal — a reporting stage.
  repro <- tryCatch({
    if (!is.null(pg_run) && nrow(pg_run)) {
      rs <- reproducibility_stage(pg_run[pg_canonical %in% c(shared_pgs, only_a, only_b)],
              imp_mat = if (!is.null(imp)) imp$matrix else NULL,
              LA, LB, if (is.null(cfg$reproducibility)) list() else cfg$reproducibility, have_reps)
      .PROGRESS(sprintf("Stage 6 — well quantified (%s)", rs$mode), 0.97)
      rs
    } else NULL
  }, error = function(e) { notes <<- c(notes, paste0("reproducibility stage skipped: ", conditionMessage(e))); NULL })

  ## ---- invariant: ONE channel per protein, across every sample ----------------------------------
  # The deliverable is a per-protein decision, not a per-cell one. If a protein ever carried two
  # source channels the matrix would silently become a per-sample switch, which is a different (and
  # much weaker) claim than the one the report makes. Check it rather than trust it.
  chk <- final_matrix[, .(n_ch = uniqueN(source_channel[!is.na(source_channel)])), by = pg_canonical]
  if (any(chk$n_ch > 1L))
    stop(sprintf(paste0("internal check failed: %d protein(s) carry more than one source channel in ",
      "the final matrix. The recommendation must be constant across samples for a protein."),
      sum(chk$n_ch > 1L)))
  setcolorder(final_matrix, intersect(c("pg_canonical","gene","sample_id","value",
    "recommended_channel","source_channel","value_source","donor_channel","donor_value",
    "offset_used","recommendation_confidence","recommendation_margin","consequential"),
    names(final_matrix)))

  .PROGRESS("Done", 1)
  list(prec = prec, pg = pg, pg_run = pg_run, samples = samples_shared, have_replicates = have_reps,
       detection = detection, shared_pgs = shared_pgs,
       qual = qual, qual_sweep = sweep, metrics = mets,
       evidence = evidence, probe_flags = probe_flags, indep = indep, roll = roll,
       concordance = cls, null_rhos = nullr, threshold = thr_val,
       rec = rec, final_matrix = final_matrix, imputation = imp, reproducibility = repro,
       ingest_summary = ingest_summary(prec, pg, LA, LB), notes = notes, cfg = cfg)
}
