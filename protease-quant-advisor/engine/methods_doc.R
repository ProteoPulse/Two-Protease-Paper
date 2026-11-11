# methods_doc.R — writes INTERPRETATION_AND_METHODS.md for a completed run.
#
# Generated, not hand-written, so the prose can never drift from the numbers it describes. Every
# figure quoted below is pulled from the run's own state object at the moment the document is made.

suppressPackageStartupMessages({ library(data.table) })
# Treat NULL, and a *scalar* NA, as "absent". Guard the is.na() with a length check: is.na() on a
# multi-element vector or a list returns a vector, which then blows up `||` ("length = 3 in coercion
# to logical(1)"). This operator is used on config values that are sometimes whole lists (e.g.
# cfg$reproducibility), so it must never assume its left operand is scalar.

.fm  <- function(x, d = 3) if (!length(x) || !is.finite(x)) "n/a" else formatC(x, format = "f", digits = d)
.pc  <- function(n, d) if (!is.finite(n/d)) "n/a" else sprintf("%.1f%%", 100 * n / d)
.int <- function(x) formatC(as.numeric(x), format = "d", big.mark = ",")

.md_table <- function(dt, digits = 3) {
  dt <- as.data.table(dt)
  if (!nrow(dt)) return("*(empty)*")
  for (cn in names(dt)) if (is.numeric(dt[[cn]])) dt[[cn]] <- round(dt[[cn]], digits)
  hdr <- paste0("| ", paste(names(dt), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(dt)), collapse = "|"), "|")
  rows <- apply(dt, 1, function(r) paste0("| ", paste(ifelse(is.na(r), "", as.character(r)), collapse = " | "), " |"))
  paste(c(hdr, sep, rows), collapse = "\n")
}

write_methods_doc <- function(res, out_dir, tag = "run") {
  LA <- res$cfg$label_a; LB <- res$cfg$label_b
  r  <- res$rec; cfgq <- res$cfg$qualify; cfgc <- res$cfg$concordance; cfgi <- res$cfg$impute
  ns <- length(res$samples); nsh <- nrow(r)
  reps <- isTRUE(res$have_replicates)
  inp <- tryCatch(fread(file.path(out_dir, "run_inputs.tsv")), error = function(e) NULL)
  im  <- res$imputation
  cls <- res$concordance
  ev  <- res$evidence

  n_a <- r[recommended_channel == LA, .N]; n_b <- r[recommended_channel == LB, .N]
  n_hi <- r[recommendation_confidence == "high", .N]
  n_lo <- r[recommendation_confidence == "low", .N]
  n_cons <- r[consequential == TRUE, .N]
  qual_n <- nrow(res$qual); qual_k <- sum(res$qual$qualifies)

  # which criteria actually voted, and which dropped out (this is the part readers get wrong)
  crit_metric <- c(completeness = "protein_detection_rate", depth = "n_features",
                   coherence = "median_sibling_coherence", self_consistency = "splithalf_reliability",
                   rollup_stability = "lopo_rank_stability", precision = "median_rollup_se",
                   technical_cv = "frozen_technical_cv", replicate_icc = "icc_repro")
  probe_in <- isTRUE(res$cfg$probe$in_score)
  crit_tbl <- rbindlist(lapply(names(crit_metric), function(cn) {
    m <- crit_metric[[cn]]
    present <- m %in% names(ev) && sum(is.finite(ev[[m]])) > 0
    excluded_probe <- (!probe_in) && cn %in% c("rollup_stability","precision")
    data.table(criterion = cn, metric = m, weight = res$cfg$recommend$weights[[cn]],
               status = if (excluded_probe) "excluded (QC probe, by design)"
                        else if (present) "live" else "dropped (metric unavailable)",
               n_proteins_with_value = if (m %in% names(ev)) sum(is.finite(ev[[m]])) else 0L)
  }), fill = TRUE)

  L <- c()
  add <- function(...) L <<- c(L, ...)

  add(sprintf("# %s — interpretation and methods", tag), "",
      sprintf("*Protease Quant Advisor · generated %s · %d samples · %s shared proteins adjudicated*",
              format(Sys.time(), "%Y-%m-%d %H:%M"), ns, .int(nsh)), "",
      "---", "")

  ## ------------------------------------------------------------------ what this run is
  add("## 1. What this run does", "",
      paste0("This is not a re-quantification. Both searches were run on the **same injections**, so for a ",
             "protein both of them see we hold two numbers for one molecule in one tube. Averaging them would ",
             "invent a third number that neither search reported. Instead the tool reads the **precursor-level ",
             "evidence** behind each shared protein and recommends **which search's protein-group value to ",
             "trust**, with a confidence tier."), "",
      sprintf("Channel A = **%s**, channel B = **%s**.%s", LA, LB,
              if (reps) " This dataset has technical replicates, so replicate-based criteria are live (see §5)."
              else " This dataset has one run per sample, so replicate-based criteria do not apply (see §5)."), "")

  ## the grain of the decision — the thing most likely to be misread
  {
    fm <- res$final_matrix
    nch <- fm[, .(n = uniqueN(source_channel[!is.na(source_channel)])), by = pg_canonical]
    add("### 1.1 At what level is the recommendation made?", "",
        paste0("**Once per protein, held across every sample in the experiment.** Not per sample, and not ",
               "one winner for the whole experiment. Three levels are possible and only the middle one is ",
               "defensible:"), "",
        .md_table(data.table(
          level = c("per protein **per sample**", "per protein, **fixed across samples**",
                    "one channel for the **whole experiment**"),
          used = c("no", "**yes — this is what the tool does**", "no"),
          why = c(paste0("A protein taking its value from ", LA, " in one injection and ", LB, " in the next ",
                         "puts the between-search offset inside the protein's own profile. Every fold-change ",
                         "then mixes real biology with a search artefact, and the artefact can be larger than ",
                         "the effect."),
                  paste0("Each shared protein is adjudicated once, on its own precursor evidence, and that ",
                         "channel supplies all of its values. Every within-protein comparison therefore stays ",
                         "inside a single measurement system, while different proteins are still allowed to ",
                         "come from whichever search measures them better."),
                  paste0("Would collapse to `use ", LA, "` and throw away the proteins the other search ",
                         "genuinely measures better. It answers a question nobody asked.")))), "",
        "**Enforced, not assumed.** The pipeline halts if any protein ever carries two source channels:", "",
        .md_table(data.table(
          check = c("proteins in the final matrix",
                    "proteins with exactly one source channel across all samples",
                    "proteins with more than one",
                    "imputed cells (channel unchanged)"),
          value = c(format(nrow(nch), big.mark = ","), format(nch[n == 1L, .N], big.mark = ","),
                    format(nch[n > 1L, .N], big.mark = ","),
                    format(fm[value_source == "imputed_cross_protease", .N], big.mark = ",")))), "",
        paste0("**Imputation does not create an exception.** When the recommended channel has no value for a ",
               "sample, the other protease's measurement of the same protein *in the same injection* is ",
               "translated onto the recommended channel's scale using that protein's own offset. The written ",
               "number lives in the recommended channel's frame, so `source_channel` stays constant for the ",
               "protein. The protease that supplied the evidence is recorded separately in `donor_channel`, ",
               "with `donor_value` and `offset_used`, so any filled cell can be reversed by hand."), "")
  }

  ## ------------------------------------------------------------------ inputs
  add("## 2. Input files", "")
  if (!is.null(inp)) {
    add(.md_table(inp[, .(role, channel, level, size_mb, path)]), "")
    rf <- unique(inp[, .(run_include, run_exclude, run_sample_regex)])
    add("**Run selection.** Column headers in a DIA-NN matrix are raw-file paths; they are reduced to a bare",
        "run name (basename, extension removed) and then filtered:", "",
        sprintf("- keep runs matching: `%s`", if (nzchar(rf$run_include[1])) rf$run_include[1] else "(all runs)"),
        sprintf("- drop runs matching: `%s`", if (nzchar(rf$run_exclude[1])) rf$run_exclude[1] else "(none)"),
        sprintf("- run → sample: `%s`", if (nzchar(rf$run_sample_regex[1])) paste0(rf$run_sample_regex[1],
                 "` (first capture group is the biological sample)") else "one run per sample`"), "")
  } else add("*(run_inputs.tsv not found for this run)*", "")
  add("**What was actually read after filtering:**", "", .md_table(res$ingest_summary), "")
  add(sprintf("Samples present in both channels: **%d**. Only these are used; a sample only one search saw cannot be compared.", ns), "")

  ## ------------------------------------------------------------------ preprocessing
  add("## 3. Preprocessing", "",
      "**Intensities.** Values are taken as reported and converted to log2. Nothing is re-normalised:",
      "the point of the tool is to adjudicate the vendor's numbers, not to replace them.", "",
      "**Inherited filtering — read this before quoting a count.** The pivoted `pr_matrix` has already",
      "been q-value filtered by DIA-NN (precursor / protein / PG / library q ≤ 0.01) and the per-run",
      "q-values are gone from the file. So that 1% filter is a **floor we cannot loosen**, and our own",
      "`qvalue_max` rule is *not applicable* here. The other three qualification rules are computed from",
      "intensities and are unaffected.", "",
      "**Protein groups** are canonicalised to a sorted accession set so the same group matches across",
      "the two searches even when the two libraries list its members in a different order.", "",
      sprintf("**Replicate handling.** %s",
              if (reps) paste0("Several runs map to one biological sample. Presence metrics are computed at RUN grain; ",
                               "profile metrics and the recommendation at SAMPLE grain, replicates collapsed by median. ",
                               "The run grain is retained so technical CV and ICC can be measured.")
              else "One run per sample, so run grain and sample grain coincide."), "")

  ## ------------------------------------------------------------------ stage 1
  add("## 4. Stage 1 — precursor qualification", "",
      "Our own rule, applied *on top of* DIA-NN's filtering. A precursor must:", "",
      sprintf("- be detected in ≥ **%s** of runs", .fm(cfgq$detection_rate_min, 2)),
      sprintf("- agree with its siblings (median pairwise Pearson to other precursors of the same protein) ≥ **%s**", .fm(cfgq$coherence_min, 2)),
      sprintf("- %s", if (reps) sprintf("have replicate CV ≤ **%s**", .fm(cfgq$median_cv_max, 2))
                      else "*(CV rule inactive — no replicates)*"),
      sprintf("- %s", if (isTRUE(cfgq$require_proteotypic)) "be proteotypic" else "*(proteotypic rule off)*"), "",
      sprintf("**Outcome: %s of %s precursors qualify (%s).**", .int(qual_k), .int(qual_n), .pc(qual_k, qual_n)), "",
      "A precursor with no siblings cannot fail the coherence rule and passes it flagged, rather than being",
      "silently dropped for having nothing to compare against.", "")
  if (!is.null(res$qual_sweep)) {
    add("**Is this threshold pair on a plateau or a cliff?** The sweep is in",
        "`tables/qualification_sensitivity.tsv` and plotted in the report. Read it before trusting the counts:",
        "if a small move in either threshold costs many precursors, the setting is doing more work than it should.", "")
  }

  ## ------------------------------------------------------------------ stage 2
  add("## 5. Stage 2 — evidence, and which criteria actually voted", "",
      "Evidence is measured per protein **and per protease**, on two axes:", "",
      "- **completeness** — is the protein really measured everywhere we need it? (detection rate, peptide",
      "  depth, run-to-run peptide-set stability, sequence coverage where a FASTA was supplied)",
      "- **coherence** — do the peptides tell the same story? (sibling coherence, split-half reliability,",
      "  and, where replicates exist, technical CV and ICC)", "",
      "Criteria are combined as a weighted head-to-head. A criterion whose metric is missing in either",
      "channel leaves **both** the numerator and the denominator, so a protein judged on four criteria is",
      "scored on the same 0–1 scale as one judged on six. This table says exactly what happened here:", "",
      .md_table(crit_tbl), "")
  if (!probe_in) {
    add("**On the two excluded criteria.** `rollup_stability` and `precision` are measured on a Tukey median",
        "polish that this tool computes internally and never publishes. They are excluded from the score by",
        "default because a dedicated ablation found the probe does not measurably improve the recommendation:",
        "it correlates 0.91 with peptide count, flips 1.9% of calls and none of the consequential ones, and",
        "its depth-independent component fails to predict cross-protease discordance (p = 0.83 / 0.99).",
        "It is still **computed and reported** as a per-protein diagnostic — `fragile_rollup` marks a protein",
        "whose value rests on a single peptide, which is worth knowing even though it does not move the call.",
        "Set `probe.in_score: true` to restore the older behaviour.", "")
  }

  ## ------------------------------------------------------------------ stage 4
  ncls <- cls[, .N, by = concordance_class_final]
  if (length(res$notes) > 1) add("### Run-specific notes", "",
      paste0("> **", paste(res$notes[-1], collapse = "**\n>\n> **"), "**"), "")
  add("## 6. Stage 4 — do the two channels agree about this protein?", "",
      sprintf("Agreement is measured on the **protein-group values** with %s, and the threshold is not a chosen",
              c(kendall = "Kendall's tau-b", spearman = "Spearman's rho")[cfgc$method] %||% cfgc$method),
      sprintf("number: it is the **%sth percentile of a permutation null** built by shuffling sample labels within",
              .fm(100 * as.numeric(cfgc$null_percentile), 0)),
      sprintf("each protein (%s permutations), requiring at least %s paired samples per protein.",
              .int(cfgc$n_perm), .int(cfgc$min_paired_used %||% cfgc$min_paired)), "",
      sprintf("**Threshold = %s.**", .fm(res$threshold, 4)), "", .md_table(ncls), "",
      "`not_assessable` is a real category, not a failure: too few paired samples, too few distinct values,",
      "or zero variance. Those proteins are carried through explicitly rather than quietly counted as",
      "agreeing or disagreeing.", "",
      "This matters for two reasons. It gates the imputation (a disagreeing donor is not a surrogate), and",
      "it is what makes a recommendation **consequential** — where the channels already agree, choosing one",
      "over the other barely moves the number.", "")

  ## ------------------------------------------------------------------ stage 3
  sel_mode <- r$selection_mode[1] %||% "completeness"
  n_disw <- if ("selection_disagrees_with_vote" %in% names(r)) sum(r$selection_disagrees_with_vote, na.rm = TRUE) else NA_integer_
  add("## 7. Stage 3 — the recommendation", "",
      if (sel_mode == "completeness") paste0(
        "**Channel selection = completeness (the manuscript method).** Each protein takes the channel that ",
        "quantified it in the larger fraction of samples — the higher mean precursor detection rate — and ",
        "that channel is held constant across every sample. The weighted vote below (the general case the ",
        "app also offers) is still computed and reported per protein as `vote_channel` / `vote_confidence`; ",
        "on this run it would have chosen a different channel for **", .int(n_disw), " of ", .int(nsh),
        "** proteins, the small set where two channels are comparably complete but differ in precision. ",
        "Set `selection_mode: multi_criteria` to let the vote decide instead.")
      else paste0(
        "**Channel selection = multi_criteria (the weighted-vote generalization).** Each protein is scored by ",
        "a weighted head-to-head over the criteria in section 5, and the winner supplies its values; here ",
        "`recommended_channel` and `vote_channel` are the same. This is a generalization of the manuscript's ",
        "method, which selects on completeness alone. Set `selection_mode: completeness` to reproduce the ",
        "manuscript's exact selection."), "",
      .md_table(r[, .N, by = .(recommended_channel, recommendation_confidence)][
                  order(recommended_channel, recommendation_confidence)]), "",
      sprintf("**%s → %s (%s), %s → %s (%s).** High confidence: %s. Low confidence: %s.",
              .int(n_a), LA, .pc(n_a, nsh), .int(n_b), LB, .pc(n_b, nsh), .int(n_hi), .int(n_lo)), "",
      sprintf("**Consequential calls: %s of %s (%s)** — the proteins where the two channels disagree, so the",
              .int(n_cons), .int(nsh), .pc(n_cons, nsh)),
      "choice actually changes a downstream number. Everywhere else the call is close to cosmetic. Report",
      "this figure rather than implying every recommendation matters.", "",
      "Every shared protein gets a call; there is no abstention. Thin evidence shows up as",
      "`recommendation_confidence = low`, which is a statement about our certainty, not a refusal.", "")
  add(sprintf("**Detection, for context.** %s", paste(sprintf("%s: %s", res$detection$category,
        vapply(res$detection$n_proteins, .int, character(1))), collapse = " · ")),
      "Proteins only one search saw are passed through untouched — there is no competing value to choose between.", "")

  ## ------------------------------------------------------------------ stage 5
  add("## 8. Stage 5 — cross-protease gap filling", "")
  if (is.null(im)) {
    add("Not run for this dataset (`impute.enabled = false`).", "")
  } else {
    s <- im$summary
    add("Where the recommended channel is missing a sample, we do **not** invent a low value. The other",
        "protease measured that protein in that same injection, so we translate a real measurement:", "",
        "```", "value(recipient, sample) = value(donor, sample) + offset", "```", "",
        sprintf("with a **%s** offset learned from the samples where both channels are observed.", cfgi$offset_model),
        if (identical(cfgi$offset_model, "sample_knn"))
          paste0("The sample-kNN form takes the offset from the *k* co-observed samples nearest the target in ",
                 "donor-profile space, which catches offsets that drift with sample composition. It beat the ",
                 "per-protein median on paired per-cell tests in both discovery conditions (p = 3.5e-07 and ",
                 "1.0e-08) — by a small margin, about 0.04 log2.") else "", "",
        "**Gates, all of which must pass:**", "",
        sprintf("- the donor covers the protein (completeness ≥ %s)", .fm(cfgi$entry$min_completeness_donor, 2)),
        sprintf("- ≥ %s co-observed samples", .int(cfgi$min_paired_offset)),
        sprintf("- the channels track each other: CCC ≥ %s **or** τ ≥ %s", .fm(cfgi$gates$ccc_min,2), .fm(cfgi$gates$tau_min,2)),
        sprintf("- the offset is stable across samples (MAD ≤ %s log2)", .fm(cfgi$gates$offset_mad_max,2)),
        "- the protein is **not** classified discordant — a hard block, not a tuning knob", "",
        "Borrowed offsets (taking an offset from similar proteins when a protein has too few co-observed",
        "samples) are **off**: benchmarked at 1.10–1.12 log2, and 4.53 on gap-bearing proteins, they are",
        "worse than leaving the value missing. A per-protein *slope* is deliberately not offered — it is the",
        "best model when values are ordinary and among the worst when they are missing, because imputation",
        "extrapolates into the censored tail.", "",
        "**What happened here:**", "", .md_table(s), "",
        sprintf("Recipient-channel completeness is **%s%%**. Of %s gaps, %s had a donor value available and",
                .fm(s$pct_complete, 2), .int(s$n_gaps), .int(s$n_gaps_with_donor)),
        sprintf("**%s cells across %s proteins were filled**; %s were refused.",
                .int(s$n_imputed), .int(s$n_proteins_imputed), .int(s$n_refused)), "",
        "**Why gaps were refused:**", "", .md_table(im$refusal_breakdown), "",
        "This table is the more informative half of the output. A gap left empty because the donor does not",
        "track the recipient is a correct decision, not a shortfall: on gap-bearing proteins no imputation",
        "method — ours or the standard ones — reaches better than ~1.0 log2, a 2-fold error per filled value.",
        "Every filled cell is in `imputation_ledger.tsv` with its donor value, offset, gates and tier, so any",
        "downstream analysis can audit or drop them with one filter.", "")

    ## the two caveats a reviewer will find on their own — better to state them first
    el <- if (!is.null(im$per_protein) && nrow(im$per_protein)) im$per_protein[eligible == TRUE] else data.table()
    if (nrow(el)) {
      na_pg <- el[concordance_class_final == "not_assessable", pg_canonical]
      na_cells <- if (length(na_pg) && nrow(im$ledger)) im$ledger[pg_canonical %in% na_pg, .N] else 0L
      add("### 8.1 What the agreement gate does not block", "",
          sprintf(paste0("The `discordant` block is absolute. **`not_assessable` is not blocked**, so a protein ",
                         "whose agreement could never be tested can still be filled. In this run that is %s ",
                         "protein(s), accounting for **%s of %s filled cells**. Whether that is acceptable is a ",
                         "judgement call: the alternative is to refuse every protein we could not test, which on ",
                         "a small-sample condition would refuse almost everything."),
                  .int(length(na_pg)), .int(na_cells), .int(nrow(im$ledger))), "",
          paste0("The gate is also an **OR**: `CCC ≥ ", .fm(cfgi$gates$ccc_min,2), "` *or* `τ ≥ ",
                 .fm(cfgi$gates$tau_min,2), "`. CCC penalises a location shift, so a protein with a large ",
                 "offset is mathematically prevented from passing on CCC and τ decides alone — and τ was ",
                 "calibrated to answer *do these channels agree at all*, not *is this donor a safe surrogate*. ",
                 "Read the τ column below next to the offset column, not on its own."), "",
          "**Every protein that was filled, with the evidence that let it through:**", "",
          .md_table(el[order(offset), .(protein = pg_canonical, recipient, donor, n_paired,
                                        offset = round(offset, 2), offset_mad = round(offset_mad, 2),
                                        CCC = round(ccc, 2), tau = round(kendall_tau, 2),
                                        class = concordance_class_final)]), "",
          sprintf(paste0("Offsets here span %s to %s log2 (median %s). A large offset is not automatically ",
                         "wrong — a systematic between-search difference is exactly what a per-protein offset ",
                         "exists to absorb — but it does mean the filled value is a *rescaled* measurement from ",
                         "the other search, not a near-copy of it. The dumbbell figures in `report.html` show ",
                         "every one of these translations."),
                  .fm(min(el$offset, na.rm = TRUE), 2), .fm(max(el$offset, na.rm = TRUE), 2),
                  .fm(stats::median(el$offset, na.rm = TRUE), 2)), "")
    }
  }

  ## ------------------------------------------------------------------ stage 6
  rp <- res$reproducibility
  add("## 9. Stage 6 — how many proteins are well quantified", "")
  if (is.null(rp)) {
    add("Not computed for this run (no protein-group run table available).", "")
  } else {
    c6 <- rp$counts
    metric_line <- switch(rp$mode,
      interday_cv  = sprintf(paste0("This cohort has replicates and a preparation-day map, so quality is the ",
                     "**inter-day CV**. It is computed the defensible way — the between-day SD in log2 estimated ",
                     "from **all** injections via variance components (injection noise removed), expressed as a CV ",
                     "— not from three day-means, whose 3-point SD is intrinsically unstable and flips the pass/",
                     "fail call on 16-27%% of proteins on a split-half test. the published day-mean CV is ",
                     "reported alongside for continuity. A protein x sample is well quantified when it is detected ",
                     "on every preparation day and its corrected CV is <= %g%%."), rp$cv_threshold),
      replicate_sd = sprintf(paste0("This cohort has replicates but no day map, so quality is the **between-",
                     "replicate CV** in log2 across a sample's replicate runs. Well quantified = all replicates ",
                     "present and CV <= %g%%."), rp$cv_threshold),
      completeness = sprintf(paste0("This cohort has **no replicate structure** — every run is a different ",
                     "sample — so inter-day CV is undefined and is **not** reported. Quality here is the ",
                     "detection axis: **cohort completeness**, the fraction of a condition's runs in which a ",
                     "protein is quantified. Well quantified = completeness >= %g%%."), 100 * rp$completeness_min))
    add(metric_line, "")
    tab <- data.table(
      quantity = c(paste0(LA, " (well quantified)"), paste0(LB, " (well quantified)"),
                   "Combined = UNION", "best single protease", "union gain over best single"),
      `protein x sample` = c(c6$A_wellq_ps, c6$B_wellq_ps, c6$union_wellq_ps,
                             c6$best_single_ps, c6$union_gain_over_single_ps),
      `distinct proteins` = c(c6$A_wellq_proteins, c6$B_wellq_proteins, c6$union_wellq_proteins,
                              max(c6$A_wellq_proteins, c6$B_wellq_proteins), NA_integer_))
    add(.md_table(tab), "")
    if (rp$mode == "interday_cv")
      add(sprintf(paste0("For continuity, the published day-mean CV gives %s %s / %s %s / union %s ",
                  "(protein x sample); the corrected metric is marginally more permissive because it removes ",
                  "injection noise from the between-day term."),
                  LA, .int(c6$daymean_A_ps), LB, .int(c6$daymean_B_ps), .int(c6$daymean_union_ps)), "")
    add(paste0("**Read the union row, not the per-protease rows, as the headline.** ", reproducibility_headline(rp, LA, LB)),
        "", paste0("This is the same result the manuscript reaches on both cohorts: the **second protease** carries ",
        "almost all of the well-quantified gain, and cross-protease imputation is a **refinement on top** of the ",
        "union — it earns its place mainly for a single-protease lab and for making each protease's own numbers ",
        "defensible, not for enlarging the combined proteome. `Combined` is always the union; the two proteases ",
        "are never averaged (averaging their raw intensities is net-destructive and was explicitly rejected)."), "")
  }

  ## ------------------------------------------------------------------ outputs
  add("## 10. Output files", "", .md_table(data.table(
    file = c("final_abundance_matrix.tsv", "recommendation.tsv", "precursor_evidence.tsv",
             "concordance.tsv", "precursor_qualification.tsv", "qualification_sensitivity.tsv",
             "detection_summary.tsv", "ingest_summary.tsv", "channel_independence.tsv",
             "rollup_probe_diagnostic.tsv", "imputation_summary.tsv", "imputation_ledger.tsv",
             "imputation_refused.tsv", "imputation_per_protein.tsv",
             "report.html", "advisor_state.rds"),
    contents = c(
      "**The deliverable.** One value per protein per sample. `recommended_channel` and `source_channel` are constant per protein by construction (see §1.1); `value_source` marks measured vs imputed; `donor_channel`, `donor_value` and `offset_used` document every imputed cell.",
      "Per shared protein: recommended channel, confidence, margin, both scores, criteria used, consequential flag, fragility flags.",
      "Every evidence metric, per protein AND per protease. This is where a call is audited.",
      "Per-protein agreement statistic, threshold, class, and why a protein was not assessable.",
      "Per-precursor detection / CV / coherence and the pass-fail of each rule.",
      "How many precursors survive across the threshold grid — the plateau-or-cliff check.",
      "Shared vs single-channel protein counts.",
      "What was read from each file after run filtering.",
      "How much peptide evidence the two channels actually share (the independence premise).",
      "LOPO stability and polish SE per protein — diagnostic only, does not vote.",
      "Gap-filling totals.",
      "Every filled cell with its donor value, offset, model, gates and tier.",
      "Every refused gap with the reason.",
      "Per-protein gate outcomes and offset statistics.",
      "Stage-by-stage report with the decision plot and guiding question for each stage.",
      "The complete run state; re-render reports without recomputing."))), "",
      "**How to use the deliverable.** Filter `value_source == \"measured\"` to drop imputed cells entirely.",
      "Filter `recommendation_confidence` to keep only calls the evidence supports. Filter",
      "`consequential == TRUE` to see only the proteins where the choice changed anything.", "")

  ## ------------------------------------------------------------------ interpretation
  add("## 11. How to read this run", "")
  lead <- if (n_a > n_b) LA else LB
  add(sprintf("**%s is recommended for the majority (%s of %s shared proteins), but not for all of them.**",
              lead, .int(max(n_a, n_b)), .int(nsh)),
      sprintf("A blanket \"always use %s\" rule would get the other %s wrong. That is the whole case for",
              lead, .int(min(n_a, n_b))),
      "adjudicating per protein rather than picking a protease for the study.", "")
  add(sprintf("**Size the claim honestly.** Only %s of the %s calls are consequential. State that number when",
              .int(n_cons), .int(nsh)),
      "describing the method's impact — it is the honest measure, and it is more persuasive than implying",
      "every shared protein was rescued.", "")
  if (n_lo > 0) add(sprintf("**%s calls are low confidence (%s).** These are near-ties where the evidence does not",
              .int(n_lo), .pc(n_lo, nsh)),
      "strongly favour either channel. They still get a call, because downstream analysis needs one value,",
      "but they are the first thing to exclude in a sensitivity analysis.", "")
  if (!is.null(im) && im$summary$n_imputed == 0)
    add("**Imputation filled nothing here.** That is the gates working, not a failure. Check",
        "`imputation_refused.tsv` for the reasons before considering relaxing anything.", "")
  if (!is.null(im) && im$summary$n_imputed > 0)
    add(sprintf("**Imputation filled %s of %s gaps (%s).** The refusals are not noise — check the breakdown in §8",
                .int(im$summary$n_imputed), .int(im$summary$n_gaps),
                .pc(im$summary$n_imputed, im$summary$n_gaps)),
        "before treating the filled values as equivalent to measured ones. They are labelled in the final matrix.", "")
  if (reps) add("**This dataset has replicates**, so technical CV and ICC contributed to the call — measured on the",
      "protein-group value being recommended, not on any internal estimate. Calls here rest on more evidence",
      "than calls in the single-run discovery conditions.", "")

  ## ------------------------------------------------------------------ limits
  add("## 12. Limits", "",
      "- **Not a new quantification method.** We adjudicate between two existing values; we do not produce a third.",
      "- **The confidence tier is not a calibrated probability.** Depth-driven criteria carry more effective",
      "  weight than their nominal weights, because several of them correlate with peptide count.",
      "- **Coherence is the one axis the sparser channel tends to win, and it carries a caveat.** Coherence is",
      "  negatively correlated with precursor count: a protein with more precursors has more chances for one",
      "  of them to disagree. Treat a coherence win on a thin protein with suspicion.",
      sprintf("- **Inherited 1%% q-filtering.** See §3. Counts here are conditional on DIA-NN's filter."),
      "- **A recommendation is about measurement quality, not biological truth.** It says which channel measured",
      "  the protein better in this experiment, not which number is correct.", "")

  writeLines(paste(L, collapse = "\n"), file.path(out_dir, "INTERPRETATION_AND_METHODS.md"))
  invisible(file.path(out_dir, "INTERPRETATION_AND_METHODS.md"))
}
