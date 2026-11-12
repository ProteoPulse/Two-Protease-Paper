# outputs.R — the tables the tool writes, and the HTML report that explains them.
#
# Design rule: every table a reader might want to act on is a plain TSV, and every claim in the HTML
# is backed by one of them. The report follows the framework's stage order and puts the decision
# question next to the decision plot, so each stage can be validated on its own.

suppressPackageStartupMessages({ library(data.table) })

advisor_tables <- function(res) {
  LA <- res$cfg$label_a; LB <- res$cfg$label_b
  keep_rec <- intersect(c("pg_canonical","recommended_channel","recommendation_confidence",
    "recommendation_margin","selection_mode","vote_channel","vote_confidence","vote_margin",
    "selection_disagrees_with_vote","score_a","score_b","n_usable_criteria","consequential",
    "concordance_class_final","flag_low_independence","rationale",
    paste0("fragile_rollup_", c(LA, LB))), names(res$rec))
  out <- list(
    `recommendation` = res$rec[, ..keep_rec],
    `final_abundance_matrix` = res$final_matrix,
    `precursor_evidence` = res$evidence,
    `precursor_qualification` = res$qual,
    `concordance` = res$concordance,
    `detection_summary` = res$detection,
    `ingest_summary` = res$ingest_summary)
  if (!is.null(res$qual_sweep))  out$qualification_sensitivity <- res$qual_sweep
  if (!is.null(res$indep))       out$channel_independence <- res$indep
  if (!is.null(res$probe_flags)) out$rollup_probe_diagnostic <- res$probe_flags
  if (!is.null(res$imputation)) {
    out$imputation_summary   <- res$imputation$summary
    out$imputation_ledger    <- res$imputation$ledger
    out$imputation_refused   <- res$imputation$refused
    out$imputation_per_protein <- res$imputation$per_protein
  }
  if (!is.null(res$reproducibility)) {
    out$reproducibility_wellquantified <- res$reproducibility$counts
    out$reproducibility_per_protein    <- res$reproducibility$per_protein
  }
  out
}

write_outputs <- function(res, out_dir, tag = "run") {
  td <- file.path(out_dir, "tables"); dir.create(td, recursive = TRUE, showWarnings = FALSE)
  tb <- advisor_tables(res)
  for (nm in names(tb)) if (!is.null(tb[[nm]]) && nrow(tb[[nm]]))
    fwrite(tb[[nm]], file.path(td, paste0(nm, ".tsv")), sep = "\t")
  saveRDS(res, file.path(out_dir, "advisor_state.rds"))
  invisible(names(tb))
}

# Headline counts used by both the app banner and the report.
advisor_headline <- function(res) {
  LA <- res$cfg$label_a; LB <- res$cfg$label_b; r <- res$rec
  data.table(
    metric = c("samples", "shared proteins adjudicated",
               paste0("recommended -> ", LA), paste0("recommended -> ", LB),
               "high confidence", "consequential (channels disagree)",
               "agreement threshold", "cells imputed"),
    value = c(length(res$samples), nrow(r),
              r[recommended_channel == LA, .N], r[recommended_channel == LB, .N],
              r[recommendation_confidence == "high", .N],
              r[consequential == TRUE, .N],
              sprintf("%.3f", res$threshold),
              if (is.null(res$imputation)) "not run" else as.character(res$imputation$summary$n_imputed)))
}

render_advisor_report <- function(res, out_dir, tag = "run") {
  LA <- res$cfg$label_a; LB <- res$cfg$label_b
  rd <- file.path(out_dir); dir.create(rd, recursive = TRUE, showWarnings = FALSE)
  # Build each plot EAGERLY into a PNG. A ggplot object only fails when it is drawn, so constructing
  # it inside tryCatch would not catch anything; rendering it here does. One unbuildable figure then
  # costs its own panel and nothing else — the report and every table still come out.
  P <- function(f) tryCatch(f(res), error = function(e) NULL)
  .tmpdir <- file.path(tempdir(), paste0("advfig_", as.integer(Sys.time())))
  dir.create(.tmpdir, recursive = TRUE, showWarnings = FALSE)
  .i <- 0L
  # Subtitles carry the interpretation, so a subtitle that runs off the right edge loses exactly the
  # part a reader needs. Wrap it to the figure's width before saving rather than trimming the text.
  .wrap_sub <- function(obj, w) {
    sub <- tryCatch(obj$labels$subtitle, error = function(e) NULL)
    if (is.null(sub) || !is.character(sub) || !nzchar(sub)) return(obj)
    obj$labels$subtitle <- paste(strwrap(sub, width = max(60, floor(w * 14.5))), collapse = "\n")
    obj
  }
  add <- function(blocks, obj, cap, w = 9, h = 5) {
    if (is.null(obj)) return(blocks)
    obj <- .wrap_sub(obj, w)
    .i <<- .i + 1L
    f <- file.path(.tmpdir, sprintf("p%02d.png", .i))
    okp <- tryCatch({ ggplot2::ggsave(f, obj, width = w, height = h, dpi = 110); TRUE },
                    error = function(e) { message("figure skipped (", cap, "): ", conditionMessage(e)); FALSE })
    if (!okp) return(blocks)
    c(blocks, list(rb_plot(f, cap, width = w, height = h)))
  }

  b <- list(
    rb_text(paste0(
      "This tool does not produce a new protein quantity. It reads the precursor evidence behind each ",
      "shared protein and recommends WHICH of the two searches' existing protein-group values to trust, ",
      "protein by protein, with a confidence tier. Proteins only one search saw are passed through ",
      "untouched. Gap filling is optional and off unless you asked for it.")),
    rb_table(advisor_headline(res), "Headline"),
    rb_table(res$ingest_summary, "What was actually read"),
    rb_text(res$notes[1]))
  # The single most common misreading of this framework is the grain of its decision. State it before
  # any number is shown, and back it with the invariant the pipeline actually enforces.
  .fm <- res$final_matrix
  .nch <- .fm[, .(n = uniqueN(source_channel[!is.na(source_channel)])), by = pg_canonical]
  b <- c(b, list(
    rb_text(paste0(
      "WHAT THE DECISION IS, AND AT WHAT LEVEL. The recommendation is made ONCE PER PROTEIN and held ",
      "across every sample in the experiment. It is not made per sample, and it is not a single ",
      "experiment-wide winner. Read the three levels in order: (1) per protein per sample — REJECTED. ",
      "Letting a protein come from trypsin in one injection and HTA in the next puts a channel offset ",
      "inside the protein's own profile and manufactures fold-changes that are pure search artefact. ",
      "(2) per protein, fixed across samples — WHAT WE DO. Each shared protein is adjudicated once on ",
      "its precursor evidence, and that channel supplies all of its values, so every comparison within ",
      "a protein stays inside one measurement system. (3) one channel for the whole experiment — too ",
      "blunt. It would be `use ", res$cfg$label_a, "', and would discard the proteins the other search ",
      "genuinely measures better.")),
    rb_table(data.table(
      check = c("proteins in the final matrix",
                "proteins carrying exactly one source channel across all samples",
                "proteins carrying more than one",
                "imputed cells (channel unchanged; donor translated onto the recommended scale)"),
      value = c(format(nrow(.nch), big.mark = ","),
                format(.nch[n == 1L, .N], big.mark = ","),
                format(.nch[n > 1L, .N], big.mark = ","),
                format(.fm[value_source == "imputed_cross_protease", .N], big.mark = ","))),
      "Grain check — enforced by the pipeline, not assumed"),
    rb_text(paste0(
      "Imputation does not break this. When the recommended channel has no value for a sample, the ",
      "other protease's measurement of the same protein in the same injection is translated onto the ",
      "recommended channel's scale by that protein's own offset. The number that is written still ",
      "lives in the recommended channel's frame, so `source_channel` stays constant for the protein; ",
      "the protease that supplied the evidence is recorded separately in `donor_channel`."))))
  if (length(res$notes) > 1) b <- c(b, list(rb_question(res$notes[-1])))

  b <- c(b, list(rb_question(
    "Do the two searches cover the same samples, and is the shared set big enough to matter?")))
  b <- add(b, P(p_detection), "Shared proteins are the only ones with a decision to make.")

  b <- c(b, list(rb_text("STAGE 1 — QUALIFICATION. Our own rule, applied on top of DIA-NN's 1% q-filtering."),
                 rb_question("Does this threshold pair hold quality without discarding most of the data?",
                             "Is the chosen setting on a plateau, or on a cliff where a small change costs many precursors?")))
  b <- add(b, P(p_qual_sweep), "Qualifying precursors across the threshold grid; the ring is the setting in use.")
  b <- add(b, P(p_qc_metrics), "Each qualification rule drawn on the distribution it cuts.", w = 10, h = 4)
  b <- add(b, P(p_qc_traces), "A representative protein: the precursors the rule keeps, and the ones it rejects.", w = 9.5, h = 5)
  b <- c(b, list(rb_table(res$qual[, .(n_precursors = .N,
                                       qualifying = sum(qualifies),
                                       pct = round(100*mean(qualifies), 1))],
                          "Stage 1 outcome")))

  b <- c(b, list(rb_text("STAGE 2 — EVIDENCE. Completeness and coherence, measured per protein and per protease."),
                 rb_question("Which channel is better supported for this protein, and by how much?")))
  b <- add(b, P(p_evidence), "Distribution of each evidence family across the shared set.")
  if (!is.null(res$probe_flags)) {
    b <- c(b, list(rb_text(paste0(
      "The median polish below is an INTERNAL QC PROBE. Its output is never published as a quantity, and ",
      "on the Aug-24 ablation it did not measurably improve the recommendation (Spearman 0.91 with peptide ",
      "count; McNemar p = 0.48). It is reported here as a per-protein diagnostic and, by default, is kept ",
      "out of the score."))))
    b <- add(b, P(p_probe), "Proteins whose value rests on a single peptide.")
  }

  b <- c(b, list(rb_text("STAGE 4 — AGREEMENT. Computed from the protein-group values, with a permutation null."),
                 rb_question("Is this threshold defensible as a null percentile rather than a chosen cut?")))
  b <- add(b, P(p_concordance), "Observed statistic against the permutation null.")

  .sel <- res$rec$selection_mode[1] %||% "completeness"
  .disw <- if ("selection_disagrees_with_vote" %in% names(res$rec))
    sum(res$rec$selection_disagrees_with_vote, na.rm = TRUE) else NA_integer_
  b <- c(b, list(rb_text(paste0(
      "STAGE 3 — THE RECOMMENDATION. Channel selection = ",
      if (.sel == "completeness") paste0(
        "COMPLETENESS, the manuscript method: each protein takes the channel with the higher mean ",
        "precursor detection rate, held across all samples. The weighted vote is still computed and ",
        "reported as `vote_channel`; on this run it would have changed ", .disw, " of ", nrow(res$rec),
        " calls. Switch `selection_mode` to `multi_criteria` to let the weighted vote decide instead.")
      else paste0(
        "MULTI_CRITERIA, the weighted-vote generalization the app offers (completeness, coherence and ",
        "split-half reliability weight 2; peptide depth weight 1; with replicates technical CV weight 2 ",
        "and replicate ICC weight 1). Switch `selection_mode` to `completeness` to reproduce the ",
        "manuscript's exact selection."))),
                 rb_question("Is the split between the two channels plausible for this chemistry?",
                             "How many of these calls actually change a downstream number?")))
  b <- add(b, P(p_recommendation), "Calls by channel and confidence tier.")
  b <- add(b, P(p_consequential), "Only where the channels disagree does the choice change the answer.")

  if (!is.null(res$imputation)) {
    .g <- res$cfg$impute$gates
    .pp <- res$imputation$per_protein
    .el <- if (!is.null(.pp) && nrow(.pp)) .pp[eligible == TRUE] else data.table()
    .na_cells <- if (nrow(.el) && !is.null(res$imputation$ledger) && nrow(res$imputation$ledger))
      res$imputation$ledger[pg_canonical %in% .el[concordance_class_final == "not_assessable", pg_canonical], .N] else 0L
    b <- c(b, list(rb_text(paste0(
      "STAGE 5 — OPTIONAL GAP FILLING. Gaps in the recommended channel are filled from the other ",
      "protease's measurement of the same protein in the same injection, using a ",
      res$cfg$impute$offset_model, " offset. The channel does not change: the donor's value is ",
      "translated onto the recommended channel's scale, and the protein keeps its recommendation.")),
      rb_text(paste0(
        "WHAT THE AGREEMENT GATE DOES AND DOES NOT BLOCK — read this before quoting the filled values. ",
        "A protein classified `discordant` is never filled; that block is absolute. But `not_assessable` ",
        "is NOT blocked, so a protein whose agreement could never be tested can still be used. In this ",
        "run that accounts for ", .na_cells, " of ", nrow(res$imputation$ledger), " filled cells. ",
        "Separately, the gate is `CCC >= ", .g$ccc_min, " OR Kendall tau >= ", .g$tau_min, "`. CCC is ",
        "mathematically depressed by a large location shift, so for a protein with a big offset CCC can ",
        "essentially never pass and tau decides alone — and tau is being asked a question it was not ",
        "calibrated for. The table below gives every filled protein's offset, CCC, tau and class so a ",
        "reader can judge each fill rather than take the total on trust.")),
      rb_question("Are the refusals concentrated on proteins where the two channels genuinely diverge?",
                  "For each filled protein: is the offset small enough that the donor is a surrogate, or large enough that it is a different measurement being rescaled?",
                  "Would relaxing the agreement gate buy accuracy, or just more numbers?")))
    b <- add(b, P(p_imputation), "Every gap, filled or refused, with the reason.")
    b <- add(b, P(p_impute_dumbbell),
             "Every imputed cell: the donor's own value, and the value written after translation.", w = 9, h = 6)
    b <- add(b, P(p_impute_protein_dumbbell),
             "Where the filled values land relative to each protein's own measured range.", w = 9, h = 6.5)
    b <- c(b, list(rb_table(res$imputation$summary, "Gap-filling outcome"),
                   rb_table(res$imputation$refusal_breakdown, "Why gaps were refused")))
    if (nrow(.el)) b <- c(b, list(rb_table(
      .el[order(offset), .(protein = pg_canonical, recipient, donor, n_paired,
                           offset = round(offset, 2), offset_mad = round(offset_mad, 2),
                           CCC = round(ccc, 2), tau = round(kendall_tau, 2),
                           agreement_class = concordance_class_final,
                           completeness_recipient = round(completeness_recipient, 2))],
      "Every protein that was filled, with the evidence that let it through")))
  }

  b <- c(b, list(rb_text(paste0(
    "OUTPUT. `final_abundance_matrix.tsv` carries one value per protein per sample. `recommended_channel` ",
    "and `source_channel` are constant for a protein by construction; `value_source` says whether the cell ",
    "was measured or imputed, and for imputed cells `donor_channel`, `donor_value` and `offset_used` show ",
    "exactly how the number was produced. Filter on `value_source` to drop imputed cells, or on ",
    "`recommendation_confidence` to keep only confident calls. `recommendation.tsv` is the decision itself: ",
    "one row per shared protein, which is the grain at which this framework commits."))))

  render_stage_report(sprintf("Protease Quant Advisor — %s", tag), b,
                      file.path(rd, "report.html"),
                      subtitle = sprintf("%d samples | %d shared proteins | threshold %.3f | %s",
                                         length(res$samples), nrow(res$rec), res$threshold,
                                         format(Sys.time(), "%Y-%m-%d %H:%M")))
}
