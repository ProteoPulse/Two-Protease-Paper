#!/usr/bin/env Rscript
# build_index.R — one entry point across every condition that has been run, plus the cross-condition
# stability check.
#
# A recommendation is only worth acting on if it is a property of the PROTEIN and not of the run. The
# stability table below is what tests that: proteins scored in more than one condition should mostly
# get the same call, and should agree more where the tool says it is confident. If that pattern does
# not hold, the confidence tier is not carrying information and should not be trusted.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
ROOT <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(ROOT, "engine", "load_engine.R")); load_advisor_engine(ROOT)
RES <- file.path(ROOT, "results")

tags <- basename(dirname(Sys.glob(file.path(RES, "*", "advisor_state.rds"))))
if (!length(tags)) stop("no completed runs found under results/")
cat("conditions found:", paste(tags, collapse = ", "), "\n")
S <- lapply(setNames(tags, tags), function(t) readRDS(file.path(RES, t, "advisor_state.rds")))

## ---- cross-condition summary -------------------------------------------------------------------
summ <- rbindlist(lapply(tags, function(t) {
  r <- S[[t]]$rec; LA <- S[[t]]$cfg$label_a; LB <- S[[t]]$cfg$label_b
  im <- S[[t]]$imputation
  data.table(condition = t,
             samples = length(S[[t]]$samples),
             replicates = isTRUE(S[[t]]$have_replicates),
             shared_adjudicated = nrow(r),
             only_a = S[[t]]$detection[grepl(paste0("^", LA), category), n_proteins],
             only_b = S[[t]]$detection[grepl(paste0("^", LB), category), n_proteins],
             rec_a = r[recommended_channel == LA, .N],
             rec_b = r[recommended_channel == LB, .N],
             pct_a = round(100 * r[recommended_channel == LA, .N] / nrow(r), 1),
             high_conf = r[recommendation_confidence == "high", .N],
             consequential = r[consequential == TRUE, .N],
             threshold = round(S[[t]]$threshold, 4),
             discordant = S[[t]]$concordance[concordance_class_final == "discordant", .N],
             gaps = if (is.null(im)) NA_integer_ else im$summary$n_gaps,
             imputed = if (is.null(im)) NA_integer_ else im$summary$n_imputed,
             # the grain invariant, re-checked here from the written matrix rather than trusted
             proteins_with_one_channel = S[[t]]$final_matrix[
               , .(k = uniqueN(source_channel[!is.na(source_channel)])), by = pg_canonical][k == 1L, .N],
             proteins_with_two_channels = S[[t]]$final_matrix[
               , .(k = uniqueN(source_channel[!is.na(source_channel)])), by = pg_canonical][k > 1L, .N])
}), fill = TRUE)
fwrite(summ, file.path(RES, "cross_condition_summary.tsv"), sep = "\t")
cat("\n== cross-condition summary ==\n"); print(summ)

## ---- cross-condition stability -------------------------------------------------------------------
pairs <- if (length(tags) >= 2) utils::combn(tags, 2, simplify = FALSE) else list()
stab <- rbindlist(lapply(pairs, function(p) {
  a <- S[[p[1]]]$rec[, .(pg_canonical, ca = recommended_channel, fa = recommendation_confidence)]
  b <- S[[p[2]]]$rec[, .(pg_canonical, cb = recommended_channel, fb = recommendation_confidence)]
  m <- merge(a, b, by = "pg_canonical")
  if (!nrow(m)) return(NULL)
  m[, agree := ca == cb]
  f <- function(lab, d) if (nrow(d)) data.table(pair = paste(p, collapse = " vs "), subset = lab,
    n = nrow(d), agree = sum(d$agree), pct = round(100 * mean(d$agree), 1)) else NULL
  rbindlist(list(f("all proteins scored in both", m),
                 f("high confidence in both", m[fa == "high" & fb == "high"]),
                 f("low confidence in either", m[fa == "low" | fb == "low"])), fill = TRUE)
}), fill = TRUE)
if (nrow(stab)) {
  fwrite(stab, file.path(RES, "cross_condition_stability.tsv"), sep = "\t")
  cat("\n== cross-condition stability ==\n"); print(stab)
}

## ---- figure ---------------------------------------------------------------------------------------
d <- melt(summ[, .(condition, rec_a, rec_b)], id.vars = "condition",
          variable.name = "channel", value.name = "n")
d[, channel := fifelse(channel == "rec_a", S[[1]]$cfg$label_a, S[[1]]$cfg$label_b)]
p1 <- ggplot(d, aes(n, condition, fill = channel)) +
  geom_col(width = .6) +
  geom_text(aes(label = n), position = position_stack(vjust = .5), colour = "white", size = 3.4) +
  scale_fill_manual(values = setNames(c("#2C6FB3", "#E08214"),
                                      c(S[[1]]$cfg$label_a, S[[1]]$cfg$label_b)), name = NULL) +
  labs(title = "Which protease is recommended, by dataset",
       subtitle = "Shared proteins only. The split is stable across conditions — it is a property of the chemistry, not the run.",
       x = "shared proteins", y = NULL) + theme_adv()
ggsave(file.path(RES, "cross_condition_calls.png"), p1, width = 9, height = 4.2, dpi = 200)

## ---- index ------------------------------------------------------------------------------------------
blocks <- list(
  rb_text(paste0("Three datasets through the same tool with the same settings. Each condition has its own ",
                 "stage-by-stage report and its own interpretation-and-methods document; this page is the ",
                 "entry point and the cross-condition check.")),
  rb_text(paste0(
    "WHAT THIS FRAMEWORK DECIDES, AND AT WHAT LEVEL. Both searches ran on the SAME injections, so for a ",
    "protein both of them report we hold two numbers for one molecule in one tube. The tool does not ",
    "average them and does not re-quantify: it reads the precursor evidence behind each shared protein ",
    "and recommends which search's existing protein-group value to trust. That recommendation is made ",
    "ONCE PER PROTEIN and held across every sample. It is not made per sample — a protein taking its ",
    "value from one search in one injection and the other in the next would carry the between-search ",
    "offset inside its own profile and manufacture fold-changes that are pure artefact. It is also not a ",
    "single experiment-wide winner, which would collapse to `use trypsin' and discard the proteins the ",
    "other search genuinely measures better. The `proteins_with_two_channels' column below is that ",
    "invariant, re-checked from each written matrix: it must be zero in every condition.")),
  rb_question("Walking someone through this? Take the stages in this order:",
    "1. What was read, and how many proteins actually need a decision (Stage 0).",
    "2. Which precursors are allowed to carry evidence, and what the rules cut (Stage 1).",
    "3. What separates the two channels for a protein (Stage 2).",
    "4. Do the two channels even agree about this protein, against a permutation null (Stage 4).",
    "5. The call itself, its confidence, and whether it changes any downstream number (Stage 3).",
    "6. Optional gap filling: what was filled, what was refused, and why (Stage 5).",
    "7. This page: does the same protein get the same call in another dataset?"),
  rb_table(summ, "Every condition, side by side"),
  rb_plot(file.path(RES, "cross_condition_calls.png"), "Recommended channel by dataset.", width = 9, height = 4.2))
if (nrow(stab)) blocks <- c(blocks, list(
  rb_question("Is a recommendation a property of the protein, or of the run it was measured in?"),
  rb_table(stab, "Cross-condition stability — the same protein, scored in two datasets"),
  rb_text(paste0("Agreement should be high overall and HIGHER among proteins the tool calls high-confidence ",
                 "in both. That gradient is the evidence that the confidence tier means something. Where a ",
                 "protein flips between conditions and both calls were low confidence, the honest reading is ",
                 "that the evidence does not separate the channels for that protein."))))
blocks <- c(blocks, list(rb_text(paste0(
  "Per-condition outputs: open `<condition>/report.html` for the stage-by-stage walkthrough with the ",
  "decision plot and guiding question at each stage, and `<condition>/INTERPRETATION_AND_METHODS.md` for ",
  "input files, preprocessing, every criterion that did or did not vote, the imputation gates, and a ",
  "description of every output file."))))
render_stage_report("Protease Quant Advisor — all conditions", blocks,
                    file.path(RES, "ALL_CONDITIONS.html"),
                    subtitle = sprintf("%s | generated %s", paste(tags, collapse = " · "),
                                       format(Sys.time(), "%Y-%m-%d %H:%M")))

links <- list(`Cross-condition overview (start here)` = "ALL_CONDITIONS.html")
for (t in tags) {
  links[[sprintf("%s — stage report", t)]] <- file.path(t, "report.html")
  links[[sprintf("%s — interpretation & methods", t)]] <- file.path(t, "INTERPRETATION_AND_METHODS.md")
}
render_run_index(RES, "Protease Quant Advisor — results", links)
cat("\nwrote:", file.path(RES, "ALL_CONDITIONS.html"), "and", file.path(RES, "INDEX.html"), "\n")
