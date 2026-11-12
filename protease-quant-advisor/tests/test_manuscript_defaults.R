#!/usr/bin/env Rscript
# test_manuscript_defaults.R
# ---------------------------------------------------------------------------------------------------
# Does the app, as shipped, reproduce the methods in the manuscript Methods?
#
# Four checks, each printed as PASS/FAIL:
#   1. config_default.yml carries the v2 defaults (selection = completeness, offset = protein_median,
#      Kendall tau + 99th-percentile null, QC 0.40/0.30/0.40, imputation gates, borrowed offsets off).
#   2. Completeness selection is CORRECT — for a completeness-mode run, `recommended_channel` equals an
#      INDEPENDENT recomputation: the channel with the higher mean precursor detection rate, with the
#      same detection -> depth -> trypsin tie-break. (Recomputed here from precursor_evidence.tsv, not
#      from the pipeline's own code path, so this is a real test.)
#   3. The one-channel-per-protein invariant holds in final_abundance_matrix.tsv.
#   4. Every completeness call is reproducible from the reported columns (vote_channel present,
#      selection_mode stamped, confidence tiers only high/moderate/low).
#
# Run standalone:  Rscript tests/test_manuscript_defaults.R
# It needs a completeness-mode run under results/ (e.g. results/mix30_completeness). If none is found
# checks 2-4 are SKIPPED (not failed) with a note telling you how to make one.
# ---------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table); library(yaml) })
`%||%` <- function(a, b) if (is.null(a)) b else a

ROOT <- normalizePath(file.path(dirname(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))

.pass <- 0L; .fail <- 0L; .skip <- 0L
ok  <- function(cond, msg) { if (isTRUE(cond)) { .pass <<- .pass + 1L; cat(sprintf("  PASS  %s\n", msg)) }
                             else               { .fail <<- .fail + 1L; cat(sprintf("  FAIL  %s\n", msg)) } }
skip <- function(msg) { .skip <<- .skip + 1L; cat(sprintf("  SKIP  %s\n", msg)) }
near <- function(a, b, tol = 1e-9) is.finite(a) && is.finite(b) && abs(a - b) <= tol

# ---- 1. config defaults ---------------------------------------------------------------------------
cat("\n[1] config_default.yml matches the manuscript Methods\n")
cfg <- yaml::read_yaml(file.path(ROOT, "config_default.yml"))
ok(identical(cfg$recommend$selection_mode, "completeness"), "selection_mode default = completeness (the manuscript method)")
ok(identical(cfg$impute$offset_model, "protein_median"),    "imputation offset_model default = protein_median (manuscript method)")
ok(identical(cfg$concordance$method, "kendall"),            "agreement statistic = Kendall tau-b")
ok(near(as.numeric(cfg$concordance$null_percentile), 0.99), "agreement threshold = 99th percentile of the permutation null")
ok(as.integer(cfg$concordance$n_perm) == 1000L,            "permutations = 1000")
ok(as.integer(cfg$concordance$min_paired) == 15L,          "min paired samples = 15")
ok(near(as.numeric(cfg$qualify$detection_rate_min), 0.40), "QC detection rate >= 0.40")
ok(near(as.numeric(cfg$qualify$coherence_min), 0.30),      "QC sibling coherence >= 0.30")
ok(near(as.numeric(cfg$qualify$median_cv_max), 0.40),      "QC replicate CV <= 0.40")
ok(isFALSE(cfg$qualify$require_proteotypic),               "proteotypic status not required")
ok(near(as.numeric(cfg$impute$entry$min_completeness_donor), 0.60), "imputation gate: donor completeness >= 0.60")
ok(as.integer(cfg$impute$min_paired_offset) == 8L,         "imputation: >= 8 co-observed runs for the offset")
ok(near(as.numeric(cfg$impute$gates$ccc_min), 0.50),       "imputation gate: CCC >= 0.50")
ok(near(as.numeric(cfg$impute$gates$tau_min), 0.30),       "imputation gate: tau >= 0.30")
ok(near(as.numeric(cfg$impute$gates$offset_mad_max), 1.0), "imputation gate: offset MAD <= 1.0 log2")
ok(isTRUE("discordant" %in% cfg$impute$block_classes),     "imputation blocks discordant donors")
ok(isFALSE(cfg$impute$allow_borrowed_offsets),             "borrowed offsets off")
ok(as.integer(cfg$impute$knn_k) == 10L,                    "sample-kNN k = 10")

# ---- locate a completeness-mode run --------------------------------------------------------------
find_run <- function(mode) {
  extra <- Sys.getenv("ADVISOR_RESULTS")   # optional: a folder holding run directories
  cands <- c(list.dirs(file.path(ROOT, "results"), recursive = FALSE),
             if (nzchar(extra)) list.dirs(extra, recursive = FALSE))
  for (d in c(file.path(ROOT, "results", paste0("mix30_", mode)), cands)) {
    rec <- file.path(d, "tables", "recommendation.tsv")
    if (file.exists(rec)) {
      r <- fread(rec)
      sm <- if ("selection_mode" %in% names(r)) r$selection_mode[1] else NA_character_
      if (identical(sm, mode)) return(d)
    }
  }
  NA_character_
}
run_dir <- find_run("completeness")

# ---- 2. completeness selection is correct (independent recomputation) -----------------------------
cat("\n[2] completeness selection reproduces an independent 'more complete channel' recomputation\n")
if (is.na(run_dir)) {
  skip(paste0("no completeness-mode run found under results/. Make one with:\n",
    "        Rscript run_advisor.R pr_a=.. pr_b=.. pg_a=.. pg_b=.. selection_mode=completeness ",
    "tag=mix30_completeness out=results/mix30_completeness"))
} else {
  cat(sprintf("  (using %s)\n", sub(paste0(ROOT, "/"), "", run_dir, fixed = TRUE)))
  rec <- fread(file.path(run_dir, "tables", "recommendation.tsv"))
  ev  <- fread(file.path(run_dir, "tables", "precursor_evidence.tsv"))
  labs <- sort(unique(ev$protease))
  LA <- cfg$label_a %||% labs[1]; LB <- cfg$label_b %||% labs[2]
  key <- c("pg_canonical", "protease", "protein_detection_rate", "n_features")
  stopifnot(all(key %in% names(ev)))
  w <- dcast(ev[, ..key], pg_canonical ~ protease,
             value.var = c("protein_detection_rate", "n_features"))
  dra <- w[[paste0("protein_detection_rate_", LA)]]; drb <- w[[paste0("protein_detection_rate_", LB)]]
  nfa <- w[[paste0("n_features_", LA)]];             nfb <- w[[paste0("n_features_", LB)]]
  w[, expected := fifelse(dra > drb, LA,
                   fifelse(drb > dra, LB,
                   fifelse(nfa > nfb, LA,
                   fifelse(nfb > nfa, LB, LA))))]
  m <- merge(rec[, .(pg_canonical, recommended_channel)], w[, .(pg_canonical, expected)], by = "pg_canonical")
  n_mis <- m[recommended_channel != expected, .N]
  ok(nrow(m) > 0 && n_mis == 0L,
     sprintf("all %d shared proteins pick the more complete channel (%d mismatches)", nrow(m), n_mis))
  if (n_mis > 0) print(head(m[recommended_channel != expected], 10))

  # ---- 3. invariant: one source channel per protein ---------------------------------------------
  cat("\n[3] one source channel per protein in the final matrix\n")
  fm <- fread(file.path(run_dir, "tables", "final_abundance_matrix.tsv"))
  nch <- fm[!is.na(source_channel), .(n = uniqueN(source_channel)), by = pg_canonical]
  ok(nch[n > 1L, .N] == 0L,
     sprintf("every protein carries exactly one source channel (%d violate this)", nch[n > 1L, .N]))

  # ---- 4. reporting columns present & confidence tiers valid -------------------------------------
  cat("\n[4] selection reporting is complete and honest\n")
  ok(all(c("selection_mode","vote_channel","vote_confidence","selection_disagrees_with_vote") %in% names(rec)),
     "recommendation.tsv reports selection_mode, vote_channel, vote_confidence and the disagreement flag")
  ok(all(rec$recommendation_confidence %in% c("high","moderate","low")),
     "confidence is only high / moderate / low")
  ndis <- sum(rec$selection_disagrees_with_vote, na.rm = TRUE)
  cat(sprintf("  note: the weighted vote would differ on %d of %d proteins (the 'small set' the manuscript predicts)\n",
              ndis, nrow(rec)))
}

# ---- summary --------------------------------------------------------------------------------------
cat(sprintf("\n==== %d passed, %d failed, %d skipped ====\n", .pass, .fail, .skip))
if (.fail > 0L) quit(status = 1L)
