#!/usr/bin/env Rscript
# test_example_end_to_end.R — the whole pipeline on the bundled SYNTHETIC dataset (example_data/).
#
# Self-contained: needs no MassIVE data, so it runs on a fresh clone and in CI. Three runs:
#   A. replicate design (16 samples x 3 days), manuscript defaults, imputation on
#      -> every stage completes, Stage 6 picks the inter-day CV metric, gaps are filled
#   B. as A with offset_model = sample_knn -> the alternative offset model also fills, reversibly
#   C. discovery design (each run its own sample, no replicates) -> Stage 6 falls back to completeness
# and checks the invariants the tool promises: one source channel per protein, the chosen channel is
# the more complete one, filled values are exactly donor + offset, discordant proteins are never
# filled, measured values are never altered, and every output file is written.
#
# Usage: Rscript tests/test_example_end_to_end.R

suppressPackageStartupMessages({ library(data.table); library(yaml) })
ROOT <- normalizePath(file.path(dirname(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
source(file.path(ROOT, "engine", "load_engine.R")); load_advisor_engine(ROOT)
options(advisor.progress = function(msg, frac) invisible(NULL))

EX <- file.path(ROOT, "example_data")
f <- function(ch, lvl) file.path(EX, sprintf("%s_report.%s_matrix.tsv", ch, lvl))
stopifnot("example data missing; run example_data/make_example_data.R" =
            all(file.exists(c(f("trypsin","pr"), f("hta","pr"), f("trypsin","pg"), f("hta","pg")))))

.pass <- 0L; .fail <- 0L
ok <- function(cond, msg) {
  if (isTRUE(cond)) { .pass <<- .pass + 1L; cat(sprintf("  PASS  %s\n", msg)) }
  else              { .fail <<- .fail + 1L; cat(sprintf("  FAIL  %s\n", msg)) }
}

base_cfg <- function() {
  cfg <- yaml::read_yaml(file.path(ROOT, "config_default.yml"))
  cfg$impute$enabled <- TRUE
  cfg$concordance$n_perm <- 200L          # fast; the default 1000 only sharpens the threshold
  cfg$evidence$splithalf_reps <- 50L
  cfg
}
run <- function(cfg) run_advisor(f("trypsin","pr"), f("hta","pr"), f("trypsin","pg"), f("hta","pg"), cfg)

check_common <- function(res, label) {
  fm <- res$final_matrix; rec <- res$rec
  ok(length(res$shared_pgs) > 50L, sprintf("[%s] shared proteins adjudicated (%d)", label, length(res$shared_pgs)))
  ok(all(c("trypsin", "hta") %in% rec$recommended_channel), sprintf("[%s] each protease wins some proteins", label))
  nch <- fm[!is.na(source_channel), uniqueN(source_channel), by = pg_canonical]
  ok(nch[V1 > 1L, .N] == 0L, sprintf("[%s] exactly one source channel per protein", label))
  ok(fm[!is.na(source_channel) & recommended_channel %in% rec$recommended_channel,
        all(source_channel == recommended_channel)],
     sprintf("[%s] source_channel equals recommended_channel for every cell", label))
  # completeness selection: chosen channel's mean precursor detection rate is the higher one
  ev <- dcast(res$evidence[, .(pg_canonical, protease, protein_detection_rate)],
              pg_canonical ~ protease, value.var = "protein_detection_rate")
  m <- merge(rec[, .(pg_canonical, recommended_channel)], ev, by = "pg_canonical")
  bad <- m[(recommended_channel == "trypsin" & hta > trypsin) | (recommended_channel == "hta" & trypsin > hta), .N]
  ok(bad == 0L, sprintf("[%s] completeness mode picks the more complete channel (%d violations)", label, bad))
  # measured values are the pg value of the recommended channel, untouched
  meas <- fm[value_source == "measured", .(pg_canonical, sample_id = as.character(sample_id), value, source_channel)]
  pgv  <- res$pg[, .(pg_canonical, sample_id = as.character(sample_id), source_channel = protease, log2q)]
  mm <- merge(meas, pgv, by = c("pg_canonical", "sample_id", "source_channel"))
  ok(nrow(mm) == nrow(meas) && isTRUE(all.equal(mm$value, mm$log2q)),
     sprintf("[%s] measured cells carry the unaltered PG value (%d cells)", label, nrow(meas)))
}

check_imputation <- function(res, label, model) {
  fm <- res$final_matrix; led <- res$imputation$ledger
  imp <- fm[value_source == "imputed_cross_protease"]
  ok(nrow(imp) > 0L, sprintf("[%s] imputation filled gaps (%d cells, %d proteins)", label, nrow(imp),
                             uniqueN(imp$pg_canonical)))
  ok(isTRUE(all.equal(imp$value, imp$donor_value + imp$offset_used)),
     sprintf("[%s] every filled value == donor value + offset (reversible)", label))
  ok(all(led$offset_model_used == model), sprintf("[%s] ledger records offset model '%s'", label, model))
  disc <- res$concordance[concordance_class_final == "discordant", pg_canonical]
  ok(!any(imp$pg_canonical %in% disc), sprintf("[%s] no discordant protein is filled (%d discordant)", label, length(disc)))
  ok(all(imp$donor_channel != imp$recommended_channel), sprintf("[%s] donor is always the other protease", label))
}

# ---- A. replicate design, manuscript defaults ----------------------------------------------------
cat("\n[A] replicate design, defaults (completeness selection, protein_median offset)\n")
cfgA <- base_cfg()
cfgA$ingest$run_sample_regex <- "^[0-9]+_([0-9]+)_[0-9]+$"
cfgA$reproducibility$run_day_regex <- "^([0-9]+)_.*$"
resA <- run(cfgA)
ok(isTRUE(resA$have_replicates) && length(resA$samples) == 16L, "16 biological samples with replicate runs")
check_common(resA, "A")
check_imputation(resA, "A", "protein_median")
ok(resA$concordance[concordance_class_final == "discordant", .N] > 0L, "[A] the planted discordant proteins are found")
ok(identical(resA$reproducibility$mode, "interday_cv"), "[A] Stage 6 uses inter-day CV when replicates + a day map exist")
ok(!any(grepl("skipped", resA$notes)), "[A] no stage was skipped")

out <- tempfile("advisor_example_"); dir.create(file.path(out, "tables"), recursive = TRUE)
write_outputs(resA, out, "example"); invisible(render_advisor_report(resA, out, "example")); write_methods_doc(resA, out, "example")
need <- c("report.html", "INTERPRETATION_AND_METHODS.md", "advisor_state.rds",
          file.path("tables", c("recommendation.tsv", "final_abundance_matrix.tsv", "concordance.tsv",
                                "imputation_ledger.tsv", "imputation_per_protein.tsv", "imputation_refused.tsv",
                                "reproducibility_wellquantified.tsv")))
miss <- need[!file.exists(file.path(out, need))]
ok(!length(miss), sprintf("[A] all outputs written%s", if (length(miss)) paste0(" — missing: ", paste(miss, collapse = ", ")) else ""))

# ---- B. sample-kNN offset ------------------------------------------------------------------------
cat("\n[B] replicate design, offset_model = sample_knn\n")
cfgB <- cfgA; cfgB$impute$offset_model <- "sample_knn"
resB <- run(cfgB)
check_imputation(resB, "B", "sample_knn")
ok(identical(resA$rec$recommended_channel, resB$rec$recommended_channel),
   "[B] the offset model changes fills only, never the channel choice")

# ---- C. discovery design -------------------------------------------------------------------------
cat("\n[C] discovery design (one run per sample, no replicates)\n")
cfgC <- base_cfg()
resC <- run(cfgC)
ok(!isTRUE(resC$have_replicates) && length(resC$samples) == 48L, "48 runs treated as 48 samples")
check_common(resC, "C")
ok(identical(resC$reproducibility$mode, "completeness"), "[C] Stage 6 falls back to completeness without replicates")

cat(sprintf("\n==== %d passed, %d failed ====\n", .pass, .fail))
quit(status = if (.fail > 0L) 1L else 0L)
