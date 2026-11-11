# test_reproducibility_stage.R — the Stage 6 reproducibility module must
#   (1) reproduce Panel A exactly in interday_cv mode (the published linear day-mean CV),
#   (2) fall back to completeness mode on discovery data (no replicates),
# so the app's well-quantified numbers cannot silently drift from the manuscript.
#
# Needs TWO_PROTEASE_DATA pointing at the DIA-NN search results (MassIVE MSV000103457; layout in
# ../../figure-scripts/DATA_LAYOUT.md). Without it the test prints SKIP and exits 0.
#
# Usage: TWO_PROTEASE_DATA=/path/to/search_results Rscript protease-quant-advisor/tests/test_reproducibility_stage.R

suppressPackageStartupMessages({ library(data.table) })
ROOT <- normalizePath(file.path(dirname(sub("^--file=", "",
  commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])), ".."))
DATA <- Sys.getenv("TWO_PROTEASE_DATA")
REPRO <- file.path(DATA, "Reproducibility_Search")
DISC_TR <- file.path(DATA, "Trypsin_Library_Search", "trypsin_library_search_report.pg_matrix.tsv")
DISC_KR <- file.path(DATA, "Krakatoa_Library_Search", "report.pg_matrix.tsv")
need <- c(file.path(REPRO, c("Trypsin Library", "Krakatoa Library"), "report.pg_matrix.tsv"), DISC_TR, DISC_KR)
if (!nzchar(DATA) || !all(file.exists(need))) {
  cat("  SKIP  Stage 6 on real data: set TWO_PROTEASE_DATA to the DIA-NN search results\n"); quit(status = 0)
}
source(file.path(ROOT, "engine", "reproducibility.R"))

ANN <- c("Genes","Protein.Group","Protein.Ids","Protein.Names",
         "First.Protein.Description","N.Sequences","N.Proteotypic.Sequences")

# build pg_run (run grain) the way the app's ingest would, from a pg_matrix
load_pg_run <- function(path, protease, sample_tok, day_tok = NULL) {
  dt <- fread(path, sep = "\t", header = TRUE, showProgress = FALSE)
  rc <- setdiff(names(dt), ANN)
  L <- melt(dt, id.vars = "Protein.Group", measure.vars = rc,
            variable.name = "Sample", value.name = "Intensity", variable.factor = FALSE)
  L[, Intensity := suppressWarnings(as.numeric(Intensity))]
  L[!is.finite(Intensity) | Intensity <= 0, Intensity := NA_real_]
  L <- L[, .(Intensity = if (all(is.na(Intensity))) NA_real_ else sum(Intensity, na.rm = TRUE)),
         by = .(pg_canonical = Protein.Group, run = Sample)]
  L <- L[!is.na(pg_canonical) & pg_canonical != "" & is.finite(Intensity)]
  b <- sub("\\.raw$", "", basename(gsub("\\\\", "/", L$run)), ignore.case = TRUE)
  L[, run := b]
  p <- tstrsplit(b, "_", fixed = TRUE)
  L[, sample_id := p[[sample_tok]]]
  L[, protease := protease][, log2q := log2(Intensity)][, Intensity := NULL][]
}

fail <- function(msg) { cat("FAIL:", msg, "\n"); quit(status = 1) }

# ============================ (1) repro: interday_cv ==========================
tr <- load_pg_run(file.path(REPRO, "Trypsin Library",  "report.pg_matrix.tsv"), "trypsin", 2)
kr <- load_pg_run(file.path(REPRO, "Krakatoa Library", "report.pg_matrix.tsv"), "hta", 2)
pg_run <- rbind(tr, kr)

cfg_rep <- list(run_day_regex = "^([0-9]+)_.*$", cv_threshold = 20,
                completeness_min = 0.50)
rs <- reproducibility_stage(pg_run, imp_mat = NULL, "trypsin", "hta", cfg_rep, have_reps = TRUE)
cat("mode:", rs$mode, "\n")
if (rs$mode != "interday_cv") fail("expected interday_cv mode on the reproducibility cohort")

# per-sample day-mean CV counts must equal Panel A
so <- as.character(sort(as.integer(unique(pg_run$sample_id))))
per_sample <- function(pr) {
  w <- rs$per_protein[protease == pr & well_daymean == TRUE, .N, by = sample_id]
  v <- w$N[match(so, w$sample_id)]; v[is.na(v)] <- 0L; as.integer(v)
}
PANEL_A_TRY <- c(120L,107L,127L,123L,123L,123L,146L,131L)
PANEL_A_HTA <- c(92L,80L,88L,89L,108L,101L,91L,82L)
bt <- per_sample("trypsin"); bh <- per_sample("hta")
cat(sprintf("  Trypsin  day-mean CV counts: %s  %s\n", paste(bt, collapse=" "),
            if (identical(bt, PANEL_A_TRY)) "MATCH" else "MISMATCH"))
cat(sprintf("  Krakatoa day-mean CV counts: %s  %s\n", paste(bh, collapse=" "),
            if (identical(bh, PANEL_A_HTA)) "MATCH" else "MISMATCH"))
if (!identical(bt, PANEL_A_TRY) || !identical(bh, PANEL_A_HTA))
  fail("Stage 6 interday_cv mode did NOT reproduce Panel A")

# corrected CV totals (protein x sample) and published-CV continuity
cat(sprintf("  corrected-CV well (protein x sample): trypsin %d  hta %d  union %d (+%d over best single)\n",
            rs$counts$A_wellq_ps, rs$counts$B_wellq_ps, rs$counts$union_wellq_ps,
            rs$counts$union_gain_over_single_ps))
cat(sprintf("  Published day-mean CV totals:      trypsin %d  hta %d  union %d\n",
            rs$counts$daymean_A_ps, rs$counts$daymean_B_ps, rs$counts$daymean_union_ps))
if (rs$counts$daymean_A_ps != 1000L || rs$counts$daymean_B_ps != 731L)
  fail("Published day-mean CV protein x sample totals drifted from 1000 (trypsin) / 731 (hta)")
cat("  ", reproducibility_headline(rs, "Trypsin", "Krakatoa"), "\n")

# ============================ (2) discovery: completeness ====================
# for discovery, treat one condition (Mix30) — each run its own sample, no reps
dtr <- load_pg_run(DISC_TR, "trypsin", 1); dkr <- load_pg_run(DISC_KR, "hta", 1)
mix30 <- function(d) d[grepl("Plasma_Mix", run) & !grepl("20min", run)]
dtr <- mix30(dtr); dkr <- mix30(dkr)
# every run is its own sample -> collapse sample_id to run so have_reps = FALSE
dtr[, sample_id := run]; dkr[, sample_id := run]
pg_run_d <- rbind(dtr, dkr)
have_reps_d <- pg_run_d[, uniqueN(run)] > pg_run_d[, uniqueN(sample_id)]
rs_d <- reproducibility_stage(pg_run_d, imp_mat = NULL, "trypsin", "hta",
                              list(cv_threshold = 20, completeness_min = 0.50), have_reps_d)
cat("\ndiscovery mode:", rs_d$mode, "(have_reps =", have_reps_d, ")\n")
if (rs_d$mode != "completeness") fail("expected completeness mode on discovery (no replicates)")
# NOTE: with sample_id==run, completeness per run is 0/1; the meaningful cohort
# completeness is computed by figure-scripts/06_supplementary/discovery_imputation.R. Here we only assert
# the stage SELECTS the detection axis rather than fabricating a CV.
cat("  discovery correctly routed to the detection axis (completeness), CV suppressed.\n")

cat("\nALL TESTS PASSED\n")
