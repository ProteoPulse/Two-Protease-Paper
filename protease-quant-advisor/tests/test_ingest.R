#!/usr/bin/env Rscript
# test_ingest.R — ingest of the real Mix30 discovery matrices (MassIVE MSV000103457).
#
# Needs TWO_PROTEASE_DATA pointing at the DIA-NN search results (layout: ../../figure-scripts/DATA_LAYOUT.md).
# Without it the test prints SKIP and exits 0. With it, it asserts the numbers reported in the manuscript:
# 54 Mix30 runs in both searches, 392 protein groups quantified by both proteases.
#
# Usage: TWO_PROTEASE_DATA=/path/to/search_results Rscript tests/test_ingest.R

suppressPackageStartupMessages({ library(data.table) })
ROOT <- normalizePath(file.path(dirname(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
S <- Sys.getenv("TWO_PROTEASE_DATA")
PR_A <- file.path(S, "Trypsin_Library_Search",  "report.pr_matrix.tsv")
PR_B <- file.path(S, "Krakatoa_Library_Search", "report.pr_matrix.tsv")
PG_A <- file.path(S, "Trypsin_Library_Search",  "trypsin_library_search_report.pg_matrix.tsv")
PG_B <- file.path(S, "Krakatoa_Library_Search", "report.pg_matrix.tsv")
if (!nzchar(S) || !all(file.exists(c(PR_A, PR_B, PG_A, PG_B)))) {
  cat("  SKIP  real-data ingest: set TWO_PROTEASE_DATA to the DIA-NN search results\n"); quit(status = 0)
}
source(file.path(ROOT, "engine", "load_engine.R")); load_advisor_engine(ROOT)

fails <- 0L
ok <- function(cond, msg) { cat(sprintf("  %s  %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg))
                            if (!isTRUE(cond)) fails <<- fails + 1L }

t0 <- Sys.time()
prA <- load_pr_matrix(PR_A, "trypsin", run_include = "Plasma_Mix", run_exclude = "20min")
prB <- load_pr_matrix(PR_B, "hta",     run_include = "Plasma_Mix", run_exclude = "20min")
pgA <- load_pg_matrix(PG_A, "trypsin", run_include = "Plasma_Mix", run_exclude = "20min")
pgB <- load_pg_matrix(PG_B, "hta",     run_include = "Plasma_Mix", run_exclude = "20min")
cat(sprintf("ingest took %.1f s\n\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
prec <- rbind(prA, prB); pg <- rbind(pgA, pgB)
print(ingest_summary(prec, pg, "trypsin", "hta"))

n_ch <- pg[, uniqueN(protease), by = pg_canonical]
sA <- unique(pgA$sample_id); sB <- unique(pgB$sample_id)
cat("\n")
ok(length(sA) == 54L && length(sB) == 54L,         sprintf("54 Mix30 runs in each search (%d / %d)", length(sA), length(sB)))
ok(setequal(sA, sB),                                "the two searches cover the same samples")
ok(n_ch[V1 == 2L, .N] == 392L,                      sprintf("392 protein groups in both searches (got %d)", n_ch[V1 == 2L, .N]))
ok(!any(grepl("[/\\\\]|\\.raw$", prA$run)),          "run names reduced to basenames without extension")
ok(all(is.finite(prec$log2q)),                      "precursor values are finite log2")
ok(setequal(unique(prA$sample_id), sA),             "precursor and protein matrices agree on samples")
cat(sprintf("\n%s\n", if (fails) sprintf("%d check(s) FAILED", fails) else "ALL TESTS PASSED"))
quit(status = if (fails) 1L else 0L)
