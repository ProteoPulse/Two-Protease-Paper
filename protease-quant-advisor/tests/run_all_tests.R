#!/usr/bin/env Rscript
# run_all_tests.R — run every tests/test_*.R in its own R process and summarise.
#
#   Rscript tests/run_all_tests.R
#
# Each test script is self-contained and prints its own PASS/FAIL lines; a test "passes" here when its
# process exits 0. Exit status is non-zero if any test fails, so this is usable in CI.
#
# test_example_end_to_end.R and test_shiny_app.R need nothing but this repository (they use
# example_data/; the Shiny test skips without shiny/DT). The others SKIP,
# not fail, when their inputs aren't there:
#   - test_ingest.R                 needs TWO_PROTEASE_DATA = the DIA-NN search results (MassIVE MSV000103457)
#   - test_reproducibility_stage.R  needs TWO_PROTEASE_DATA as above
#   - test_manuscript_defaults.R    checks 2-4 need a completeness-mode run under results/ or $ADVISOR_RESULTS
# Read each script's header for how to satisfy it. This harness reports skips distinctly from failures.

ROOT <- normalizePath(file.path(dirname(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
tdir  <- file.path(ROOT, "tests")
tests <- sort(list.files(tdir, pattern = "^test_.*\\.R$", full.names = TRUE))
rscript <- file.path(R.home("bin"), "Rscript")

cat(sprintf("Running %d test file(s) with %s\n\n", length(tests), rscript))
res <- data.table::data.table(test = basename(tests), status = NA_integer_, note = "")
for (i in seq_along(tests)) {
  tf <- tests[i]
  cat(sprintf("======== %s ========\n", basename(tf)))
  out <- suppressWarnings(system2(rscript, shQuote(tf), stdout = TRUE, stderr = TRUE))
  st  <- attr(out, "status"); st <- if (is.null(st)) 0L else as.integer(st)
  cat(paste(out, collapse = "\n"), "\n\n")
  skipped <- any(grepl("\\bSKIP\\b", out)) && !any(grepl("\\bFAIL\\b", out)) && st == 0L
  res$status[i] <- st
  res$note[i]   <- if (st != 0L) "FAILED" else if (skipped) "passed (some checks skipped)" else "passed"
}

cat("================ SUMMARY ================\n")
for (i in seq_len(nrow(res)))
  cat(sprintf("  %-32s %s\n", res$test[i], res$note[i]))
nfail <- sum(res$status != 0L)
cat(sprintf("\n%d/%d test files passed%s\n", nrow(res) - nfail, nrow(res),
            if (nfail) sprintf(" — %d FAILED", nfail) else ""))
quit(status = if (nfail) 1L else 0L)
