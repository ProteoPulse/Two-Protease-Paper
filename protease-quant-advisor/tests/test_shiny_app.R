#!/usr/bin/env Rscript
# test_shiny_app.R — the GUI's server logic, driven headlessly with shiny::testServer on the bundled
# synthetic data: load the four matrices, press "Run advisor", and render every table and plot output.
# Skips (exit 0) if the GUI packages (shiny, DT) are not installed.
#
# Usage: Rscript tests/test_shiny_app.R

if (!all(vapply(c("shiny", "DT"), requireNamespace, logical(1), quietly = TRUE))) {
  cat("  SKIP  shiny/DT not installed (Rscript install.R installs them)\n"); quit(status = 0)
}
suppressPackageStartupMessages({ library(shiny); library(data.table) })
ROOT <- normalizePath(file.path(dirname(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
EX <- file.path(ROOT, "example_data")
CFG0 <- yaml::read_yaml(file.path(ROOT, "config_default.yml"))   # what the GUI shows by default
fi <- function(name) data.frame(name = name, size = file.size(file.path(EX, name)), type = "text/tab-separated-values",
                                datapath = file.path(EX, name), stringsAsFactors = FALSE)

.pass <- 0L; .fail <- 0L
check <- function(cond, msg) {
  if (isTRUE(cond)) { .pass <<- .pass + 1L; cat(sprintf("  PASS  %s\n", msg)) }
  else              { .fail <<- .fail + 1L; cat(sprintf("  FAIL  %s\n", msg)) }
}

app <- shinyAppDir(file.path(ROOT, "app"))
testServer(app, {
  # testServer does not apply the UI's initial values, so set every sidebar input to the value the UI
  # would show (the UI takes them from config_default.yml, i.e. CFG0)
  session$setInputs(
    pr_a = fi("trypsin_report.pr_matrix.tsv"), pr_b = fi("hta_report.pr_matrix.tsv"),
    pg_a = fi("trypsin_report.pg_matrix.tsv"), pg_b = fi("hta_report.pg_matrix.tsv"),
    label_a = CFG0$label_a, label_b = CFG0$label_b, run_include = "", run_exclude = "",
    run_sample_regex = "^[0-9]+_([0-9]+)_[0-9]+$",
    det_min = CFG0$qualify$detection_rate_min, coh_min = CFG0$qualify$coherence_min,
    cv_max = CFG0$qualify$median_cv_max, sel_mode = CFG0$recommend$selection_mode,
    cc_method = CFG0$concordance$method, cc_pct = CFG0$concordance$null_percentile,
    cc_minpaired = CFG0$concordance$min_paired, cc_nperm = 200,
    probe_on = CFG0$probe$enabled, probe_score = CFG0$probe$in_score,
    imp_on = TRUE, imp_model = CFG0$impute$offset_model)
  session$setInputs(run = 1)
  r <- RES()
  check(!inherits(r, "advisor_error"), if (inherits(r, "advisor_error")) paste("run failed:", r$error) else "run completes from the GUI")
  check(nrow(r$rec) > 50L, sprintf("GUI run adjudicated %d shared proteins", nrow(r$rec)))
  check(identical(r$cfg$impute$offset_model, "protein_median") &&
          identical(r$cfg$recommend$selection_mode, "completeness"), "GUI defaults = manuscript method")
  check(sum(r$final_matrix$value_source == "imputed_cross_protease") > 0L, "GUI run with imputation on fills gaps")
  outs <- c("status", "headline", "ingest", "qual_sum", "conc_sum",
            "p_detect", "p_sweep", "p_evi", "p_conc", "p_rec", "p_conseq", "p_imp", "p_qcmet")
  bad <- character()
  for (o in outs) {
    v <- tryCatch(output[[o]], error = function(e) structure(conditionMessage(e), class = "err"))
    if (inherits(v, "err") && !grepl("^$", v)) bad <- c(bad, sprintf("%s (%s)", o, v))
  }
  check(!length(bad), sprintf("%d outputs render%s", length(outs),
                           if (length(bad)) paste0(" — failed: ", paste(bad, collapse = "; ")) else ""))
  rep <- tryCatch(output$dl_report, error = function(e) NULL)
  check(!is.null(rep) && file.exists(rep) && file.size(rep) > 10000, "HTML report downloads")
})

cat(sprintf("\n==== %d passed, %d failed ====\n", .pass, .fail))
quit(status = if (.fail > 0L) 1L else 0L)
