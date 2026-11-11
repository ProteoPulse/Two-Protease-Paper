#!/usr/bin/env Rscript
# render_reports.R — re-render every condition's report.html, INTERPRETATION_AND_METHODS.md and
# tables from its saved advisor_state.rds, WITHOUT recomputing the analysis.
#
# Stage 2 takes minutes and is deterministic, so once a run's state is on disk there is no reason to
# pay for it again just to change how the results are described. Use this after editing outputs.R,
# methods_doc.R or plots.R; use run_all_conditions.sh only when the ANALYSIS changes.
#
# Optional args: tag=mix30 (repeatable) to re-render a subset.

suppressPackageStartupMessages({ library(data.table) })
ROOT <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(ROOT, "engine", "load_engine.R")); load_advisor_engine(ROOT)
RES <- file.path(ROOT, "results")

args <- commandArgs(TRUE)
want <- sub("^tag=", "", grep("^tag=", args, value = TRUE))
tags <- basename(dirname(Sys.glob(file.path(RES, "*", "advisor_state.rds"))))
if (length(want)) tags <- intersect(tags, want)
if (!length(tags)) stop("no saved runs to re-render under results/")

for (t in tags) {
  cat("── re-rendering", t, "\n")
  d <- file.path(RES, t)
  res <- readRDS(file.path(d, "advisor_state.rds"))
  write_outputs(res, d, tag = t)
  render_advisor_report(res, d, tag = t)
  write_methods_doc(res, d, tag = t)
  cat("   report.html + INTERPRETATION_AND_METHODS.md + tables/ refreshed\n")
}
cat("\ndone:", paste(tags, collapse = ", "), "\n")
