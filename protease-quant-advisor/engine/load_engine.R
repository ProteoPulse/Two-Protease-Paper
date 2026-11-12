# load_engine.R — single entry point that sources the engine in dependency order.
# Both app/app.R and run_advisor.R use this, so neither can drift from the other.
load_advisor_engine <- function(root) {
  for (f in c("utils_ids.R","utils_stats.R","utils_compare.R","harmonize.R","feature_matrix.R",
              "threshold.R","classify.R","amplitude.R","reliability.R","diagnostics_peptide.R",
              "detection.R","qualify_precursors.R","rollup_tmp.R","icc.R","seq_coverage.R",
              "precursor_evidence.R","channel_independence.R","recommend_rollup.R",
              "validate_recommendation.R","report_stage.R"))
    source(file.path(root, "vendor", f))
  for (f in c("utils.R","ingest_matrices.R","impute_pg.R","reproducibility.R","pipeline.R","outputs.R","plots.R","methods_doc.R"))
    source(file.path(root, "engine", f))
  invisible(TRUE)
}
