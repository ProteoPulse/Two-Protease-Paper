#!/usr/bin/env Rscript
# run_advisor.R — CLI over the same engine the app uses.
#
# Usage:
#   Rscript run_advisor.R pr_a=<trypsin pr_matrix> pr_b=<hta pr_matrix> \
#     pg_a=<trypsin pg_matrix> pg_b=<hta pg_matrix> \
#     [run_include=Plasma_Mix] [run_exclude=20min] [run_sample_regex=...] \
#     [selection_mode=completeness] [impute=true] [offset_model=protein_median] [probe_in_score=false] \
#     [run_day_regex=...] [cv_threshold=20] [completeness_min=0.5] [min_paired=15] [n_perm=1000] \
#     [fasta=<uniprot.fasta>] [config=my_overrides.yml] out=results/<tag> tag=<tag>
#
# Only the keys above are accepted on the command line. Anything else in config_default.yml (QC and
# imputation-gate thresholds, labels, seed, ...) is changed with config=<file>: a YAML file holding just
# the keys to change, layered over config_default.yml.
#
# Example on the bundled synthetic data:
#   Rscript run_advisor.R pr_a=example_data/trypsin_report.pr_matrix.tsv \
#     pr_b=example_data/hta_report.pr_matrix.tsv pg_a=example_data/trypsin_report.pg_matrix.tsv \
#     pg_b=example_data/hta_report.pg_matrix.tsv run_sample_regex='^[0-9]+_([0-9]+)_[0-9]+$' \
#     run_day_regex='^([0-9]+)_.*$' impute=true tag=example
suppressPackageStartupMessages({ library(data.table); library(yaml) })
args <- commandArgs(trailingOnly = TRUE)
kv <- list(); for (a in args) { p <- regmatches(a, regexpr("=", a), invert = TRUE)[[1]]
  if (length(p) == 2) kv[[trimws(p[1])]] <- trimws(p[2]) }
ROOT <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(ROOT, "engine", "load_engine.R")); load_advisor_engine(ROOT)

# config= is layered OVER the defaults, so a partial file only needs the keys it changes.
cfg <- yaml::read_yaml(file.path(ROOT, "config_default.yml"))
if (!is.null(kv$config)) cfg <- utils::modifyList(cfg, yaml::read_yaml(kv$config))
for (k in c("run_include","run_exclude","run_sample_regex")) if (!is.null(kv[[k]])) cfg$ingest[[k]] <- kv[[k]]
if (!is.null(kv$selection_mode)) cfg$recommend$selection_mode <- kv$selection_mode
if (!is.null(kv$impute))         cfg$impute$enabled <- as.logical(kv$impute)
if (!is.null(kv$offset_model))   cfg$impute$offset_model <- kv$offset_model
if (!is.null(kv$probe_in_score)) cfg$probe$in_score <- as.logical(kv$probe_in_score)
if (!is.null(kv$fasta))          cfg$fasta <- kv$fasta
if (!is.null(kv$min_paired))     cfg$concordance$min_paired <- as.integer(kv$min_paired)
if (!is.null(kv$n_perm))         cfg$concordance$n_perm <- as.integer(kv$n_perm)
if (!is.null(kv$run_day_regex))  cfg$reproducibility$run_day_regex <- kv$run_day_regex
if (!is.null(kv$cv_threshold))   cfg$reproducibility$cv_threshold <- as.numeric(kv$cv_threshold)
if (!is.null(kv$completeness_min)) cfg$reproducibility$completeness_min <- as.numeric(kv$completeness_min)
for (r in c("pr_a","pr_b","pg_a","pg_b")) if (is.null(kv[[r]])) stop(r, " is required")

tag <- kv$tag %||% "run"
out_dir <- kv$out %||% file.path(ROOT, "results", tag)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
logf <- file.path(out_dir, "run_log.txt")
options(advisor.progress = function(msg, frac) {
  m <- sprintf("[%3.0f%%] %s", 100*frac, msg); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) })

# Record exactly what was fed in, so the methods document never has to be told separately.
inp <- data.table(
  role = c("pr_matrix_a","pr_matrix_b","pg_matrix_a","pg_matrix_b"),
  channel = c(cfg$label_a, cfg$label_b, cfg$label_a, cfg$label_b),
  level = c("precursor","precursor","protein","protein"),
  path = c(kv$pr_a, kv$pr_b, kv$pg_a, kv$pg_b))
inp[, `:=`(exists = file.exists(path),
           size_mb = round(file.info(path)$size / 1024^2, 2),
           mtime = as.character(file.info(path)$mtime))]
inp[, `:=`(run_include = cfg$ingest$run_include %||% "",
           run_exclude = cfg$ingest$run_exclude %||% "",
           run_sample_regex = cfg$ingest$run_sample_regex %||% "")]
fwrite(inp, file.path(out_dir, "run_inputs.tsv"), sep = "\t")

t0 <- Sys.time()
res <- run_advisor(kv$pr_a, kv$pr_b, kv$pg_a, kv$pg_b, cfg)
write_outputs(res, out_dir, tag)
render_advisor_report(res, out_dir, tag)
write_methods_doc(res, out_dir, tag)
cat(sprintf("\nDONE in %.1f min -> %s\n", as.numeric(difftime(Sys.time(), t0, units="mins")), out_dir))
