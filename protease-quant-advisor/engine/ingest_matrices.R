# ingest_matrices.R — read a DIA-NN **pivoted** report.pr_matrix.tsv and report.pg_matrix.tsv.
#
# WHY THIS EXISTS. The Aug-19/Aug-20 framework read the long `report.parquet` because it keeps
# per-run Q.Value and lets us apply our own qualification from scratch. This tool deliberately reads
# the PIVOTED pr_matrix instead, because that is the file people actually have and can hand to an app.
#
# WHAT THAT COSTS, STATED PLAINLY. The pr_matrix has already been q-value filtered by DIA-NN
# (Precursor / Protein / PG / Lib Q <= 0.01) and the per-run q-values are gone. So:
#   * we INHERIT DIA-NN's 1% filtering as a floor and cannot loosen it;
#   * the `qvalue_max` qualification rule becomes vacuous and is reported as "not applicable";
#   * the remaining three rules — detection rate, replicate CV, sibling coherence — are computed from
#     intensities and are completely unaffected.
# The tool records this in the run header so no reader has to guess which regime produced a number.
#
# Column layout of a pr_matrix: a fixed block of annotation columns, then ONE COLUMN PER RUN whose
# header is the raw file path (e.g. "A:\...\251117_DIA_Plasma_Mix_BRH000000_20min.raw"). Run columns
# are identified positionally as "everything that is not an annotation column", so a search that adds
# an annotation column does not silently become a run.
#
# Emits the SAME long contract as the parquet loader in the earlier framework, so every vendored
# module downstream works unchanged:
#   pg_canonical, protease, precursor_id, stripped_sequence, modified_sequence, charge,
#   sample_id, run, log2q, q_value (NA here), pg_q_value (NA here), proteotypic
#
# Depends on: data.table, utils_ids.R (canonicalize_pg), utils_stats.R (to_log2).

suppressPackageStartupMessages({ library(data.table) })

# Annotation columns DIA-NN writes before the run columns. Anything else is a run.
.PR_ANNOT <- c("Protein.Group","Protein.Ids","Protein.Names","Genes","First.Protein.Description",
               "Proteotypic","Stripped.Sequence","Modified.Sequence","Precursor.Charge",
               "Precursor.Id","Precursor.Lib.Index","Lib.Q.Value","Lib.PG.Q.Value")
.PG_ANNOT <- c("Protein.Group","Protein.Ids","Protein.Names","Genes","First.Protein.Description",
               "N.Sequences","N.Proteotypic.Sequences","Molecular.Weight")

# Run column headers are full raw-file paths. Reduce to a bare run name: basename, drop extension.
normalize_run <- function(x) {
  b <- sub("^.*[\\\\/]", "", as.character(x))     # strip Windows or POSIX directories
  sub("\\.(raw|d|mzML|mzml|wiff|dia)$", "", b, ignore.case = TRUE)
}

# Which columns of a matrix are runs, after keep/drop regexes on the NORMALIZED run name.
select_run_columns <- function(cols, annot, run_include = NULL, run_exclude = NULL) {
  runs <- setdiff(cols, annot)
  if (!length(runs)) stop("no run columns found — is this really a DIA-NN matrix?")
  nm <- normalize_run(runs)
  keep <- rep(TRUE, length(runs))
  if (!is.null(run_include) && nzchar(run_include)) keep <- keep & grepl(run_include, nm)
  if (!is.null(run_exclude) && nzchar(run_exclude)) keep <- keep & !grepl(run_exclude, nm)
  data.table(col = runs[keep], run = nm[keep])
}

# ---- precursor matrix -> long ------------------------------------------------------------------
#   path        : report.pr_matrix.tsv
#   protease    : channel label; set by WHICH SEARCH the file came from, never inferred
#   run_regex   : optional regex with one capture group; first capture becomes sample_id
#                 (reproducibility runs encode day/sample/replicate; discovery is one run per sample)
#   min_runs_obs: drop precursors seen in fewer than this many of the selected runs (0 = keep all)
load_pr_matrix <- function(path, protease, run_include = NULL, run_exclude = NULL,
                           run_regex = NULL, min_runs_obs = 1L, already_log = NA) {
  hdr <- names(fread(path, nrows = 0L, sep = "\t", check.names = FALSE))
  sel <- select_run_columns(hdr, .PR_ANNOT, run_include, run_exclude)
  if (!nrow(sel)) stop("no run columns survived run_include/run_exclude for ", protease)
  annot_present <- intersect(.PR_ANNOT, hdr)
  d <- fread(path, sep = "\t", select = c(annot_present, sel$col), check.names = FALSE,
             na.strings = c("", "NA", "NaN"))

  id_cols <- intersect(c("Protein.Group","Genes","Stripped.Sequence","Modified.Sequence",
                         "Precursor.Charge","Precursor.Id","Proteotypic"), names(d))
  L <- melt(d, id.vars = id_cols, measure.vars = sel$col,
            variable.name = "run_col", value.name = "intensity", variable.factor = FALSE)
  L <- L[is.finite(intensity) & intensity > 0]
  L[sel, run := i.run, on = .(run_col = col)]

  out <- data.table(
    pg_canonical      = canonicalize_pg(L[["Protein.Group"]]),
    protease          = protease,
    precursor_id      = as.character(L[["Precursor.Id"]]),
    stripped_sequence = as.character(L[["Stripped.Sequence"]]),
    modified_sequence = if ("Modified.Sequence" %in% names(L)) as.character(L[["Modified.Sequence"]]) else NA_character_,
    charge            = if ("Precursor.Charge" %in% names(L)) as.integer(L[["Precursor.Charge"]]) else NA_integer_,
    run               = L$run,
    log2q             = to_log2(L$intensity, already_log = already_log),
    q_value           = NA_real_,     # not in a pr_matrix — see the header note
    pg_q_value        = NA_real_,
    proteotypic       = if ("Proteotypic" %in% names(L)) as.integer(L[["Proteotypic"]]) else NA_integer_
  )
  out[, sample_id := derive_sample_id(run, run_regex)]
  out <- out[is.finite(log2q) & !is.na(pg_canonical) & nzchar(pg_canonical)]
  if (min_runs_obs > 0L) {
    keep <- out[, .N, by = precursor_id][N >= min_runs_obs, precursor_id]
    out <- out[precursor_id %in% keep]
  }
  setcolorder(out, c("pg_canonical","protease","precursor_id","stripped_sequence","modified_sequence",
                     "charge","sample_id","run","log2q","q_value","pg_q_value","proteotypic"))
  setattr(out, "source_file", path); setattr(out, "n_runs_selected", nrow(sel))
  out[]
}

# Same regex-capture rule as the parquet loader, kept identical on purpose.
derive_sample_id <- function(run, run_regex = NULL) {
  if (is.null(run_regex) || !nzchar(run_regex)) return(as.character(run))
  m <- regmatches(run, regexec(run_regex, run))
  vapply(seq_along(run), function(i) {
    g <- m[[i]]; if (length(g) >= 2 && nzchar(g[2])) g[2] else as.character(run[i])
  }, character(1))
}

# ---- protein-group matrix -> long ----------------------------------------------------------------
# This is the QUANT the tool adjudicates. We never recompute it; we only decide which channel's copy
# to trust. Returned at sample grain (replicate runs collapsed by median, matching the evidence grain).
load_pg_matrix <- function(path, protease, run_include = NULL, run_exclude = NULL,
                           run_regex = NULL, already_log = NA) {
  hdr <- names(fread(path, nrows = 0L, sep = "\t", check.names = FALSE))
  sel <- select_run_columns(hdr, .PG_ANNOT, run_include, run_exclude)
  if (!nrow(sel)) stop("no run columns survived run_include/run_exclude for ", protease)
  annot_present <- intersect(.PG_ANNOT, hdr)
  d <- fread(path, sep = "\t", select = c(annot_present, sel$col), check.names = FALSE,
             na.strings = c("", "NA", "NaN"))
  id_cols <- intersect(c("Protein.Group","Genes"), names(d))
  L <- melt(d, id.vars = id_cols, measure.vars = sel$col,
            variable.name = "run_col", value.name = "intensity", variable.factor = FALSE)
  L <- L[is.finite(intensity) & intensity > 0]
  L[sel, run := i.run, on = .(run_col = col)]
  L[, sample_id := derive_sample_id(run, run_regex)]
  out <- L[, .(pg_canonical = canonicalize_pg(Protein.Group),
               protein_group_original = as.character(Protein.Group),
               gene = if ("Genes" %in% names(L)) as.character(Genes) else NA_character_,
               protease = protease, run, sample_id,
               log2q = to_log2(intensity, already_log = already_log))]
  out <- out[is.finite(log2q) & !is.na(pg_canonical) & nzchar(pg_canonical)]
  # RUN grain is returned. Replicate structure is what makes technical CV and replicate ICC
  # measurable, and those must be measured on the value we actually recommend — the vendor PG number
  # — not on an internal roll-up. The caller collapses to sample grain with pg_to_sample_grain().
  setattr(out, "source_file", path); setattr(out, "n_runs_selected", nrow(sel))
  out[]
}

# Collapse a run-grain PG table to the sample grain the evidence and recommendation work at.
pg_to_sample_grain <- function(pg_run) {
  pg_run[, .(log2q = stats::median(log2q), n_runs = .N),
         by = .(pg_canonical, protein_group_original, gene, protease, sample_id)]
}

# Technical CV (linear scale) of the RECOMMENDED quantity across replicate runs within a sample,
# median over samples. NA without replicates.
pg_technical_cv <- function(pg_run) {
  per <- pg_run[, {
    lin <- 2^log2q
    .(cv = if (.N >= 2L && mean(lin) > 0) stats::sd(lin) / mean(lin) else NA_real_)
  }, by = .(pg_canonical, protease, sample_id)]
  per[, .(technical_cv_measured = stats::median(cv, na.rm = TRUE),
          n_samples_with_reps = sum(is.finite(cv))), by = .(pg_canonical, protease)]
}

# Sanity report on a loaded pair — printed into the run log and the HTML header so the reader can see
# what was actually read before any statistic is quoted.
ingest_summary <- function(prec, pg, labelA, labelB) {
  f <- function(d, lab, what) data.table(
    channel = lab, level = what,
    n_rows = d[protease == lab, .N],
    n_proteins = d[protease == lab, uniqueN(pg_canonical)],
    n_features = if ("precursor_id" %in% names(d)) d[protease == lab, uniqueN(precursor_id)] else NA_integer_,
    n_runs = if ("run" %in% names(d)) d[protease == lab, uniqueN(run)] else NA_integer_,
    n_samples = d[protease == lab, uniqueN(sample_id)])
  rbindlist(list(f(prec, labelA, "precursor"), f(prec, labelB, "precursor"),
                 f(pg, labelA, "protein"),    f(pg, labelB, "protein")), fill = TRUE)
}
