# harmonize.R — generalized front-end for the concordance engine.
# Takes TWO DIA-NN-style pg matrices (channel A and channel B), selects a shared set of
# sample columns, canonicalizes protein groups, and returns aligned per-sample long log2
# tables plus the A_only / B_only / shared disposition. Dataset-agnostic: no Mix20/Mix30
# tokens are hardwired — the sample universe is chosen by pattern/explicit list/auto-intersect.
#
# Depends on: utils_ids.R (canonicalize_pg, leading_accession, is_single_accession),
#             utils_stats.R (to_log2, detect_is_log).

suppressPackageStartupMessages({ library(data.table) })

# Metadata (non-sample) columns commonly emitted by DIA-NN pg/pr matrices.
.PG_META_DEFAULT <- c("Protein.Group","Protein.Ids","Protein.Names","Genes",
                      "First.Protein.Description","N.Sequences","N.Proteotypic.Sequences",
                      "Stripped.Sequence","Modified.Sequence","Precursor.Id","Precursor.Charge",
                      "Proteotypic")

read_pg_matrix <- function(path) {
  fread(path, sep = "\t", na.strings = c("","NA","NaN"), check.names = FALSE)
}

# Which columns are sample (quantity) columns? Everything that is not a known metadata column.
sample_columns <- function(mat, extra_meta = character()) {
  setdiff(names(mat), union(.PG_META_DEFAULT, extra_meta))
}

# Choose the shared sample universe across the two matrices.
#   sample_pattern : optional regex; keep only sample columns whose name matches it (e.g. "30min_Mix").
#   sample_cols    : optional explicit character vector of column names to use.
# Returns the intersection of (matching) sample columns present in BOTH matrices.
select_shared_samples <- function(matA, matB, sample_pattern = NULL, sample_cols = NULL) {
  sa <- sample_columns(matA); sb <- sample_columns(matB)
  shared <- intersect(sa, sb)
  if (!is.null(sample_cols)) shared <- intersect(shared, sample_cols)
  if (!is.null(sample_pattern)) shared <- shared[grepl(sample_pattern, shared)]
  if (!length(shared))
    stop("no shared sample columns after selection — check inputs / sample_pattern / sample_cols")
  shared
}

# Wide pg matrix -> long log2 table keyed by canonical protein group over the given sample cols.
# `id_fn` optionally maps a raw column name to a biological sample id (default: identity).
# Canonicalization-duplicate guard: identical rows deduped; conflicting rows stop the run.
harmonize_to_long <- function(mat, sample_cols, already_log = NA, id_fn = identity) {
  pgcol <- if ("Protein.Group" %in% names(mat)) "Protein.Group" else names(mat)[1]
  keep  <- c(pgcol, sample_cols)
  m <- mat[, ..keep]
  setnames(m, pgcol, "protein_group_original")
  m[, pg_canonical := canonicalize_pg(protein_group_original)]
  m <- m[!is.na(pg_canonical) & nzchar(pg_canonical)]
  dup <- m[, .N, by = pg_canonical][N > 1, pg_canonical]
  if (length(dup)) {
    for (d in dup) {
      sub <- m[pg_canonical == d, ..sample_cols]
      if (nrow(unique(sub)) > 1)
        stop(sprintf("canonicalization created a conflicting duplicate for '%s' — stop & inspect", d))
    }
    m <- m[!duplicated(pg_canonical)]
  }
  long <- melt(m, id.vars = c("protein_group_original","pg_canonical"),
               measure.vars = sample_cols, variable.name = "original_column", value.name = "raw")
  long[, original_column := as.character(original_column)]
  long[, sample_id := vapply(original_column, function(x) as.character(id_fn(x)), character(1))]
  long[, log2q := to_log2(raw, already_log = already_log)]
  long[is.finite(log2q), .(pg_canonical, protein_group_original, sample_id, log2q)]
}

# Per-protein-group annotation (genes, accessions, proteotypic count, single-accession flag).
pg_annotation <- function(mat) {
  pgcol <- if ("Protein.Group" %in% names(mat)) "Protein.Group" else names(mat)[1]
  a <- data.table(protein_group_original = mat[[pgcol]])
  a[, pg_canonical := canonicalize_pg(protein_group_original)]
  a[, genes := if ("Genes" %in% names(mat)) mat[["Genes"]] else NA_character_]
  a[, accessions := if ("Protein.Ids" %in% names(mat)) mat[["Protein.Ids"]] else protein_group_original]
  a[, n_proteotypic := if ("N.Proteotypic.Sequences" %in% names(mat)) as.integer(mat[["N.Proteotypic.Sequences"]]) else NA_integer_]
  a[, single_accession := is_single_accession(pg_canonical)]
  a[, leading := leading_accession(pg_canonical)]
  a[!is.na(pg_canonical) & !duplicated(pg_canonical)]
}

# Top-level: read both matrices, select samples, harmonize, and summarize disposition.
# Returns a list with A_long, B_long, annotations, the shared sample vector, and counts.
harmonize_two <- function(pathA, pathB, sample_pattern = NULL, sample_cols = NULL,
                          already_log = NA, id_fn = identity,
                          labelA = "channel_A", labelB = "channel_B") {
  matA <- read_pg_matrix(pathA); matB <- read_pg_matrix(pathB)
  shared_samples <- select_shared_samples(matA, matB, sample_pattern, sample_cols)
  A_long <- harmonize_to_long(matA, shared_samples, already_log, id_fn)
  B_long <- harmonize_to_long(matB, shared_samples, already_log, id_fn)
  a_det <- unique(A_long$pg_canonical); b_det <- unique(B_long$pg_canonical)
  shared_pg <- intersect(a_det, b_det)
  disposition <- data.table(
    class = c(paste0(labelA, "_only"), "shared", paste0(labelB, "_only")),
    n = c(length(setdiff(a_det, b_det)), length(shared_pg), length(setdiff(b_det, a_det))))
  list(A_long = A_long, B_long = B_long,
       annA = pg_annotation(matA), annB = pg_annotation(matB),
       shared_samples = shared_samples, n_samples = length(shared_samples),
       shared_pg = shared_pg, disposition = disposition,
       labelA = labelA, labelB = labelB)
}
