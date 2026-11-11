# feature_matrix.R — generalized (unit x feature) matrix builder.
#
# Generalizes core_matrix() (diagnostics_peptide.R) so the SAME dcast+filter logic serves three
# consumers at different grains without copy-paste:
#   Stage 1 qualification  : (run x precursor)  matrix for one protein x protease
#   Stage 2 roll-up input  : (sample x feature)  matrix fed to Tukey median polish
#   Stage 7 diagnostics    : (sample x peptide)  matrix (the original core_matrix use)
#
# A long table `d` with an id column (rows), a feature column (cols), and a value column is cast to
# a numeric matrix. Optionally keep only "core" features (observed in >= max(min_obs, core_frac*n)
# units, >= min_distinct distinct values, variance > 0) — the same quality gate core_matrix() uses.
#
# Depends on: data.table.

suppressPackageStartupMessages({ library(data.table) })

# Build a (unit x feature) numeric matrix from a long table.
#   d          : data.table with columns named by id_col / feature_col / value_col
#   id_col     : rows of the matrix (e.g. "sample_id" or "run")
#   feature_col: columns of the matrix (e.g. "precursor_id" or "stripped")
#   value_col  : cell values (e.g. "log2q")
#   core_only  : if TRUE, drop features failing the core quality gate
# Returns a matrix with rownames = id values, or NULL if empty / no surviving features.
feature_matrix <- function(d, id_col = "sample_id", feature_col = "stripped",
                           value_col = "log2q", core_only = FALSE,
                           core_frac = 0.50, min_obs = 20, min_distinct = 5) {
  if (is.null(d) || !nrow(d)) return(NULL)
  n_unit <- d[, uniqueN(get(id_col))]
  f <- stats::as.formula(sprintf("%s ~ %s", id_col, feature_col))
  w <- dcast(d, f, value.var = value_col, fun.aggregate = function(v) {
    v <- v[is.finite(v)]; if (!length(v)) NA_real_ else stats::median(v)
  })
  ids <- w[[1]]
  m <- as.matrix(w[, -1, drop = FALSE]); rownames(m) <- as.character(ids)
  if (!ncol(m)) return(NULL)
  if (core_only) {
    keep <- apply(m, 2, function(v) {
      o <- is.finite(v)
      sum(o) >= max(min_obs, core_frac * n_unit) &&
        length(unique(v[o])) >= min_distinct && stats::var(v[o]) > 0
    })
    m <- m[, keep, drop = FALSE]
    if (!ncol(m)) return(NULL)
  }
  m
}

# Median pairwise correlation among the COLUMNS (features) of a matrix, on pairwise-complete
# overlaps. `method` selects Pearson (Stage-1 sibling coherence; amplitude-aware, per the plan) or
# Spearman (rank coherence). Returns NA for <2 columns. This is the method-parameterized
# generalization of median_pairwise_coherence() (diagnostics_peptide.R, Spearman-only).
median_pairwise_cor <- function(m, method = "pearson") {
  if (is.null(m) || ncol(m) < 2) return(NA_real_)
  cc <- suppressWarnings(stats::cor(m, method = method, use = "pairwise.complete.obs"))
  stats::median(cc[upper.tri(cc)], na.rm = TRUE)
}

