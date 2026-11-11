# seq_coverage.R — FASTA-gated residue coverage for a protein x protease.
#
# Adapted from build_hta_residue_coverage.R (2026-07-31), which established that HTA "sparsity" must
# be judged on RESIDUES covered, not on a peptide count: a protein with one long peptide and one with
# one short peptide are not equally supported. Here the same idea supplies two completeness criteria:
#   seq_coverage_frac         union of all qualifying peptides' residues / protein length
#   median_unit_coverage_frac per-unit (run/sample) residue coverage, median over units
#
# DEGRADES GRACEFULLY: with no FASTA every column is NA, and NA criteria drop out of the weighted
# recommendation score's numerator AND denominator (see recommend_rollup.R). Nothing else changes.
#
# Depends on: data.table.

suppressPackageStartupMessages({ library(data.table) })

# Read a UniProt-style FASTA into data.table(accession, seq, len). Handles `sp|P12345|NAME` and
# bare-accession headers.
read_fasta <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  lines <- readLines(path, warn = FALSE)
  h <- which(startsWith(lines, ">"))
  if (!length(h)) return(NULL)
  starts <- h + 1L
  ends <- c(h[-1] - 1L, length(lines))
  seqs <- vapply(seq_along(h), function(i) {
    if (starts[i] > ends[i]) "" else paste0(lines[starts[i]:ends[i]], collapse = "")
  }, character(1))
  hdr <- sub("^>", "", lines[h])
  acc <- ifelse(grepl("^[a-z]{2}\\|", hdr), sub("^[a-z]{2}\\|([^|]+)\\|.*$", "\\1", hdr),
                sub("[ \t].*$", "", hdr))
  dt <- data.table(accession = acc, seq = toupper(gsub("[^A-Za-z]", "", seqs)))
  dt <- dt[nzchar(accession) & nzchar(seq)]
  dt[, len := nchar(seq)]
  unique(dt, by = "accession")
}

# Sequence for a canonical protein group: first accession in the group that the FASTA knows.
# pg_canonical is a sorted ";"-joined accession set (utils_ids.R::canonicalize_pg).
pg_sequence <- function(pg, fa_index) {
  if (is.null(fa_index)) return(NULL)
  for (a in strsplit(pg, ";", fixed = TRUE)[[1]]) {
    a <- trimws(a)
    s <- fa_index[[a]]
    if (!is.null(s)) return(s)
    a2 <- sub("-[0-9]+$", "", a)            # isoform suffix -> canonical entry
    s <- fa_index[[a2]]
    if (!is.null(s)) return(s)
  }
  NULL
}

# Residue mask (logical of length nchar(seq)) covered by a set of peptide sequences.
.residue_mask <- function(seq, peptides) {
  mask <- logical(nchar(seq))
  for (p in unique(peptides)) {
    if (is.na(p) || !nzchar(p)) next
    st <- gregexpr(p, seq, fixed = TRUE)[[1]]
    if (st[1] == -1L) next
    L <- nchar(p)
    for (s in st) mask[s:(s + L - 1L)] <- TRUE
  }
  mask
}

# Coverage for one protein x protease.
#   peptides_all  : character vector of stripped sequences (qualifying, all units pooled)
#   peptides_byunit: list of character vectors, one per unit (run/sample)
# Returns list(seq_coverage_frac, median_unit_coverage_frac, protein_len, n_peptides_mapped).
protein_coverage <- function(pg, peptides_all, peptides_byunit, fa_index) {
  s <- pg_sequence(pg, fa_index)
  if (is.null(s) || !nzchar(s))
    return(list(seq_coverage_frac = NA_real_, median_unit_coverage_frac = NA_real_,
                protein_len = NA_integer_, n_peptides_mapped = NA_integer_))
  L <- nchar(s)
  mapped <- vapply(unique(peptides_all), function(p) grepl(p, s, fixed = TRUE), logical(1))
  m <- .residue_mask(s, peptides_all)
  per_unit <- if (length(peptides_byunit)) {
    vapply(peptides_byunit, function(pep) sum(.residue_mask(s, pep)) / L, numeric(1))
  } else numeric(0)
  list(seq_coverage_frac = sum(m) / L,
       median_unit_coverage_frac = if (length(per_unit)) stats::median(per_unit) else NA_real_,
       protein_len = L, n_peptides_mapped = sum(mapped))
}

# Index a FASTA table as an environment for O(1) accession lookup.
fasta_index <- function(fa) {
  if (is.null(fa) || !nrow(fa)) return(NULL)
  e <- new.env(hash = TRUE, parent = emptyenv(), size = nrow(fa))
  for (i in seq_len(nrow(fa))) assign(fa$accession[i], fa$seq[i], envir = e)
  structure(list(env = e), class = "fasta_index")
}
`[[.fasta_index` <- function(x, name) {
  if (exists(name, envir = x$env, inherits = FALSE)) get(name, envir = x$env, inherits = FALSE) else NULL
}
