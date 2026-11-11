# config.R — the one place the figure scripts learn where the data are and where outputs go.
#
# Every script in figure-scripts/0*/ sources this file first. Set two environment variables (or accept
# the defaults):
#
#   TWO_PROTEASE_DATA   folder holding the DIA-NN search results from MassIVE MSV000103457, laid out as
#                       described in DATA_LAYOUT.md.            default: <repo>/data
#   TWO_PROTEASE_OUT    where figures and tables are written.   default: <repo>/output
#
# Outputs land in <TWO_PROTEASE_OUT>/<script folder>/{figures,tables}, e.g.
# output/05_cross_protease/figures/. Scripts that read another script's tables (see the run order in
# README.md) read them from the same place.

.cfg_dir <- function() {
  # directory of this file, whether it was reached via Rscript --file= or source()
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(dirname(normalizePath(of)))
  }
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) return(normalizePath(file.path(dirname(f[1]), "..")))
  normalizePath("figure-scripts", mustWork = FALSE)
}

FIGSCRIPTS_DIR <- .cfg_dir()
REPO_DIR       <- normalizePath(file.path(FIGSCRIPTS_DIR, ".."), mustWork = FALSE)
ADVISOR_DIR    <- file.path(REPO_DIR, "protease-quant-advisor")

DATA_DIR <- Sys.getenv("TWO_PROTEASE_DATA", file.path(REPO_DIR, "data"))
OUT_ROOT <- Sys.getenv("TWO_PROTEASE_OUT",  file.path(REPO_DIR, "output"))

# The searches used by the manuscript (see DATA_LAYOUT.md)
SEARCH       <- DATA_DIR
TRYPSIN_DISC <- file.path(DATA_DIR, "Trypsin_Library_Search")
KRAKATOA_DISC <- file.path(DATA_DIR, "Krakatoa_Library_Search")
REPRO        <- file.path(DATA_DIR, "Reproducibility_Search")

# Output folders for the calling script: <OUT_ROOT>/<step>/{figures,tables}
step_dirs <- function(step) {
  d <- list(root = file.path(OUT_ROOT, step))
  d$figures <- file.path(d$root, "figures"); d$tables <- file.path(d$root, "tables")
  for (x in d) dir.create(x, recursive = TRUE, showWarnings = FALSE)
  d
}

if (!dir.exists(DATA_DIR))
  warning("TWO_PROTEASE_DATA does not exist: ", DATA_DIR,
          "\n  Download the DIA-NN search results (MassIVE MSV000103457) and point TWO_PROTEASE_DATA at them.",
          call. = FALSE)
