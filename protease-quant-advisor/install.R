#!/usr/bin/env Rscript
# install.R — install the R packages the Protease Quant Advisor needs.
#
#   Rscript install.R          # engine + GUI packages
#   Rscript install.R --engine # engine only (skip the Shiny GUI packages)
#
# Idempotent: only packages that are not already installed are fetched from CRAN.

args    <- commandArgs(trailingOnly = TRUE)
engine  <- c("data.table", "ggplot2", "yaml", "base64enc")
gui     <- c("shiny", "DT")
want    <- if ("--engine" %in% args) engine else c(engine, gui)

# R version check --------------------------------------------------------------------------------
rv <- getRversion()
cat(sprintf("R version: %s\n", rv))
if (rv < "4.5")
  warning(sprintf("This app was developed and tested on R >= 4.5 (you have %s). It may still work, ",
                  "but please upgrade if you hit trouble.", rv))

# Install what is missing ------------------------------------------------------------------------
have    <- rownames(installed.packages())
missing <- setdiff(want, have)
if (length(missing)) {
  cat(sprintf("Installing from CRAN: %s\n", paste(missing, collapse = ", ")))
  install.packages(missing, repos = "https://cloud.r-project.org")
} else {
  cat("All required packages are already installed.\n")
}

# Report -----------------------------------------------------------------------------------------
have <- rownames(installed.packages())
ok   <- want %in% have
cat("\nStatus:\n")
for (p in want)
  cat(sprintf("  [%s] %s%s\n", if (p %in% have) "x" else " ", p,
              if (p %in% gui) "  (GUI only)" else ""))
if (all(ok)) {
  cat("\nReady. Try:  Rscript tests/run_all_tests.R\n")
} else {
  cat(sprintf("\nStill missing: %s. Install these manually and re-run.\n",
              paste(want[!ok], collapse = ", ")))
  quit(status = 1L)
}
