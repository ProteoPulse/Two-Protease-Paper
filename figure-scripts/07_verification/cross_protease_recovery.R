# =============================================================================
# Which protease recovers the other's missing values?   2026-09-28
# Hypothesis under test: "Krakatoa is not really recovering trypsin's missed values; trypsin is
# helping Krakatoa." Tested as recovery RATES (not raw counts), in the discovery
# conditions (Mix30, Mix20; 80:20 trypsin:Krakatoa) and the reproducibility cohort
# (50:50, 8 samples x 3 days x 8 replicates = 192 runs).
#
# For each cohort and each recipient protease:
#   gaps          = protein x run cells missing in the recipient, over proteins it quantifies >= once
#   donor_has     = gaps where the other protease measured that protein in that run (ceiling)
#   filled        = gaps filled under the strict gate (>= 5 co-observed runs; recipient seen >= once)
#   in_final      = filled cells that land in the final one-channel-per-protein matrix
# Same fill rule as 05_cross_protease (impute_runs); totals are asserted against it.
# =============================================================================
suppressPackageStartupMessages(library(data.table))
.sd <- function() { a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f)) normalizePath(dirname(f)) else getwd() }
source(file.path(.sd(), "..", "config.R"))   # DATA_DIR, OUT_ROOT, step_dirs()
source(file.path(.sd(), "helpers.R"))
TBL <- step_dirs("07_verification")$tables

# protein x run presence (runs as columns)
run_mat <- function(L) {
  w <- dcast(L[is.finite(Intensity), .(FeatureID, Run, p = TRUE)], FeatureID ~ Run, value.var = "p", fill = FALSE)
  M <- as.matrix(w[, -1]); rownames(M) <- w$FeatureID; M
}
recovery <- function(Tp, Kp, cohort) {
  runs <- union(colnames(Tp), colnames(Kp)); prots <- union(rownames(Tp), rownames(Kp))
  ex <- function(M) { o <- matrix(FALSE, length(prots), length(runs), dimnames = list(prots, runs))
    o[rownames(M), colnames(M)] <- M; o }
  Tp <- ex(Tp); Kp <- ex(Kp)
  fill <- function(tgt, don) {                       # strict gate
    ok <- rowSums(tgt & don) >= MIN_PAIRED & rowSums(tgt) > 0
    !tgt & don & ok
  }
  fT <- fill(Tp, Kp); fK <- fill(Kp, Tp)
  chosen <- ifelse(rowSums(Tp) >= rowSums(Kp), "Trypsin", "Krakatoa")
  one <- function(rec, R, D, F) {
    seen <- rowSums(R) > 0; R <- R[seen, , drop = FALSE]; D <- D[seen, , drop = FALSE]; F <- F[seen, , drop = FALSE]
    gaps <- !R
    data.table(cohort = cohort, recipient = rec, runs = ncol(R), proteins = nrow(R),
               shared_proteins = sum(rowSums(D) > 0),
               measured = sum(R), gaps = sum(gaps), donor_has = sum(gaps & D), filled = sum(F),
               in_final = sum(F[chosen[rownames(F)] == rec, , drop = FALSE]))
  }
  rbind(one("Trypsin", Tp, Kp, fT), one("Krakatoa", Kp, Tp, fK))
}

TR <- load_disc(PG_TRY); KR <- load_disc(PG_KRA)
disc <- rbindlist(lapply(c("Mix30", "Mix20"), function(cd)
  recovery(run_mat(TR[Cond == cd]), run_mat(KR[Cond == cd]), cd)))
ld_rep <- function(p) melt_pg(fread(p, sep = "\t", header = TRUE, showProgress = FALSE))
rep <- recovery(run_mat(ld_rep(RP_TRY)), run_mat(ld_rep(RP_HTA)), "Reproducibility (50:50)")
out <- rbind(disc, rep)
out[, `:=`(pct_gaps_other_measured = round(100 * donor_has / gaps, 1),
           pct_gaps_filled = round(100 * filled / gaps, 1),
           pct_measured_before = round(100 * measured / (proteins * runs), 1),
           pct_after_fill = round(100 * (measured + filled) / (proteins * runs), 1))]

# guard: must reproduce the Sep-11 figure set (impact_*_datapoints.csv)
chk <- out[, .(cohort, recipient, filled)]
exp <- data.table(cohort = rep(c("Mix30","Mix20","Reproducibility (50:50)"), each = 2),
                  recipient = rep(c("Trypsin","Krakatoa"), 3), filled = c(494L, 4541L, 298L, 4152L, 719L, 2598L))
stopifnot(identical(merge(exp, chk, by = c("cohort","recipient"))[, filled.x == filled.y], rep(TRUE, 6L)))
cat("Fill counts reproduce the Sep-11 figure set (494/4541, 298/4152, 719/2598).\n")
cat("Cells in the final matrix: Mix30", sum(out[cohort == "Mix30", in_final]), "| Mix20", sum(out[cohort == "Mix20", in_final]),
    "| Repro", sum(out[cohort == "Reproducibility (50:50)", in_final]), "\n\n")
print(out[, .(cohort, recipient, proteins, measured, gaps, donor_has, pct_gaps_other_measured, filled, pct_gaps_filled,
              in_final, pct_measured_before, pct_after_fill)])
fwrite(out, file.path(TBL, "cross_protease_recovery.csv"))
