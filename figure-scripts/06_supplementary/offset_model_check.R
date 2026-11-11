# =============================================================================
# offset_model_check.R — does the offset model matter for repro Panel C?
#
# The Aug-24 benchmark that chose sample_knn over the per-protein median ran on
# Mix30/Mix20 (both DISCOVERY). The reproducibility cohort was never tested, and
# Panel C used the per-protein median as a TRANSPARENCY choice, not a benchmarked
# one. Hypothesis: in repro, every run of a sample is a near-replicate, so the
# k nearest runs carry almost the same offset as the whole protein — the two
# models should give near-identical counts.
#
# This reruns the repro Panel C imputation under BOTH offset models and compares
# the well-quantified counts (published inter-day CV, strict tier). Panel A is
# asserted first.
#
# Usage: Rscript figure-scripts/06_supplementary/offset_model_check.R
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

.sd <- function() { a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f)) normalizePath(dirname(f)) else getwd() }
source(file.path(.sd(), "..", "config.R"))   # DATA_DIR, OUT_ROOT, step_dirs()
TBL  <- step_dirs("06_supplementary")$tables

PG_TRY <- file.path(REPRO, "Trypsin Library",  "report.pg_matrix.tsv")
PG_HTA <- file.path(REPRO, "Krakatoa Library", "report.pg_matrix.tsv")
ANN <- c("Genes","Protein.Group","Protein.Ids","Protein.Names",
         "First.Protein.Description","N.Sequences","N.Proteotypic.Sequences")
CV_THRESHOLD <- 20; MIN_PAIRED <- 5L; KNN_K <- 10L

# ---- load (identical to panel_c_imputation.R) -------------------------------
load_pg <- function(path, search) {
  dt <- fread(path, sep = "\t", header = TRUE, showProgress = FALSE)
  run_cols <- setdiff(names(dt), ANN)
  L <- melt(dt, id.vars = "Protein.Group", measure.vars = run_cols,
            variable.name = "Sample", value.name = "Intensity", variable.factor = FALSE)
  L[, Intensity := as.numeric(Intensity)]
  L[!is.finite(Intensity) | Intensity <= 0, Intensity := NA_real_]
  L <- L[, .(Intensity = if (all(is.na(Intensity))) NA_real_ else sum(Intensity, na.rm = TRUE)),
         by = .(FeatureID = Protein.Group, Sample)]
  L <- L[!is.na(FeatureID) & FeatureID != ""]
  b <- sub("\\.raw$", "", basename(gsub("\\\\", "/", L$Sample)), ignore.case = TRUE)
  p <- tstrsplit(b, "_", fixed = TRUE); n <- length(p)
  L[, Day := as.integer(p[[n - 2]])]; L[, SampleID := p[[n - 1]]]
  L[, Replicate := as.integer(p[[n]])]; L[, Search := search][, imputed := FALSE][]
}
inter_cv <- function(L) {
  dm <- L[is.finite(Intensity) & Intensity > 0,
          .(DayMean = mean(Intensity), n = .N, n_imp = sum(imputed)), by = .(FeatureID, SampleID, Day)]
  ed <- L[, .(ExpectedDays = uniqueN(Day)), by = SampleID]
  out <- dm[, .(CV = if (.N < 2) NA_real_ else 100 * sd(DayMean) / mean(DayMean),
                MeasuredDays = .N, days_fully_imputed = sum(n_imp == n)), by = .(FeatureID, SampleID)]
  out <- merge(out, ed, by = "SampleID")
  out[, CompleteDays := MeasuredDays == ExpectedDays][]
}
wq <- function(cv, tier = "strict") {
  keep <- switch(tier, strict = cv$days_fully_imputed == 0, ungated = rep(TRUE, nrow(cv)))
  cv[keep & CompleteDays %in% TRUE & MeasuredDays >= 2 & !is.na(CV) & CV <= CV_THRESHOLD]
}

# ---- offset models ----------------------------------------------------------
# per-protein median offset (the transparent incumbent)
offsets_median <- function(target, donor) {
  m <- merge(target[, .(FeatureID, Sample, tI = Intensity)],
             donor[,  .(FeatureID, Sample, dI = Intensity)], by = c("FeatureID","Sample"))
  m <- m[is.finite(tI) & tI > 0 & is.finite(dI) & dI > 0][, d := log2(tI) - log2(dI)]
  o <- m[, .(n_paired = .N), by = FeatureID][n_paired >= MIN_PAIRED]
  list(pairs = m, elig = o$FeatureID)
}
# run-run similarity from the donor channel (Spearman over well-detected proteins)
run_similarity <- function(donor) {
  ns <- uniqueN(donor$Sample)
  det <- donor[is.finite(Intensity), .(det = uniqueN(Sample)/ns), by = FeatureID]
  keep <- det[det >= 0.90, FeatureID]
  if (length(keep) < 10L) keep <- det[order(-det)][seq_len(min(50L, .N)), FeatureID]
  w <- dcast(donor[FeatureID %in% keep & is.finite(Intensity)], FeatureID ~ Sample, value.var = "Intensity")
  M <- log2(as.matrix(w[, -1, drop = FALSE]))
  suppressWarnings(stats::cor(M, use = "pairwise.complete.obs", method = "spearman"))
}

# impute gaps in target from donor, using per-(FeatureID,Sample) offset supplied by `model`
impute_model <- function(target, donor, model = c("median","sample_knn")) {
  model <- match.arg(model)
  om <- offsets_median(target, donor); elig <- om$elig; pairs <- om$pairs
  m <- merge(target, donor[, .(FeatureID, Sample, dI = Intensity)],
             by = c("FeatureID","Sample"), all.x = TRUE)
  gaps <- m[is.na(Intensity) & is.finite(dI) & dI > 0 & FeatureID %in% elig]
  if (model == "median") {
    med <- pairs[, .(offset = median(d)), by = FeatureID]
    gaps <- merge(gaps, med, by = "FeatureID", all.x = TRUE)
  } else {
    S <- run_similarity(donor)
    by_pg <- split(pairs[, .(Sample, off_s = d)], pairs$FeatureID)
    gaps[, offset := {
      vapply(seq_len(.N), function(i) {
        o <- by_pg[[FeatureID[i]]]; if (is.null(o) || !nrow(o)) return(NA_real_)
        s <- Sample[i]
        if (is.null(S) || !(s %in% rownames(S))) return(median(o$off_s))
        sim <- S[s, o$Sample]; ok <- is.finite(sim)
        if (!any(ok)) return(median(o$off_s))
        idx <- order(-sim[ok])[seq_len(min(KNN_K, sum(ok)))]
        median(o$off_s[which(ok)][idx])
      }, numeric(1))
    }]
  }
  gaps <- gaps[is.finite(offset)]
  m[, imputed := FALSE]
  key <- paste(gaps$FeatureID, gaps$Sample)
  mk  <- paste(m$FeatureID, m$Sample)
  hit <- match(mk, key)
  fill <- !is.na(hit)
  m[fill, Intensity := 2^(log2(gaps$dI[hit[fill]]) + gaps$offset[hit[fill]])]
  m[fill, imputed := TRUE]
  m[, dI := NULL][]
}

# ---- run --------------------------------------------------------------------
cat("========================================\n")
cat("OFFSET MODEL CHECK — reproducibility Panel C\n")
cat("========================================\n\n")
try_L <- load_pg(PG_TRY, "Trypsin"); hta_L <- load_pg(PG_HTA, "Krakatoa")
so <- as.character(sort(as.integer(unique(try_L$SampleID))))
cnt <- function(cv, tier = "strict") {
  w <- wq(cv, tier)[, .N, by = SampleID]; v <- w$N[match(so, w$SampleID)]; v[is.na(v)] <- 0L; as.integer(v)
}
PANEL_A_TRY <- c(120L,107L,127L,123L,123L,123L,146L,131L)
PANEL_A_HTA <- c(92L,80L,88L,89L,108L,101L,91L,82L)
bt <- cnt(inter_cv(try_L), "strict"); bh <- cnt(inter_cv(hta_L), "strict")
# baseline strict == baseline measured here (nothing imputed yet)
cat("--- Panel A guardrail ---\n")
cat(sprintf("  Trypsin  : %s  %s\n", paste(bt, collapse=" "), if (identical(bt, PANEL_A_TRY)) "MATCH" else "MISMATCH"))
cat(sprintf("  Krakatoa : %s  %s\n", paste(bh, collapse=" "), if (identical(bh, PANEL_A_HTA)) "MATCH" else "MISMATCH"))
if (!identical(bt, PANEL_A_TRY) || !identical(bh, PANEL_A_HTA)) stop("Panel A did not reproduce.")
cat("  PASSED.\n\n")

res <- rbindlist(lapply(c("median","sample_knn"), function(md) {
  tr_i <- impute_model(try_L, hta_L, md)
  hta_i <- impute_model(hta_L, try_L, md)
  data.table(model = md, SampleID = so,
             Trypsin_base = bt, Trypsin_imp = cnt(inter_cv(tr_i), "strict"),
             Krakatoa_base = bh, Krakatoa_imp = cnt(inter_cv(hta_i), "strict"),
             tr_filled = tr_i[imputed == TRUE, .N], hta_filled = hta_i[imputed == TRUE, .N])
}))
fwrite(res, file.path(TBL, "offset_model_repro_counts.csv"))

cat("--- strict-tier well-quantified totals (sum over 8 samples) ---\n")
tot <- res[, .(Trypsin_base = sum(Trypsin_base), Trypsin_imp = sum(Trypsin_imp),
               Krakatoa_base = sum(Krakatoa_base), Krakatoa_imp = sum(Krakatoa_imp),
               tr_filled = sum(tr_filled)/8, hta_filled = sum(hta_filled)/8), by = model]
print(tot)

cat("\n--- per-sample Krakatoa (the channel that gains most) ---\n")
print(dcast(res, SampleID ~ model, value.var = "Krakatoa_imp"))
cat("\n--- per-sample Trypsin ---\n")
print(dcast(res, SampleID ~ model, value.var = "Trypsin_imp"))

d <- dcast(res, SampleID ~ model, value.var = c("Trypsin_imp","Krakatoa_imp"))
cat(sprintf("\nCounts identical between models?  Trypsin: %s   Krakatoa: %s\n",
    identical(d$Trypsin_imp_median, d$Trypsin_imp_sample_knn),
    identical(d$Krakatoa_imp_median, d$Krakatoa_imp_sample_knn)))
cat(sprintf("Max per-sample abs difference: Trypsin %d, Krakatoa %d\n",
    max(abs(d$Trypsin_imp_median - d$Trypsin_imp_sample_knn)),
    max(abs(d$Krakatoa_imp_median - d$Krakatoa_imp_sample_knn))))
cat("\n========================================\nDONE.\n========================================\n")
