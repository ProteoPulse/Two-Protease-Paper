# =============================================================================
# panel_c_imputation.R — Panel C: does cross-protease imputation increase the
# number of well-quantified proteins?
#
# Mirrors Panel A exactly (the definitions in 04_reproducibility/Reproducibility_Data_Processing.R):
#   inter-day CV = linear day-mean per day, then 100*sd/mean across days
#   well quantified = CompleteDays (3/3 measured) AND inter-day CV <= 20%
# Panel A numbers are reproduced exactly before any imputation is applied.
#
# Imputation = cross-protease gap fill ONLY. Missing cells in the target
# protease are filled from the other protease rescaled by a per-protein median
# log2 offset. Measured values are never altered, and the two proteases are
# never averaged (no channel merge).
#
# Usage: Rscript figure-scripts/05_cross_protease/panel_c_imputation.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

.sd <- function() { a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f)) normalizePath(dirname(f)) else getwd() }
source(file.path(.sd(), "..", "config.R"))   # DATA_DIR, OUT_ROOT, step_dirs()
D <- step_dirs("05_cross_protease"); FIG <- D$figures; TBL <- D$tables

PG_TRY <- file.path(REPRO, "Trypsin Library",  "report.pg_matrix.tsv")
PG_HTA <- file.path(REPRO, "Krakatoa Library", "report.pg_matrix.tsv")

CV_THRESHOLD <- 20
MIN_PAIRED   <- 5L        # runs needed to estimate a protein's cross-protease offset
PUB_DPI      <- 600

TRY_COL  <- "#3C7DA6"; TRY_LIGHT  <- "#AECBDE"
HTA_COL  <- "#E1812C"; HTA_LIGHT  <- "#F4CDA4"

# ---------------------------------------------------------------- load ------
# Matches read_protein_matrix() in Reproducibility_Data_Processing.R: Protein.Group ids, zero -> NA,
# duplicate ids summed, empty ids dropped, Day_SampleID_Replicate from filename.
load_pg <- function(path, search) {
  dt <- fread(path, sep = "\t", header = TRUE, showProgress = FALSE)
  ann <- c("Genes", "Protein.Group", "Protein.Ids", "Protein.Names",
           "First.Protein.Description", "N.Sequences", "N.Proteotypic.Sequences")
  run_cols <- setdiff(names(dt), ann)
  L <- melt(dt, id.vars = "Protein.Group", measure.vars = run_cols,
            variable.name = "Sample", value.name = "Intensity", variable.factor = FALSE)
  L[, Intensity := as.numeric(Intensity)]
  L[!is.finite(Intensity) | Intensity <= 0, Intensity := NA_real_]
  L <- L[, .(Intensity = if (all(is.na(Intensity))) NA_real_ else sum(Intensity, na.rm = TRUE)),
         by = .(FeatureID = Protein.Group, Sample)]
  L <- L[!is.na(FeatureID) & FeatureID != ""]
  b <- sub("\\.raw$", "", basename(gsub("\\\\", "/", L$Sample)), ignore.case = TRUE)
  p <- tstrsplit(b, "_", fixed = TRUE); n <- length(p)
  L[, Day := as.integer(p[[n - 2]])]
  L[, SampleID := p[[n - 1]]]
  L[, Replicate := as.integer(p[[n]])]
  L[, Search := search][, imputed := FALSE]
  L[]
}

gene_map <- function(path) {
  dt <- fread(path, sep = "\t", header = TRUE, showProgress = FALSE,
              select = c("Protein.Group", "Genes"))
  unique(dt[, .(FeatureID = Protein.Group, Gene = sub(";.*", "", Genes))])[FeatureID != ""]
}

# ------------------------------------------------------------- inter-day ----
# calc_inter_day_cv() + identify_well_quantified_interday() from Reproducibility_Data_Processing.R, plus
# bookkeeping on how much of each day-mean came from imputed values.
inter_cv <- function(L) {
  dm <- L[is.finite(Intensity) & Intensity > 0,
          .(DayMean = mean(Intensity), n = .N, n_imp = sum(imputed)),
          by = .(FeatureID, SampleID, Day)]
  ed <- L[, .(ExpectedDays = uniqueN(Day)), by = SampleID]
  er <- L[, .(ExpectedRuns = uniqueN(Sample)), by = SampleID]
  out <- dm[, .(CV = if (.N < 2) NA_real_ else 100 * sd(DayMean) / mean(DayMean),
                MeasuredDays = .N, MeasuredRuns = sum(n),
                days_fully_imputed = sum(n_imp == n),
                frac_imputed = sum(n_imp) / sum(n)),
            by = .(FeatureID, SampleID)]
  out <- merge(merge(out, ed, by = "SampleID"), er, by = "SampleID")
  out[, CompleteDays := MeasuredDays == ExpectedDays]
  out[, CompleteRuns := MeasuredRuns == ExpectedRuns][]
}

# Tier controls how much of a protein's evidence may be imputed:
#   strict   — every day-mean keeps >=1 real measurement from this protease
#   partial  — at least one day-mean is real (whole days may be imputed)
#   ungated  — no restriction (upper bound; includes proteins with no real data)
# `completeness` chooses what "100% complete" means:
#   days — >=1 detection on each of the 3 preparation days (what Reproducibility_Data_Processing.R
#          actually requires, and what reproduces Panel A)
#   runs — every one of the 24 injections detected (the literal reading of the
#          "100% completeness" label on Panel A)
well_quantified <- function(cv, tier = "strict", completeness = "days") {
  gate <- if (completeness == "runs") cv$CompleteRuns %in% TRUE else cv$CompleteDays %in% TRUE
  keep <- switch(tier,
    strict  = cv$days_fully_imputed == 0,
    onefill = cv$days_fully_imputed <= 1,
    partial = cv$days_fully_imputed <= 2,
    ungated = rep(TRUE, nrow(cv)),
    stop("unknown tier"))
  cv[keep & gate & MeasuredDays >= 2 & !is.na(CV) & CV <= CV_THRESHOLD]
}

# ------------------------------------------------------------- imputation ---
make_offsets <- function(target, donor) {
  m <- merge(target[, .(FeatureID, Sample, tI = Intensity)],
             donor[,  .(FeatureID, Sample, dI = Intensity)], by = c("FeatureID", "Sample"))
  m <- m[is.finite(tI) & tI > 0 & is.finite(dI) & dI > 0]
  m[, d := log2(tI) - log2(dI)]
  o <- m[, .(offset = median(d), offset_mad = mad(d), n_paired = .N), by = FeatureID]
  o[n_paired >= MIN_PAIRED]
}

impute <- function(target, donor, off) {
  m <- merge(target, donor[, .(FeatureID, Sample, dI = Intensity)],
             by = c("FeatureID", "Sample"), all.x = TRUE)
  m <- merge(m, off[, .(FeatureID, offset)], by = "FeatureID", all.x = TRUE)
  m[, imputed := is.na(Intensity) & is.finite(dI) & dI > 0 & !is.na(offset)]
  m[imputed == TRUE, Intensity := 2^(log2(dI) + offset)]
  m[, c("dI", "offset") := NULL]
  m[]
}

# ------------------------------------------------------------------ run -----
cat("========================================\n")
cat("PANEL C — cross-protease imputation\n")
cat("========================================\n\nLoading...\n")
try_L <- load_pg(PG_TRY, "Trypsin")
hta_L <- load_pg(PG_HTA, "Krakatoa")
genes <- unique(rbind(gene_map(PG_TRY), gene_map(PG_HTA)))

cat(sprintf("  Trypsin %d proteins | Krakatoa %d proteins | shared %d | runs %d\n",
            uniqueN(try_L$FeatureID), uniqueN(hta_L$FeatureID),
            length(intersect(unique(try_L$FeatureID), unique(hta_L$FeatureID))),
            uniqueN(try_L$Sample)))

cv_try_b <- inter_cv(try_L); cv_hta_b <- inter_cv(hta_L)

# --- gate 1: Panel A must reproduce exactly ----------------------------------
PANEL_A_TRY <- c(120L, 107L, 127L, 123L, 123L, 123L, 146L, 131L)
PANEL_A_HTA <- c(92L, 80L, 88L, 89L, 108L, 101L, 91L, 82L)
sample_order <- as.character(sort(as.integer(unique(try_L$SampleID))))
cnt <- function(cv, tier = "strict", completeness = "days") {
  w <- well_quantified(cv, tier, completeness)[, .N, by = SampleID]
  v <- w$N[match(sample_order, w$SampleID)]; v[is.na(v)] <- 0L; as.integer(v)
}
base_try <- cnt(cv_try_b); base_hta <- cnt(cv_hta_b)
cat("\n--- Panel A reproduction check ---\n")
cat(sprintf("  Trypsin  : %s  %s\n", paste(base_try, collapse = " "),
            if (identical(base_try, PANEL_A_TRY)) "MATCH" else "MISMATCH"))
cat(sprintf("  Krakatoa : %s  %s\n", paste(base_hta, collapse = " "),
            if (identical(base_hta, PANEL_A_HTA)) "MATCH" else "MISMATCH"))
if (!identical(base_try, PANEL_A_TRY) || !identical(base_hta, PANEL_A_HTA))
  warning("Panel A did not reproduce — Panel C numbers are not comparable to the published panel.")

# --- impute ------------------------------------------------------------------
off_try <- make_offsets(try_L, hta_L)   # trypsin  filled from krakatoa
off_hta <- make_offsets(hta_L, try_L)   # krakatoa filled from trypsin
cat(sprintf("\nOffsets estimated (>=%d paired runs): %d proteins each direction; median MAD %.2f log2\n",
            MIN_PAIRED, nrow(off_try), median(c(off_try$offset_mad, off_hta$offset_mad), na.rm = TRUE)))

try_I <- impute(try_L, hta_L, off_try)
hta_I <- impute(hta_L, try_L, off_hta)
cat(sprintf("Cells filled: Trypsin %d / %d missing | Krakatoa %d / %d missing\n",
            try_I[imputed == TRUE, .N], try_L[is.na(Intensity), .N],
            hta_I[imputed == TRUE, .N], hta_L[is.na(Intensity), .N]))

cv_try_i <- inter_cv(try_I); cv_hta_i <- inter_cv(hta_I)

# --- tier table --------------------------------------------------------------
tiers <- rbindlist(lapply(c("days", "runs"), function(cm)
  rbindlist(lapply(c("Trypsin", "Krakatoa"), function(pr) {
    cb <- if (pr == "Trypsin") cv_try_b else cv_hta_b
    ci <- if (pr == "Trypsin") cv_try_i else cv_hta_i
    data.table(Protease = pr, Completeness = cm, SampleID = sample_order,
               Baseline = cnt(cb, "ungated", cm),
               Imputed_strict  = cnt(ci, "strict",  cm),
               Imputed_onefill = cnt(ci, "onefill", cm),
               Imputed_partial = cnt(ci, "partial", cm),
               Imputed_ungated = cnt(ci, "ungated", cm))
  }))))
tiers[, `:=`(gain_strict  = Imputed_strict  - Baseline,
             gain_onefill = Imputed_onefill - Baseline,
             gain_partial = Imputed_partial - Baseline,
             gain_ungated = Imputed_ungated - Baseline)]
fwrite(tiers, file.path(TBL, "panel_c_tier_counts.csv"))

for (cm in c("days", "runs")) {
cat(sprintf("\n=== Well-quantified per sample — completeness = %s (%s) ===\n", cm,
            if (cm == "days") "all 3 days detected; the published rule" else "all 24 injections detected"))
for (pr in c("Krakatoa", "Trypsin")) {
  d <- tiers[Protease == pr & Completeness == cm]
  cat(sprintf("\n%s\n", pr))
  cat(sprintf("  baseline        : %s\n", paste(sprintf("%3d", d$Baseline), collapse = " ")))
  cat(sprintf("  imputed strict  : %s   mean %+.1f\n",
              paste(sprintf("%3d", d$Imputed_strict), collapse = " "), mean(d$gain_strict)))
  cat(sprintf("  imputed <=1 day : %s   mean %+.1f\n",
              paste(sprintf("%3d", d$Imputed_onefill), collapse = " "), mean(d$gain_onefill)))
  cat(sprintf("  imputed partial : %s   mean %+.1f\n",
              paste(sprintf("%3d", d$Imputed_partial), collapse = " "), mean(d$gain_partial)))
  cat(sprintf("  imputed ungated : %s   mean %+.1f\n",
              paste(sprintf("%3d", d$Imputed_ungated), collapse = " "), mean(d$gain_ungated)))
}
}

# --- gained / lost detail ----------------------------------------------------
delta_detail <- function(cb, ci, pr, tier) {
  a <- well_quantified(cb, "ungated")[, .(k = paste(FeatureID, SampleID), FeatureID, SampleID)]
  b <- well_quantified(ci, tier)[, .(k = paste(FeatureID, SampleID), FeatureID, SampleID,
                                     CV_after = CV, frac_imputed, days_fully_imputed)]
  cbk <- cb[, .(k = paste(FeatureID, SampleID), CV_before = CV,
                CompleteDays_before = CompleteDays, MeasuredDays_before = MeasuredDays)]
  g <- merge(b[!k %in% a$k], cbk, by = "k", all.x = TRUE)[, `:=`(Protease = pr, change = "gained")]
  l <- merge(a[!k %in% b$k], cbk, by = "k", all.x = TRUE)[, `:=`(Protease = pr, change = "lost")]
  rbindlist(list(g, l), fill = TRUE)
}
detail <- rbindlist(list(delta_detail(cv_try_b, cv_try_i, "Trypsin",  "strict"),
                         delta_detail(cv_hta_b, cv_hta_i, "Krakatoa", "strict")), fill = TRUE)
detail <- merge(detail, genes, by = "FeatureID", all.x = TRUE)
fwrite(detail[order(Protease, change, -frac_imputed)], file.path(TBL, "panel_c_gained_lost_strict.csv"))
cat(sprintf("\nStrict tier: Krakatoa gained %d / lost %d | Trypsin gained %d / lost %d\n",
            detail[Protease == "Krakatoa" & change == "gained", .N],
            detail[Protease == "Krakatoa" & change == "lost",   .N],
            detail[Protease == "Trypsin"  & change == "gained", .N],
            detail[Protease == "Trypsin"  & change == "lost",   .N]))

# ---------------------------------------------------------------- figure ----
# Violin population matches the published figure: every protein x sample pair with a finite
# inter-day CV. Count boxes are the well-quantified subset.
viol_data <- function(cvb, cvi, pr, tier = "strict") {
  keep_i <- switch(tier, strict = cvi$days_fully_imputed == 0,
                          onefill = cvi$days_fully_imputed <= 1,
                          partial = cvi$days_fully_imputed <= 2,
                          ungated = rep(TRUE, nrow(cvi)))
  rbind(
    cvb[is.finite(CV), .(SampleID, CV, cond = "Before")],
    cvi[keep_i & is.finite(CV), .(SampleID, CV, cond = "After")]
  )[, Protease := pr][]
}

panel_row <- function(cvb, cvi, pr, col, light, ylab, tier = "strict", cm = "days") {
  d <- viol_data(cvb, cvi, pr, tier)
  d[, SampleID := factor(SampleID, levels = sample_order)]
  d[, cond := factor(cond, levels = c("Before", "After"))]
  n <- tiers[Protease == pr & Completeness == cm]
  lab <- rbind(
    data.table(SampleID = n$SampleID, cond = "Before", N = n$Baseline),
    data.table(SampleID = n$SampleID, cond = "After",
               N = if (tier == "strict") n$Imputed_strict
                   else if (tier == "onefill") n$Imputed_onefill
                   else if (tier == "partial") n$Imputed_partial else n$Imputed_ungated)
  )
  lab[, SampleID := factor(SampleID, levels = sample_order)]
  lab[, cond := factor(cond, levels = c("Before", "After"))]

  ggplot(d, aes(x = SampleID, y = CV, fill = cond)) +
    geom_violin(trim = TRUE, scale = "width", color = "black",
                linewidth = 0.3, alpha = 0.9, position = position_dodge(width = 0.85)) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", color = "black",
                 linewidth = 0.25, position = position_dodge(width = 0.85)) +
    geom_hline(yintercept = CV_THRESHOLD, linetype = "dashed", linewidth = 0.5) +
    geom_label(data = lab, aes(x = SampleID, y = 152, label = N, fill = cond),
               position = position_dodge(width = 0.85), inherit.aes = FALSE,
               color = "white", fontface = "bold", size = 3.1,
               label.size = 0, label.r = unit(2, "pt"), label.padding = unit(2.2, "pt"),
               show.legend = FALSE) +
    scale_fill_manual(values = c(Before = light, After = col),
                      labels = c(Before = "Measured only",
                                 After  = "After cross-protease imputation")) +
    scale_y_continuous(limits = c(0, 162), breaks = seq(0, 150, 50), expand = c(0, 0)) +
    labs(x = NULL, y = ylab, fill = NULL) +
    theme_classic(base_size = 12) +
    theme(legend.position = "top",
          axis.text = element_text(color = "black"),
          axis.title.y = element_text(face = "bold", size = 11),
          plot.margin = margin(4, 6, 4, 6))
}

cat("\n--- Figures ---\n")
for (cm in c("days", "runs")) {
for (tier in c("strict", "onefill", "partial", "ungated")) {
  p_h <- panel_row(cv_hta_b, cv_hta_i, "Krakatoa", HTA_COL, HTA_LIGHT,
                   "Krakatoa inter-day %CV", tier, cm)
  p_t <- panel_row(cv_try_b, cv_try_i, "Trypsin", TRY_COL, TRY_LIGHT,
                   "Trypsin inter-day %CV", tier, cm) +
         labs(x = "Sample") + theme(axis.title.x = element_text(size = 11))
  p <- (p_h / p_t) +
    plot_annotation(
      title = "C) Cross-protease imputation",
      subtitle = sprintf(
        "Well-quantified proteins: inter-day CV <%d%% and %s. Gaps filled from the other protease; measured values unchanged.\nTier '%s': %s",
        CV_THRESHOLD,
        if (cm == "runs") "ALL 24 injections detected" else "all 3 days detected",
        tier, switch(tier,
          strict  = "every day keeps >=1 real measurement from this protease (no day may be wholly reconstructed)",
          onefill = "at most one day wholly reconstructed from the other protease",
          partial = "at most two days wholly reconstructed",
          ungated = "no restriction - includes proteins with no real measurement at all (upper bound only)")),
      theme = theme(plot.title = element_text(face = "bold", size = 14),
                    plot.subtitle = element_text(size = 9, color = "grey30")))
  nm <- file.path(FIG, paste0("panel_c_imputation_", cm, "_", tier))
  ggsave(paste0(nm, ".png"), p, width = 10, height = 7.5, dpi = PUB_DPI, bg = "white")
  ggsave(paste0(nm, ".pdf"), p, width = 10, height = 7.5, bg = "white")
  cat(sprintf("  Saved: panel_c_imputation_%s_%s.png + .pdf\n", cm, tier))
}
}

cat("\n========================================\nDONE.\n")
cat("Figures: ", FIG, "\nTables:  ", TBL, "\n", sep = "")
cat("========================================\n")
