# =============================================================================
# discovery_imputation.R — cross-protease imputation on the DISCOVERY cohort
#
# The reproducibility cohort measured PRECISION (inter-day CV over 8 injections
# x 3 days). The discovery cohort has NO day/replicate structure — every run is
# a different plasma sample — so inter-day CV is undefined and Panel A/C cannot
# be reused. Discovery measures the other axis: DETECTION / COMPLETENESS.
#
# METRIC (chosen because CV is unavailable):
#   completeness(protein) = fraction of a condition's runs in which the protein
#   is quantified. "Well quantified across the cohort" = completeness >= tau.
#   Counts are reported per protease and as Combined = UNION (a protein counts
#   if EITHER protease quantifies it well; the two proteases are NEVER averaged).
#
# IMPUTATION = cross-protease gap fill only, learned WITHIN each condition:
#   missing cells in one search are filled from the other search of the SAME run,
#   rescaled by a per-protein median log2 offset (>= MIN_PAIRED co-observed runs).
#   Measured values are never altered.
#
# GATING of the "after" count (honesty about pure imputation):
#   strict  — the protein must have >=1 real measurement in that channel/condition
#             (a channel that never saw the protein is never credited).
#   ungated — no restriction; includes proteins with zero real data in a channel
#             (UPPER BOUND only; never "quantified by" that protease).
#
# GUARDRAIL: this script first re-derives the reproducibility Panel A counts and
# ASSERTS the exact match (Trypsin 120/107/127/123/123/123/146/131, Krakatoa
# 92/80/88/89/108/101/91/82) before computing anything on discovery, so any drift
# in the shared load/CV machinery fails loudly.
#
# Usage: Rscript figure-scripts/06_supplementary/discovery_imputation.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

.sd <- function() { a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f)) normalizePath(dirname(f)) else getwd() }
source(file.path(.sd(), "..", "config.R"))   # DATA_DIR, OUT_ROOT, step_dirs()
D <- step_dirs("06_supplementary"); FIG <- D$figures; TBL <- D$tables
PG_TRY <- file.path(SEARCH, "Trypsin_Library_Search",  "trypsin_library_search_report.pg_matrix.tsv")
PG_KRA <- file.path(SEARCH, "Krakatoa_Library_Search", "report.pg_matrix.tsv")
REPRO  <- file.path(SEARCH, "Reproducibility_Search")
RP_TRY <- file.path(REPRO, "Trypsin Library",  "report.pg_matrix.tsv")
RP_HTA <- file.path(REPRO, "Krakatoa Library", "report.pg_matrix.tsv")

ANN <- c("Genes","Protein.Group","Protein.Ids","Protein.Names",
         "First.Protein.Description","N.Sequences","N.Proteotypic.Sequences")
MIN_PAIRED <- 5L                      # runs needed to estimate a protein's offset
TAUS       <- c(0.50, 0.75, 0.90, 1.00)
TAU_MAIN   <- 0.50
PUB_DPI    <- 600
CONDS      <- c("Mix30", "Mix20")

TRY_COL <- "#3C7DA6"; TRY_LIGHT <- "#AECBDE"
HTA_COL <- "#E1812C"; HTA_LIGHT <- "#F4CDA4"
CMB_COL <- "#4C4C4C"; CMB_LIGHT <- "#B8B8B8"

# ----------------------------------------------------------------- loaders ---
# common DIA-NN pg_matrix melt: zero/neg -> NA, duplicate ids summed, empties dropped
melt_pg <- function(dt) {
  rc <- setdiff(names(dt), ANN)
  L <- melt(dt, id.vars = "Protein.Group", measure.vars = rc,
            variable.name = "Sample", value.name = "Intensity", variable.factor = FALSE)
  L[, Intensity := suppressWarnings(as.numeric(Intensity))]
  L[!is.finite(Intensity) | Intensity <= 0, Intensity := NA_real_]
  L <- L[, .(Intensity = if (all(is.na(Intensity))) NA_real_ else sum(Intensity, na.rm = TRUE)),
         by = .(FeatureID = Protein.Group, Sample)]
  L[!is.na(FeatureID) & FeatureID != ""]
}

# discovery: keep run + assign condition from filename
load_disc <- function(path, search) {
  L <- melt_pg(fread(path, sep = "\t", header = TRUE, showProgress = FALSE))
  b <- sub("\\.raw$", "", basename(gsub("\\\\", "/", L$Sample)), ignore.case = TRUE)
  is20 <- grepl("20min", b)
  cond <- rep(NA_character_, length(b))
  cond[grepl("Plasma_Mix", b) &  is20] <- "Mix20"
  cond[grepl("Plasma_Mix", b) & !is20] <- "Mix30"
  cond[grepl("Plasma_Try", b)] <- "Try"
  cond[grepl("Plasma_Kra", b)] <- "Kra"
  cond[grepl("KreimerSimion", b)] <- "BOB"
  cond[is.na(cond)] <- "OTHER"
  L[, Condition := cond][, Search := search][]
}

# reproducibility: Day_SampleID_Replicate from filename (as in 04_reproducibility/Reproducibility_Data_Processing.R)
load_repro <- function(path, search) {
  L <- melt_pg(fread(path, sep = "\t", header = TRUE, showProgress = FALSE))
  b <- sub("\\.raw$", "", basename(gsub("\\\\", "/", L$Sample)), ignore.case = TRUE)
  p <- tstrsplit(b, "_", fixed = TRUE); n <- length(p)
  L[, Day := as.integer(p[[n - 2]])]
  L[, SampleID := p[[n - 1]]]
  L[, Replicate := as.integer(p[[n]])]
  L[, Search := search][]
}

# ------------------------------------------------ Panel A guardrail (repro) ---
inter_cv <- function(L) {
  dm <- L[is.finite(Intensity) & Intensity > 0,
          .(DayMean = mean(Intensity), n = .N), by = .(FeatureID, SampleID, Day)]
  ed <- L[, .(ExpectedDays = uniqueN(Day)), by = SampleID]
  out <- dm[, .(CV = if (.N < 2) NA_real_ else 100 * sd(DayMean) / mean(DayMean),
                MeasuredDays = .N), by = .(FeatureID, SampleID)]
  out <- merge(out, ed, by = "SampleID")
  out[, CompleteDays := MeasuredDays == ExpectedDays][]
}

assert_panel_a <- function() {
  PANEL_A_TRY <- c(120L,107L,127L,123L,123L,123L,146L,131L)
  PANEL_A_HTA <- c(92L,80L,88L,89L,108L,101L,91L,82L)
  tL <- load_repro(RP_TRY, "Trypsin"); hL <- load_repro(RP_HTA, "Krakatoa")
  so <- as.character(sort(as.integer(unique(tL$SampleID))))
  cnt <- function(L) {
    cv <- inter_cv(L)
    w <- cv[CompleteDays %in% TRUE & MeasuredDays >= 2 & !is.na(CV) & CV <= 20, .N, by = SampleID]
    v <- w$N[match(so, w$SampleID)]; v[is.na(v)] <- 0L; as.integer(v)
  }
  bt <- cnt(tL); bh <- cnt(hL)
  ok <- identical(bt, PANEL_A_TRY) && identical(bh, PANEL_A_HTA)
  cat("--- Panel A reproduction guardrail (reproducibility cohort) ---\n")
  cat(sprintf("  Trypsin  : %s  %s\n", paste(bt, collapse=" "), if (identical(bt,PANEL_A_TRY)) "MATCH" else "MISMATCH"))
  cat(sprintf("  Krakatoa : %s  %s\n", paste(bh, collapse=" "), if (identical(bh,PANEL_A_HTA)) "MATCH" else "MISMATCH"))
  if (!ok) stop("Panel A did not reproduce — shared load/CV machinery has drifted; discovery results are not trustworthy.")
  cat("  GUARDRAIL PASSED — proceeding to discovery.\n\n")
}

# ------------------------------------------------------------- imputation ----
# per-protein median log2 offset from co-observed runs WITHIN a condition
make_offsets <- function(target, donor) {
  m <- merge(target[, .(FeatureID, Sample, tI = Intensity)],
             donor[,  .(FeatureID, Sample, dI = Intensity)], by = c("FeatureID","Sample"))
  m <- m[is.finite(tI) & tI > 0 & is.finite(dI) & dI > 0]
  m[, d := log2(tI) - log2(dI)]
  o <- m[, .(offset = median(d), offset_mad = mad(d), n_paired = .N), by = FeatureID]
  o[n_paired >= MIN_PAIRED]
}

impute <- function(target, donor, off) {
  m <- merge(target, donor[, .(FeatureID, Sample, dI = Intensity)],
             by = c("FeatureID","Sample"), all.x = TRUE)
  m <- merge(m, off[, .(FeatureID, offset)], by = "FeatureID", all.x = TRUE)
  m[, imputed := is.na(Intensity) & is.finite(dI) & dI > 0 & !is.na(offset)]
  m[imputed == TRUE, Intensity := 2^(log2(dI) + offset)]
  m[, c("dI","offset") := NULL][]
}

# --------------------------------------------------------- completeness ------
# per-protein completeness within a condition, before and after imputation.
# has_real flags whether the channel ever saw the protein (for the strict gate).
completeness <- function(L) {
  nr <- uniqueN(L$Sample)
  L[, .(n_real = sum(is.finite(Intensity) & !imputed),
        n_all  = sum(is.finite(Intensity)),
        n_imp  = sum(imputed & is.finite(Intensity))),
    by = FeatureID][, `:=`(comp_real = n_real / nr, comp_all = n_all / nr,
                           has_real = n_real > 0)][]
}

# well-quantified feature set at a threshold; strict gate requires >=1 real obs
wq_set <- function(cmp, tau, phase = "before", tier = "strict") {
  if (phase == "before") return(cmp[comp_real >= tau, FeatureID])
  keep <- if (tier == "strict") cmp$has_real else rep(TRUE, nrow(cmp))
  cmp[comp_all >= tau & keep, FeatureID]
}

# ==================================================================== run ====
cat("========================================\n")
cat("DISCOVERY cross-protease imputation\n")
cat("========================================\n\n")
assert_panel_a()

cat("Loading discovery matrices...\n")
tr_all <- load_disc(PG_TRY, "Trypsin"); tr_all[, imputed := FALSE]
kr_all <- load_disc(PG_KRA, "Krakatoa"); kr_all[, imputed := FALSE]
genes <- unique(rbind(
  fread(PG_TRY, select = c("Protein.Group","Genes"))[, .(FeatureID = Protein.Group, Gene = sub(";.*","",Genes))],
  fread(PG_KRA, select = c("Protein.Group","Genes"))[, .(FeatureID = Protein.Group, Gene = sub(";.*","",Genes))]
))[FeatureID != ""]
cat(sprintf("  Trypsin %d proteins | Krakatoa %d proteins | runs %d\n\n",
            uniqueN(tr_all$FeatureID), uniqueN(kr_all$FeatureID), uniqueN(tr_all$Sample)))

# ---- per-condition imputation + completeness --------------------------------
count_rows <- list(); dist_rows <- list(); impute_stats <- list()
CMP <- list()   # keep completeness tables for the figure
for (cond in CONDS) {
  tr <- tr_all[Condition == cond]; kr <- kr_all[Condition == cond]
  nr <- uniqueN(tr$Sample)

  off_tr <- make_offsets(tr, kr)   # trypsin filled from krakatoa
  off_kr <- make_offsets(kr, tr)   # krakatoa filled from trypsin
  tr_i <- impute(tr, kr, off_tr)
  kr_i <- impute(kr, tr, off_kr)

  impute_stats[[cond]] <- data.table(
    Condition = cond, runs = nr,
    offsets_trypsin = nrow(off_tr), offsets_krakatoa = nrow(off_kr),
    offset_mad_med = median(c(off_tr$offset_mad, off_kr$offset_mad), na.rm = TRUE),
    filled_trypsin = tr_i[imputed == TRUE, .N], missing_trypsin = tr[is.na(Intensity), .N],
    filled_krakatoa = kr_i[imputed == TRUE, .N], missing_krakatoa = kr[is.na(Intensity), .N])

  ct_b <- completeness(tr); ct_a <- completeness(tr_i)
  ck_b <- completeness(kr); ck_a <- completeness(kr_i)
  CMP[[cond]] <- list(tr_b = ct_b, tr_a = ct_a, kr_b = ck_b, kr_a = ck_a, nr = nr)

  for (tau in TAUS) {
    T_b <- wq_set(ct_b, tau, "before");                 K_b <- wq_set(ck_b, tau, "before")
    T_s <- wq_set(ct_a, tau, "after", "strict");        K_s <- wq_set(ck_a, tau, "after", "strict")
    T_u <- wq_set(ct_a, tau, "after", "ungated");       K_u <- wq_set(ck_a, tau, "after", "ungated")
    count_rows[[length(count_rows)+1L]] <- data.table(
      Condition = cond, tau = tau,
      Trypsin_before = length(T_b), Trypsin_strict = length(T_s), Trypsin_ungated = length(T_u),
      Krakatoa_before = length(K_b), Krakatoa_strict = length(K_s), Krakatoa_ungated = length(K_u),
      Combined_before = length(union(T_b, K_b)),
      Combined_strict = length(union(T_s, K_s)),
      Combined_ungated = length(union(T_u, K_u)),
      best_single_before = max(length(T_b), length(K_b)))
  }
}
counts <- rbindlist(count_rows)
istat  <- rbindlist(impute_stats)
counts[, `:=`(union_gain_over_single = Combined_before - best_single_before,
              imputation_gain_combined = Combined_strict - Combined_before)]
fwrite(counts, file.path(TBL, "discovery_completeness_counts.csv"))
fwrite(istat,  file.path(TBL, "discovery_imputation_stats.csv"))

cat("--- imputation fill statistics ---\n"); print(istat)
cat("\n--- well-quantified protein counts (Combined = union; proteases never averaged) ---\n")
print(counts[, .(Condition, tau,
                 Tryp_b = Trypsin_before, Tryp_imp = Trypsin_strict,
                 Krak_b = Krakatoa_before, Krak_imp = Krakatoa_strict,
                 Comb_b = Combined_before, Comb_imp = Combined_strict,
                 unionGain = union_gain_over_single, impGain = imputation_gain_combined)])

# ---- decomposition at TAU_MAIN, echoing the reproducibility framing ----------
cat(sprintf("\n--- benefit decomposition at tau>=%.0f%% (strict imputation) ---\n", TAU_MAIN*100))
for (cond in CONDS) {
  r <- counts[Condition == cond & tau == TAU_MAIN]
  cat(sprintf("  %s: Krakatoa alone %d | Trypsin alone %d | UNION %d (+%d over best single) | + imputation %d (+%d)\n",
              cond, r$Krakatoa_before, r$Trypsin_before, r$Combined_before,
              r$union_gain_over_single, r$Combined_strict, r$imputation_gain_combined))
}

# =================================================================== figure ===
# Visual grammar mapped from the reproducibility panel:
#   samples -> conditions (Mix30, Mix20);  inter-day CV -> per-protein completeness
#   threshold line at tau;  count labels = well-quantified proteins per protease + union.
comp_violin_data <- function(cmp_b, cmp_a, cond) {
  rbind(
    data.table(Condition = cond, cond = "Before", comp = 100 * cmp_b$comp_real),
    data.table(Condition = cond, cond = "After",  comp = 100 * cmp_a$comp_all[cmp_a$has_real])
  )
}

build_row <- function(which, col, light, ylab) {
  d <- rbindlist(lapply(CONDS, function(cd)
    comp_violin_data(CMP[[cd]][[paste0(which,"_b")]], CMP[[cd]][[paste0(which,"_a")]], cd)))
  d[, Condition := factor(Condition, levels = CONDS)]
  d[, cond := factor(cond, levels = c("Before","After"))]
  pr <- if (which == "tr") "Trypsin" else "Krakatoa"
  lab <- rbindlist(lapply(CONDS, function(cd) {
    r <- counts[Condition == cd & tau == TAU_MAIN]
    data.table(Condition = cd,
               cond = c("Before","After"),
               N = if (pr == "Trypsin") c(r$Trypsin_before, r$Trypsin_strict)
                   else c(r$Krakatoa_before, r$Krakatoa_strict))
  }))
  lab[, Condition := factor(Condition, levels = CONDS)]
  lab[, cond := factor(cond, levels = c("Before","After"))]
  ggplot(d, aes(Condition, comp, fill = cond)) +
    geom_violin(trim = TRUE, scale = "width", color = "black", linewidth = 0.3,
                alpha = 0.9, position = position_dodge(0.8)) +
    geom_boxplot(width = 0.11, outlier.shape = NA, fill = "white", color = "black",
                 linewidth = 0.25, position = position_dodge(0.8)) +
    geom_hline(yintercept = 100 * TAU_MAIN, linetype = "dashed", linewidth = 0.5) +
    geom_label(data = lab, aes(x = Condition, y = 112, label = N, fill = cond),
               position = position_dodge(0.8), inherit.aes = FALSE, color = "white",
               fontface = "bold", size = 3.4, label.size = 0, label.r = unit(2,"pt"),
               label.padding = unit(2.2,"pt"), show.legend = FALSE) +
    scale_fill_manual(values = c(Before = light, After = col),
                      labels = c(Before = "Measured only",
                                 After = "After cross-protease imputation")) +
    scale_y_continuous(limits = c(0, 120), breaks = seq(0, 100, 25), expand = c(0, 0)) +
    labs(x = NULL, y = ylab, fill = NULL) +
    theme_classic(base_size = 12) +
    theme(legend.position = "top", axis.text = element_text(color = "black"),
          axis.title.y = element_text(face = "bold", size = 11),
          plot.margin = margin(4, 6, 4, 6))
}

# combined-count strip (union), before/after, per condition
comb_strip <- function() {
  d <- rbindlist(lapply(CONDS, function(cd) {
    r <- counts[Condition == cd & tau == TAU_MAIN]
    data.table(Condition = cd, x = c(-0.2, 0.2),
               lab = c(r$Combined_before, r$Combined_strict),
               f = c(CMB_LIGHT, CMB_COL))
  }))
  d[, xpos := as.integer(factor(Condition, levels = CONDS)) + x]
  ggplot(d, aes(xpos, 1)) +
    geom_tile(aes(fill = f), width = 0.36, height = 0.86) +
    geom_text(aes(label = lab), color = "white", fontface = "bold", size = 4.0) +
    scale_fill_identity() +
    scale_x_continuous(limits = c(0.5, length(CONDS) + 0.5), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void() + theme(plot.margin = margin(1, 2, 1, 2))
}
title_row <- function(txt, size = 11.5) {
  ggplot() + annotate("text", 0, 0, label = txt, fontface = "bold", size = size/2.6) +
    theme_void() + theme(plot.margin = margin(4, 2, 1, 2))
}

p_k <- build_row("kr", HTA_COL, HTA_LIGHT, "Krakatoa completeness (% of runs)")
p_t <- build_row("tr", TRY_COL, TRY_LIGHT, "Trypsin completeness (% of runs)")
p_c <- comb_strip()

panel <- (title_row("Cross-protease imputation on the discovery cohort", 15) /
          p_k / title_row(sprintf("Combined proteins (>=%.0f%% completeness), union of the two proteases", TAU_MAIN*100)) /
          p_c / p_t) +
  plot_layout(heights = c(0.28, 3.0, 0.26, 0.34, 3.0)) +
  plot_annotation(
    caption = sprintf(paste0(
      "Well-quantified = protein detected in >=%.0f%% of a condition's 54 runs; violins show the per-protein\n",
      "completeness distribution, counts are proteins above the threshold. After = cross-protease gap fill\n",
      "(per-protein median log2 offset, >=%d paired runs); strict gate requires >=1 real measurement in that\n",
      "protease. Measured values unchanged; the two proteases are never averaged."),
      TAU_MAIN*100, MIN_PAIRED),
    theme = theme(plot.caption = element_text(hjust = 0, size = 8.4, color = "grey30",
                                              margin = margin(6, 0, 0, 0))))
ggsave(file.path(FIG, "panel_discovery_imputation.png"), panel, width = 7.6, height = 8.6, dpi = PUB_DPI, bg = "white")
ggsave(file.path(FIG, "panel_discovery_imputation.pdf"), panel, width = 7.6, height = 8.6, bg = "white")
cat("\nSaved figure: panel_discovery_imputation.png + .pdf\n")

# ---- threshold-sweep line figure (counts vs tau) ----------------------------
sw <- melt(counts[, .(Condition, tau, Trypsin = Trypsin_before, Krakatoa = Krakatoa_before,
                      Combined = Combined_before, `Combined+imp` = Combined_strict)],
           id.vars = c("Condition","tau"), variable.name = "series", value.name = "N")
sw[, series := factor(series, levels = c("Krakatoa","Trypsin","Combined","Combined+imp"))]
p_sw <- ggplot(sw, aes(100*tau, N, color = series, linetype = series)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.6) +
  facet_wrap(~Condition) +
  scale_color_manual(values = c(Krakatoa = HTA_COL, Trypsin = TRY_COL,
                                Combined = CMB_COL, `Combined+imp` = "#1b7837")) +
  scale_linetype_manual(values = c(Krakatoa = 1, Trypsin = 1, Combined = 1, `Combined+imp` = 2)) +
  labs(x = "Completeness threshold (% of runs)", y = "Well-quantified proteins",
       color = NULL, linetype = NULL,
       title = "Well-quantified proteins vs completeness threshold") +
  theme_classic(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(FIG, "threshold_sweep.png"), p_sw, width = 9, height = 4.5, dpi = PUB_DPI, bg = "white")
cat("Saved figure: threshold_sweep.png\n")

cat("\n========================================\nDONE.\n")
cat("Figures: ", FIG, "\nTables:  ", TBL, "\n", sep = "")
cat("========================================\n")
