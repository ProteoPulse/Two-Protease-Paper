# =============================================================================
# masking_validation.R — is cross-protease fill trustworthy on DISCOVERY?
#
# Hide known real values, impute them, compare to the withheld truth. Two regimes:
#   MCAR  a random fraction of co-observed recipient cells — the flattering
#         regime (masked values are ordinary), where the recipient's own median
#         is hard to beat. Reported, but does not decide the question.
#   MNAR  the lowest-abundance fraction of each protein's co-observed cells —
#         the left-censored analogue of why proteomics values are actually
#         missing. This is the regime that matches real gaps and decides it.
#
# Cross-protease (per-protein median log2 offset) is benchmarked against the
# standard fallbacks: sample-median, MinDet and MinProb (the classic downshifted
# MNAR draws), protein-mean and protein-kNN. A shuffled-donor negative control
# checks the method actually uses the donor's per-sample signal.
#
# Both directions are run: recipient = Krakatoa (donor Trypsin) is the one that
# matters on discovery (Krakatoa is the sparse channel that gets filled), and
# recipient = Trypsin (donor Krakatoa) for symmetry. Mix30 and Mix20 both.
#
# Per-cell MNAR errors are additionally stratified by intensity quintile and by
# donor completeness (how well the donor sees that protein), the two axes that
# matter for whether a fill can be believed.
#
# Reuses the imputers from protease-quant-advisor/vendor/benchmark_imputation.R.
# Usage: Rscript figure-scripts/06_supplementary/masking_validation.R
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

.sd <- function() { a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f)) normalizePath(dirname(f)) else getwd() }
source(file.path(.sd(), "..", "config.R"))   # DATA_DIR, OUT_ROOT, step_dirs()
D <- step_dirs("06_supplementary"); FIG <- D$figures; TBL <- D$tables

source(file.path(ADVISOR_DIR, "vendor", "benchmark_imputation.R"))
PG_TRY <- file.path(SEARCH, "Trypsin_Library_Search",  "trypsin_library_search_report.pg_matrix.tsv")
PG_KRA <- file.path(SEARCH, "Krakatoa_Library_Search", "report.pg_matrix.tsv")
ANN <- c("Genes","Protein.Group","Protein.Ids","Protein.Names",
         "First.Protein.Description","N.Sequences","N.Proteotypic.Sequences")
MIN_PAIRED <- 8L
SEED <- 20260909

# ----------------------------------------------------------------- load ------
load_disc <- function(path) {
  dt <- fread(path, sep = "\t", header = TRUE, showProgress = FALSE)
  rc <- setdiff(names(dt), ANN)
  L <- melt(dt, id.vars = "Protein.Group", measure.vars = rc,
            variable.name = "Sample", value.name = "Intensity", variable.factor = FALSE)
  L[, Intensity := suppressWarnings(as.numeric(Intensity))]
  L[!is.finite(Intensity) | Intensity <= 0, Intensity := NA_real_]
  L <- L[, .(Intensity = if (all(is.na(Intensity))) NA_real_ else sum(Intensity, na.rm = TRUE)),
         by = .(pg_canonical = Protein.Group, Sample)]
  L <- L[!is.na(pg_canonical) & pg_canonical != ""]
  b <- sub("\\.raw$", "", basename(gsub("\\\\", "/", L$Sample)), ignore.case = TRUE)
  is20 <- grepl("20min", b)
  cond <- fifelse(grepl("Plasma_Mix", b) & is20, "Mix20",
          fifelse(grepl("Plasma_Mix", b) & !is20, "Mix30",
          fifelse(grepl("Plasma_Try", b), "Try",
          fifelse(grepl("Plasma_Kra", b), "Kra",
          fifelse(grepl("KreimerSimion", b), "BOB", "OTHER")))))
  L[, `:=`(sample_id = b, Condition = cond)][]
}
TR <- load_disc(PG_TRY); KR <- load_disc(PG_KRA)

# X for one (condition, recipient direction): theta_rec / theta_don in log2
make_X <- function(cond, rec) {
  r <- if (rec == "Krakatoa") KR else TR
  d <- if (rec == "Krakatoa") TR else KR
  r <- r[Condition == cond, .(pg_canonical, sample_id, theta_rec = log2(Intensity))]
  d <- d[Condition == cond, .(pg_canonical, sample_id, theta_don = log2(Intensity))]
  X <- merge(r, d, by = c("pg_canonical","sample_id"), all = TRUE)
  X[is.finite(theta_rec) | is.finite(theta_don)]
}

cfg_b <- list(regimes = c("mcar","mnar"), mcar_repeats = 20L, mask_frac = 0.30,
              baselines = c("protein_mean","sample_median","mindet","minprob","protein_knn"),
              minprob_shift_sd = 1.8, minprob_width_sd = 0.3, negative_control = TRUE)

# ============================ aggregate benchmark ============================
grid <- CJ(Condition = c("Mix30","Mix20"), Recipient = c("Krakatoa","Trypsin"), sorted = FALSE)
cat("========================================\n")
cat("MASKING VALIDATION — discovery cohort\n")
cat("========================================\n")
agg <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
  cond <- grid$Condition[i]; rec <- grid$Recipient[i]
  X <- make_X(cond, rec)
  cat(sprintf("\n--- %s | recipient=%s donor=%s | maskable co-observed cells: %d ---\n",
              cond, rec, if (rec=="Krakatoa") "Trypsin" else "Krakatoa",
              X[is.finite(theta_rec) & is.finite(theta_don), .N]))
  b <- run_benchmark(X, cfg_b, min_paired = MIN_PAIRED, seed = SEED)
  s <- b$summary[, .(Condition = cond, Recipient = rec, regime, method,
                     n_cells, median_abs_error, rmse, median_signed_error, spearman)]
  print(s[regime %in% c("mnar","mcar")][order(regime, median_abs_error)])
  s
}))
fwrite(agg, file.path(TBL, "masking_aggregate_discovery.csv"))

# ================= per-cell MNAR errors, stratified ==========================
# recompute the MNAR mask ourselves so we can attach intensity + donor completeness to each cell
strat_one <- function(cond, rec) {
  X <- make_X(cond, rec)
  # donor completeness per protein (fraction of samples where donor observed)
  ns <- uniqueN(X$sample_id)
  dcomp <- X[, .(don_comp = sum(is.finite(theta_don)) / ns), by = pg_canonical]
  mi <- mask_mnar(X, cfg_b$mask_frac)
  if (!length(mi)) return(NULL)
  truth <- X$theta_rec[mi]
  Y <- copy(X); Y[, train := TRUE]; Y[mi, train := FALSE]
  Yh <- copy(Y); Yh[mi, theta_rec := NA_real_]
  preds <- list(
    cross_protease = .imp_cross_protease(Yh, Yh$train, mi, MIN_PAIRED),
    sample_median  = .imp_sample_median(Yh, Yh$train, mi),
    mindet         = .imp_mindet(Yh, Yh$train, mi),
    minprob        = .imp_minprob(Yh, Yh$train, mi, cfg_b$minprob_shift_sd, cfg_b$minprob_width_sd, SEED))
  cell <- data.table(Condition = cond, Recipient = rec,
                     pg_canonical = X$pg_canonical[mi], sample_id = X$sample_id[mi],
                     truth = truth, theta_don = X$theta_don[mi])
  cell <- merge(cell, dcomp, by = "pg_canonical", all.x = TRUE)
  for (nm in names(preds)) cell[[paste0("err_", nm)]] <- abs(preds[[nm]] - truth)
  cell
}
cells <- rbindlist(lapply(seq_len(nrow(grid)),
                          function(i) strat_one(grid$Condition[i], grid$Recipient[i])), fill = TRUE)
fwrite(cells, file.path(TBL, "masking_cells_mnar_discovery.csv"))

emeth <- c("err_cross_protease","err_sample_median","err_mindet","err_minprob")
long <- melt(cells, id.vars = c("Condition","Recipient","truth","don_comp"),
             measure.vars = emeth, variable.name = "method", value.name = "abs_err")
long[, method := factor(sub("^err_", "", method),
                        levels = c("cross_protease","sample_median","mindet","minprob"))]

# stratify by recipient-truth intensity quintile (within condition x recipient)
long[, int_bin := cut(truth, breaks = quantile(truth, probs = seq(0,1,0.2), na.rm = TRUE),
                      include.lowest = TRUE, labels = c("Q1(dim)","Q2","Q3","Q4","Q5(bright)")),
     by = .(Condition, Recipient)]
# stratify by donor completeness
long[, don_bin := cut(don_comp, breaks = c(-.01, .5, .8, .95, 1.001),
                      labels = c("<50%","50-80%","80-95%",">95%"))]

cat("\n=== MNAR per-cell median abs error by method (recipient = Krakatoa, the used direction) ===\n")
print(long[Recipient == "Krakatoa", .(median_abs_err = round(median(abs_err, na.rm=TRUE),3),
                                       n = .N), by = .(Condition, method)][order(Condition, method)])
cat("\n=== stratified by intensity quintile (Krakatoa recipient, both conditions pooled) ===\n")
print(dcast(long[Recipient == "Krakatoa"], int_bin ~ method,
            value.var = "abs_err", fun.aggregate = function(x) round(median(x, na.rm=TRUE),3)))
cat("\n=== stratified by donor completeness (Krakatoa recipient) ===\n")
print(dcast(long[Recipient == "Krakatoa"], don_bin ~ method,
            value.var = "abs_err", fun.aggregate = function(x) round(median(x, na.rm=TRUE),3)))

strat_int <- long[, .(median_abs_err = median(abs_err, na.rm=TRUE), n=.N),
                  by = .(Condition, Recipient, method, int_bin)]
strat_don <- long[, .(median_abs_err = median(abs_err, na.rm=TRUE), n=.N),
                  by = .(Condition, Recipient, method, don_bin)]
fwrite(strat_int, file.path(TBL, "masking_stratified_intensity.csv"))
fwrite(strat_don, file.path(TBL, "masking_stratified_completeness.csv"))

# ------------------------------------------------------------------ figure ---
COLM <- c(cross_protease = "#1b7837", sample_median = "#762a83",
          mindet = "#E1812C", minprob = "#3C7DA6")
pA <- ggplot(agg[regime == "mnar"], aes(method, median_abs_error, fill = method)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  facet_grid(Recipient ~ Condition) +
  scale_fill_manual(values = c(cross_protease = "#1b7837", protein_mean = "#999999",
                               sample_median = "#762a83", mindet = "#E1812C",
                               minprob = "#3C7DA6", protein_knn = "#5aae61"), guide = "none") +
  labs(x = NULL, y = "MNAR median abs error (log2)",
       title = "Masking recovery — MNAR (left-censored) regime",
       subtitle = "Lower is better. Recipient = the channel whose low-abundance cells were hidden and re-imputed.") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, color = "black"))
ggsave(file.path(FIG, "masking_mnar_bars.png"), pA, width = 9, height = 6, dpi = 600, bg = "white")

pB <- ggplot(strat_int[Recipient == "Krakatoa"],
             aes(int_bin, median_abs_err, color = method, group = method)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.8) +
  facet_wrap(~Condition) +
  scale_color_manual(values = COLM) +
  labs(x = "Recipient intensity quintile (masked truth)", y = "MNAR median abs error (log2)",
       color = NULL, title = "Fill accuracy vs intensity (Krakatoa recipient)") +
  theme_classic(base_size = 11) + theme(legend.position = "top",
       axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(FIG, "masking_stratified_intensity.png"), pB, width = 9, height = 4.6, dpi = 600, bg = "white")

cat("\nSaved: masking_mnar_bars.png, masking_stratified_intensity.png\n")
cat("Tables: masking_aggregate_discovery.csv, masking_cells_mnar_discovery.csv,\n")
cat("        masking_stratified_intensity.csv, masking_stratified_completeness.csv\n")
cat("\n========================================\nDONE.\n========================================\n")
