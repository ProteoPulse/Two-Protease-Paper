# =============================================================================
# Completeness: Mix30 (cross-protease deliverable) vs Trypsin 2 (trypsin-only)
# 2026-09-27 — see PLAN.md in this folder.
#
# Paired comparison on the same 54 plasma samples. Trypsin 2 = 251117 DIA_Trypsin
# runs (30-min gradient, verified from report.parquet RT). Primary baseline = the
# trypsin-library search only; union-of-searches is a sensitivity row.
#
# Arms (protein x sample presence, 54 samples):
#   A   Try2        trypsin-library search                       (primary baseline)
#   A2  Try2-union  union of both searches (the convention of 03_ptm_coverage/sequence_coverage_comparison.Rmd)  (sensitivity)
#   B   Mix T-ch    Mix trypsin channel only                     (dilution cost)
#   C   Mix deliv   one channel per protein, chosen by completeness (measured)
#   D   Mix +imp    C + strict cross-protease imputation          (Advisor deliverable)
#   E   Mix ceiling run-level union of both channels             (upper bound of any imputation)
#   R   Try1        250807 BOB trypsin prep (noise reference; not sample-paired)
# Mix30 is primary; Mix20 (20-min mixture) is run through the same arms.
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})
set.seed(20260927)

.sd <- function() { a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f)) normalizePath(dirname(f)) else getwd() }
source(file.path(.sd(), "..", "config.R"))   # DATA_DIR, OUT_ROOT, step_dirs()
D <- step_dirs("07_verification"); FIG <- D$figures; TBL <- D$tables
source(file.path(.sd(), "helpers.R"))

# ====================================================================== run ==
cat("== Completeness: Mix30 vs Trypsin 2 ==\n")
assert_panel_a()
TR <- load_disc(PG_TRY); KR <- load_disc(PG_KRA)

runs <- TR[, .(runs = uniqueN(Run), samples = uniqueN(SampleID)), by = Cond][order(Cond)]
print(runs)
stopifnot(all(runs[Cond %in% c("Try2","Try1","Mix30","Mix20"), runs] == 54L))
S <- sort(unique(TR[Cond == "Try2", SampleID]))
stopifnot(length(S) == 54L, setequal(S, unique(TR[Cond == "Mix30", SampleID])),
          setequal(S, unique(TR[Cond == "Mix20", SampleID])))
S1 <- sort(unique(TR[Cond == "Try1", SampleID]))
cat("Try2/Mix30/Mix20 sample sets identical (54). Try1 = 54 BOB runs (not sample-paired).\n\n")

# ---- arms ----
A  <- pres_mat(TR[Cond == "Try2"], S)
AK <- pres_mat(KR[Cond == "Try2"], S)
pr <- union(rownames(A), rownames(AK)); A2 <- align(A, pr) | align(AK, pr)
R  <- pres_mat(TR[Cond == "Try1"], S1)
MX <- list(Mix30 = build_mix_arms(TR, KR, "Mix30", S), Mix20 = build_mix_arms(TR, KR, "Mix20", S))

# ---- verification against earlier numbers ----
m30 <- MX$Mix30; k50 <- ceiling(0.5 * 54)
v <- c(sep9_T = sum(rowSums(m30$Tp) >= k50), sep9_K = sum(rowSums(m30$Kp) >= k50),
       sep9_union = sum(rowSums(m30$Tp) >= k50 | rowSums(m30$Kp) >= k50),
       sep9_union_imp = sum(rowSums(m30$sep9_union_imp) >= k50 | rowSums(m30$sep9_union_impK) >= k50),
       ms_try2_pg = nrow(A2[rowSums(A2) > 0, ]), ms_try2_perrun_mean = round(mean(colSums(A2)), 1),
       try2_T_tau50 = sum(rowSums(A) >= k50), try2_T_tau100 = sum(rowSums(A) == 54))
print(v)
chk <- c(sep9_T = 1009, sep9_K = 465, sep9_union = 1276, sep9_union_imp = 1278, ms_try2_pg = 1507,
         ms_try2_perrun_mean = 1075.5, try2_T_tau50 = 1053, try2_T_tau100 = 652)
fwrite(data.table(check = names(chk), expected = chk, observed = v[names(chk)], match = v[names(chk)] == chk),
       file.path(TBL, "verification_checks.csv"))
if (!all(v[names(chk)] == chk)) warning("Some verification checks differ — see tables/verification_checks.csv")

# ---- 1. per-sample counts (paired) ----
arm_list <- function(mx) list("Try2" = A, "Try2-union" = A2, "Mix T-ch" = mx$B, "Mix deliv" = mx$C,
                              "Mix +imp" = mx$D, "Mix ceiling" = mx$E)
ps <- rbindlist(lapply(names(MX), function(cd) {
  al <- arm_list(MX[[cd]])
  rbindlist(lapply(names(al), function(a) data.table(Mix = cd, Arm = a, SampleID = S, n = colSums(al[[a]][, S]))))
}))
ps <- unique(ps[!(Arm %in% c("Try2","Try2-union") & Mix == "Mix20")])     # baseline listed once
ps[Arm %in% c("Try2","Try2-union"), Mix := "baseline"]
fwrite(ps, file.path(TBL, "per_sample_protein_counts.csv"))
base <- ps[Arm == "Try2", .(SampleID, n0 = n)]
ps_sum <- merge(ps, base, by = "SampleID")[, .(median_n = median(n), mean_n = round(mean(n), 1),
  median_paired_diff_vs_Try2 = median(n - n0), min_diff = min(n - n0), max_diff = max(n - n0),
  n_samples_higher = sum(n > n0), wilcoxon_p = if (all(n == n0)) NA_real_ else
    suppressWarnings(wilcox.test(n, n0, paired = TRUE)$p.value)), by = .(Mix, Arm)]
ps_sum <- rbind(ps_sum, data.table(Mix = "reference", Arm = "Try1", median_n = median(colSums(R)),
  mean_n = round(mean(colSums(R)), 1), median_paired_diff_vs_Try2 = NA, min_diff = NA, max_diff = NA,
  n_samples_higher = NA, wilcoxon_p = NA))
fwrite(ps_sum, file.path(TBL, "per_sample_summary.csv"))
cat("\n--- per-sample protein counts ---\n"); print(ps_sum)

# ---- 2. threshold curve + bootstrap CIs ----
K <- 54L
curve <- rbindlist(c(
  lapply(names(MX), function(cd) { al <- arm_list(MX[[cd]])
    rbindlist(lapply(names(al), function(a) data.table(Mix = cd, Arm = a, k = 1:K, n = ge_counts(al[[a]]))))
  }),
  list(data.table(Mix = "Mix30", Arm = "Try1", k = 1:K, n = ge_counts(R)),
       data.table(Mix = "Mix20", Arm = "Try1", k = 1:K, n = ge_counts(R)))))
curve[, tau := 100 * k / K]
curve[, Mix := factor(Mix, levels = c("Mix30","Mix20"))]

# paired bootstrap over the 54 samples for Mix arms - Try2; unpaired for Try1 - Try2
boot_diff <- function(M1, M0, paired = TRUE) {
  out <- matrix(NA_real_, N_BOOT, K)
  for (b in seq_len(N_BOOT)) {
    i0 <- sample.int(K, K, TRUE); i1 <- if (paired) i0 else sample.int(K, K, TRUE)
    out[b, ] <- ge_counts(M1, i1) - ge_counts(M0, i0)
  }
  data.table(k = 1:K, lo = apply(out, 2, quantile, 0.025), hi = apply(out, 2, quantile, 0.975))
}
diffs <- rbindlist(c(
  unlist(lapply(names(MX), function(cd) lapply(c("Mix T-ch","Mix deliv","Mix +imp","Mix ceiling"), function(a) {
    M1 <- arm_list(MX[[cd]])[[a]]
    d <- data.table(Mix = cd, Arm = a, k = 1:K, diff = ge_counts(M1) - ge_counts(A))
    merge(d, boot_diff(M1, A, TRUE), by = "k") })), recursive = FALSE),
  list(merge(data.table(Mix = "reference", Arm = "Try1", k = 1:K, diff = ge_counts(R) - ge_counts(A)),
             boot_diff(R, A, FALSE), by = "k"))))
diffs[, tau := 100 * k / K]
fwrite(curve, file.path(TBL, "threshold_curve_counts.csv")); fwrite(diffs, file.path(TBL, "threshold_curve_diff_vs_Try2_bootCI.csv"))

rep_tab <- dcast(curve[k %in% ceiling(TAUS_REP * K) & !(Arm == "Try1" & Mix == "Mix20")],
                 k + tau ~ Mix + Arm, value.var = "n")
fwrite(rep_tab, file.path(TBL, "headline_counts_by_threshold.csv"))
cross <- diffs[, .(first_k_negative = if (any(diff < 0)) min(k[diff < 0]) else NA_integer_,
                   first_k_CI_below0 = if (any(hi < 0)) min(k[hi < 0]) else NA_integer_,
                   last_k_CI_above0 = if (any(lo > 0)) max(k[lo > 0]) else NA_integer_), by = .(Mix, Arm)]
cross[, `:=`(first_tau_negative = round(100 * first_k_negative / K, 1),
             first_tau_CI_below0 = round(100 * first_k_CI_below0 / K, 1),
             last_tau_CI_above0 = round(100 * last_k_CI_above0 / K, 1))]
fwrite(cross, file.path(TBL, "crossover_thresholds.csv"))
cat("\n--- headline counts (proteins with completeness >= tau) ---\n"); print(rep_tab)
cat("\n--- crossover (Mix arm - Try2 turns negative) ---\n"); print(cross)
cat("\n--- diff vs Try2 with 95% CI at headline thresholds ---\n")
print(diffs[k %in% ceiling(TAUS_REP * K), .(Mix, Arm, tau = round(tau), diff, lo, hi)])

# ---- 3. datapoints + matrix completeness ----
dp <- rbindlist(c(
  lapply(names(MX), function(cd) { al <- arm_list(MX[[cd]])
    rbindlist(lapply(names(al), function(a) { M <- al[[a]]; M <- M[rowSums(M) > 0, , drop = FALSE]
      data.table(Mix = cd, Arm = a, proteins = nrow(M), cells = sum(M), matrix_pct = round(100 * mean(M), 1)) })) }),
  list(data.table(Mix = "reference", Arm = "Try1", proteins = nrow(R), cells = sum(R), matrix_pct = round(100 * mean(R), 1)))))
dp <- unique(dp[!(Mix == "Mix20" & Arm %in% c("Try2","Try2-union"))])
dp[Arm %in% c("Try2","Try2-union"), Mix := "baseline"]
fwrite(dp, file.path(TBL, "datapoints_and_matrix_completeness.csv"))
cat("\n--- datapoints ---\n"); print(dp)

# ---- 4. decomposition: Mix deliverable (+imp) vs Try2 ----
decomp <- function(M1, M0, lab) {
  all <- union(rownames(M1), rownames(M0))
  n1 <- setNames(rep(0, length(all)), all); n0 <- n1
  n1[rownames(M1)] <- rowSums(M1); n0[rownames(M0)] <- rowSums(M0)
  cls <- fcase(n0 == 0 & n1 > 0, "new (absent in Try2)", n0 > 0 & n1 == 0, "lost (absent in Mix)",
               n1 > n0, "shared: more complete", n1 < n0, "shared: less complete", default = "shared: unchanged")
  per <- data.table(FeatureID = all, n_try2 = n0, n_mix = n1, class = cls, comparison = lab)
  net <- rbindlist(lapply(TAUS_REP, function(t) { kk <- ceiling(t * K)
    data.table(comparison = lab, tau = 100 * t,
      new_above = sum(n0 == 0 & n1 >= kk), shared_crossed_up = sum(n0 > 0 & n0 < kk & n1 >= kk),
      shared_crossed_down = sum(n1 > 0 & n0 >= kk & n1 < kk), lost_above = sum(n1 == 0 & n0 >= kk),
      net = sum(n1 >= kk) - sum(n0 >= kk)) }))
  list(per = per, net = net)
}
dcs <- list(decomp(MX$Mix30$D, A, "Mix30 +imp vs Try2"), decomp(MX$Mix30$C, A, "Mix30 deliv vs Try2"),
            decomp(MX$Mix30$B, A, "Mix30 T-ch vs Try2"), decomp(MX$Mix20$D, A, "Mix20 +imp vs Try2"),
            decomp(R, A, "Try1 vs Try2 (reference)"))
per <- rbindlist(lapply(dcs, `[[`, "per")); net <- rbindlist(lapply(dcs, `[[`, "net"))
cls_tab <- per[, .N, by = .(comparison, class)][order(comparison, class)]
fwrite(cls_tab, file.path(TBL, "decomposition_classes.csv")); fwrite(net, file.path(TBL, "decomposition_net_by_threshold.csv"))
cat("\n--- decomposition classes ---\n"); print(dcast(cls_tab, class ~ comparison, value.var = "N"))
cat("\n--- net change by threshold ---\n"); print(net)

# how many runs are lost/gained for shared Try2 proteins
dn <- per[n_try2 > 0 & n_mix > 0, .(comparison, delta = n_mix - n_try2)]
dn_sum <- dn[, .(n = .N, median_delta = median(delta), mean_delta = round(mean(delta), 2),
                 frac_less = round(mean(delta < 0), 3), frac_more = round(mean(delta > 0), 3),
                 frac_loss_1to2 = round(mean(delta %in% -2:-1), 3), frac_loss_ge5 = round(mean(delta <= -5), 3)),
             by = comparison]
fwrite(dn_sum, file.path(TBL, "shared_protein_delta_runs.csv")); cat("\n--- shared-protein run deltas ---\n"); print(dn_sum)

# ---- 5. abundance stratification (Try2 median log2 intensity quintile) ----
Ai <- int_mat(TR[Cond == "Try2"], S)
ab <- data.table(FeatureID = rownames(Ai), log2med = apply(log2(Ai), 1, median, na.rm = TRUE))
ab[, quintile := paste0("Q", cut(log2med, quantile(log2med, 0:5/5), labels = FALSE, include.lowest = TRUE))]
strat <- merge(per[n_try2 > 0], ab, by = "FeatureID")[, .(proteins = .N,
  frac_less_complete = round(mean(n_mix < n_try2), 3), frac_more_complete = round(mean(n_mix > n_try2), 3),
  median_delta_runs = median(n_mix - n_try2), mean_try2_completeness = round(100 * mean(n_try2) / K, 1)),
  by = .(comparison, quintile)][order(comparison, quintile)]
fwrite(strat, file.path(TBL, "abundance_stratified_loss.csv")); cat("\n--- by Try2 abundance quintile (Q1 = lowest) ---\n"); print(strat)

# ---- 6. dilution check: Mix30 trypsin channel vs Try2 intensity, same protein & sample ----
Mi <- int_mat(TR[Cond == "Mix30"], S); sh <- intersect(rownames(Ai), rownames(Mi))
lr <- log2(Mi[sh, S]) - log2(Ai[sh, S])
dil <- data.table(median_log2_ratio_all_cells = round(median(lr, na.rm = TRUE), 3),
                  expected_if_80pct_load = round(log2(0.8), 3), co_observed_cells = sum(is.finite(lr)),
                  median_of_per_sample_medians = round(median(apply(lr, 2, median, na.rm = TRUE)), 3))
fwrite(dil, file.path(TBL, "dilution_check_mix30T_vs_try2.csv")); cat("\n--- dilution check ---\n"); print(dil)

# examples: largest losses among well-measured Try2 proteins
genes <- unique(rbind(fread(PG_TRY, select = c("Protein.Group","Genes")), fread(PG_KRA, select = c("Protein.Group","Genes"))))
setnames(genes, "Protein.Group", "FeatureID")
ex <- merge(merge(per[comparison == "Mix30 +imp vs Try2"], ab, by = "FeatureID", all.x = TRUE), genes, by = "FeatureID", all.x = TRUE)
ex[, `:=`(chosen = MX$Mix30$chosen[FeatureID], delta = n_mix - n_try2)]
fwrite(ex[order(delta)][, .(FeatureID, Genes, n_try2, n_mix, delta, class, chosen, log2med, quintile)],
       file.path(TBL, "per_protein_mix30_vs_try2.csv"))
cat("\nFilled cells (strict, chosen channel): Mix30", MX$Mix30$n_filled, "| Mix20", MX$Mix20$n_filled, "\n")

# ---- 7. half-cohort stability (without replacement; the with-replacement bootstrap
#         inflates counts at tau = 100% because duplicated samples are easier to 'always' hit) ----
N_HALF <- 1000L; H <- K %/% 2L
half <- rbindlist(lapply(seq_len(N_HALF), function(i) {
  j0 <- sample.int(K, H); j1 <- sample.int(K, H)          # paired halves use j0; Try1 reference independent
  a0 <- ge_counts(A, j0)
  rbindlist(list(
    data.table(rep = i, Arm = "Mix30 +imp", k = 1:H, diff = ge_counts(MX$Mix30$D, j0) - a0),
    data.table(rep = i, Arm = "Mix20 +imp", k = 1:H, diff = ge_counts(MX$Mix20$D, j0) - a0),
    data.table(rep = i, Arm = "Try1 (reference)", k = 1:H, diff = ge_counts(R, j1) - a0)))
}))
half_sum <- half[k %in% ceiling(TAUS_REP * H), .(median_diff = median(diff), lo = quantile(diff, 0.025),
  hi = quantile(diff, 0.975), frac_mix_below_try2 = round(mean(diff < 0), 3)), by = .(Arm, k)][, tau := round(100 * k / H)]
fwrite(half_sum, file.path(TBL, "half_cohort_stability.csv"))
cat("\n--- half-cohort (27 of 54 samples, 1000 draws) difference vs Try2 ---\n"); print(half_sum)

# ---- 8. evidence depth of Try2 proteins by completeness change (report.parquet, 1% run+PG q) ----
suppressPackageStartupMessages(library(duckdb))
con <- dbConnect(duckdb())
prec <- as.data.table(dbGetQuery(con, sprintf("
  SELECT \"Protein.Group\" AS FeatureID, Run, count(DISTINCT \"Precursor.Id\") AS n_prec
  FROM read_parquet('%s')
  WHERE Run LIKE '251117_DIA_Plasma_Try%%' AND Decoy = 0 AND \"Q.Value\" <= 0.01 AND \"PG.Q.Value\" <= 0.01
  GROUP BY 1, 2", file.path(SEARCH, "Trypsin_Library_Search", "report.parquet"))))
dbDisconnect(con, shutdown = TRUE)
stopifnot(uniqueN(prec$Run) == 54L)
pd <- prec[, .(med_prec = median(n_prec), max_prec = max(n_prec)), by = FeatureID]
ev <- merge(per[comparison %in% c("Mix30 +imp vs Try2", "Try1 vs Try2 (reference)") & n_try2 > 0], pd, by = "FeatureID", all.x = TRUE)
ev[, grp := fcase(n_mix < n_try2, "less complete", n_mix > n_try2, "more complete", default = "unchanged")]
ev_sum <- ev[, .(proteins = .N, median_precursors_per_run = median(med_prec, na.rm = TRUE),
                 frac_single_precursor = round(mean(med_prec <= 1, na.rm = TRUE), 3),
                 mean_try2_completeness = round(100 * mean(n_try2) / K, 1)), by = .(comparison, grp)][order(comparison, grp)]
# heavy losers: lose >= 10 runs
ev_heavy <- ev[comparison == "Mix30 +imp vs Try2" & n_try2 - n_mix >= 10,
               .(proteins = .N, median_precursors_per_run = median(med_prec, na.rm = TRUE),
                 frac_single_precursor = round(mean(med_prec <= 1, na.rm = TRUE), 3))]
fwrite(ev_sum, file.path(TBL, "evidence_depth_by_change.csv"))
cat("\n--- Try2 precursor depth by completeness change ---\n"); print(ev_sum)
cat("Mix30 heavy losers (>= 10 runs lost):\n"); print(ev_heavy)
all_sp <- pd[, round(mean(med_prec <= 1), 3)]
cat("All Try2 proteins, fraction single-precursor:", all_sp, "\n")

# ================================================================= figures ==
ARM_LAB <- c("Try2" = "Trypsin 2 (trypsin lib)", "Try2-union" = "Trypsin 2 (both libs)",
             "Mix T-ch" = "Mix: trypsin channel only", "Mix deliv" = "Mix: chosen channel",
             "Mix +imp" = "Mix: chosen channel + imputation", "Mix ceiling" = "Mix: either channel (ceiling)",
             "Try1" = "Trypsin 1 (other prep)")
th <- theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold", size = 11.5),
  plot.subtitle = element_text(size = 8.5, colour = "grey30"), axis.text = element_text(colour = "black"),
  panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3), legend.key.width = unit(1.4, "lines"),
  strip.background = element_blank(), strip.text = element_text(face = "bold"))

# Fig 1: threshold curves + difference with CI
ARMS_C <- c("Try1","Try2","Mix T-ch","Mix +imp","Mix ceiling")
cd1 <- curve[Arm %in% ARMS_C & tau >= 50]
cd1[, Arm := factor(Arm, levels = ARMS_C)]
p1a <- ggplot(cd1, aes(tau, n, colour = Arm, linetype = Arm)) + geom_line(linewidth = 0.8) +
  facet_wrap(~ Mix) +
  scale_colour_manual(values = COL[ARMS_C], labels = ARM_LAB[ARMS_C], name = NULL) +
  scale_linetype_manual(values = c("Try1" = "solid", "Try2" = "solid", "Mix T-ch" = "22", "Mix +imp" = "solid",
                                   "Mix ceiling" = "11"), labels = ARM_LAB[ARMS_C], name = NULL) +
  scale_x_continuous(breaks = seq(50, 100, 10), labels = function(x) paste0(x, "%")) +
  labs(x = "Completeness threshold (% of 54 samples)", y = "Proteins at or above threshold",
       title = "A  Proteins quantified at a given completeness",
       subtitle = "Same 54 samples. Mix = 20% Krakatoa + 80% trypsin peptides; chosen channel = more complete protease per protein.") +
  th + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 2))
dd <- diffs[(Arm == "Mix +imp") | Arm == "Try1"]
dd <- rbind(dd[Arm == "Mix +imp"], copy(dd[Arm == "Try1"])[, Mix := "Mix30"], copy(dd[Arm == "Try1"])[, Mix := "Mix20"])
dd <- dd[tau >= 50]; dd[, `:=`(Arm = factor(Arm, levels = c("Try1","Mix +imp")), Mix = factor(Mix, levels = c("Mix30","Mix20")))]
p1b <- ggplot(dd, aes(tau, diff, colour = Arm, fill = Arm)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) + geom_line(linewidth = 0.8) +
  facet_wrap(~ Mix) +
  scale_colour_manual(values = COL[c("Try1","Mix +imp")], labels = c("Trypsin 1 - Trypsin 2 (prep-to-prep reference)",
                      "Mix deliverable - Trypsin 2"), name = NULL, aesthetics = c("colour","fill")) +
  scale_x_continuous(breaks = seq(50, 100, 10), labels = function(x) paste0(x, "%")) +
  labs(x = "Completeness threshold (% of 54 samples)", y = "Difference in proteins vs Trypsin 2",
       title = "B  Difference from Trypsin 2 (95% bootstrap CI over samples)") +
  th + theme(legend.position = "bottom")
f1 <- p1a / p1b + plot_layout(heights = c(1.15, 1))
ggsave(file.path(FIG, "fig1_threshold_curve.png"), f1, width = 9, height = 8.5, dpi = DPI, bg = "white")
ggsave(file.path(FIG, "fig1_threshold_curve.pdf"), f1, width = 9, height = 8.5, bg = "white")

# Fig 2: paired per-sample counts
ARMS_P <- c("Try2","Mix T-ch","Mix +imp")
pp <- rbind(ps[Arm == "Try2"][, .(SampleID, grp = "Trypsin 2", n, Arm)],
            ps[Mix == "Mix30" & Arm %in% c("Mix T-ch","Mix +imp")][, .(SampleID, grp = paste("Mix30:", fifelse(Arm == "Mix T-ch", "trypsin ch.", "deliverable")), n, Arm)],
            ps[Mix == "Mix20" & Arm == "Mix +imp"][, .(SampleID, grp = "Mix20: deliverable", n, Arm)])
pp[, grp := factor(grp, levels = c("Trypsin 2","Mix30: trypsin ch.","Mix30: deliverable","Mix20: deliverable"))]
f2 <- ggplot(pp, aes(grp, n)) + geom_line(aes(group = SampleID), colour = "grey80", linewidth = 0.3) +
  geom_point(aes(colour = Arm), size = 1.6, alpha = 0.85) +
  stat_summary(fun = median, geom = "crossbar", width = 0.45, linewidth = 0.35, colour = "black") +
  scale_colour_manual(values = COL[ARMS_P], guide = "none") +
  labs(x = NULL, y = "Proteins quantified per sample",
       title = "Proteins quantified per sample (lines join the same sample)",
       subtitle = sprintf("Median paired gain, Mix30 deliverable vs Trypsin 2: +%g (range +%d to +%d); all 54 samples higher.",
                          (ps_sum[Mix == "Mix30" & Arm == "Mix +imp", median_paired_diff_vs_Try2]),
                          as.integer(ps_sum[Mix == "Mix30" & Arm == "Mix +imp", min_diff]),
                          as.integer(ps_sum[Mix == "Mix30" & Arm == "Mix +imp", max_diff]))) + th
ggsave(file.path(FIG, "fig2_per_sample_paired.png"), f2, width = 6.5, height = 5, dpi = DPI, bg = "white")
ggsave(file.path(FIG, "fig2_per_sample_paired.pdf"), f2, width = 6.5, height = 5, bg = "white")

# Fig 3: stacked completeness bins in the style of Identification_Analysis.R
bins_of <- function(M, lab) { n <- rowSums(M); n <- n[n > 0]
  data.table(Arm = lab, Bin = bin_c(100 * n / ncol(M)))[, .N, by = .(Arm, Bin)] }
BL <- c("Trypsin 1", "Trypsin 2", "Trypsin 2\n(both libs)", "Mix30\ntrypsin ch.", "Mix30\nchosen ch.", "Mix30\n+ imputation",
        "Mix20\nchosen ch.", "Mix20\n+ imputation")
bn <- rbindlist(list(bins_of(R, BL[1]), bins_of(A, BL[2]), bins_of(A2, BL[3]), bins_of(MX$Mix30$B, BL[4]),
                     bins_of(MX$Mix30$C, BL[5]), bins_of(MX$Mix30$D, BL[6]), bins_of(MX$Mix20$C, BL[7]),
                     bins_of(MX$Mix20$D, BL[8])))
bn[, `:=`(Arm = factor(Arm, levels = BL), Bin = factor(Bin, levels = COMP_LEVELS))]
fwrite(bn, file.path(TBL, "completeness_bins.csv"))
f3 <- ggplot(bn, aes(Arm, N, fill = Bin)) +
  geom_col(width = 0.78, colour = "white", linewidth = 0.4, position = position_stack(reverse = FALSE)) +
  geom_text(aes(label = fifelse(N >= 40, as.character(N), "")), position = position_stack(vjust = 0.5),
            size = 2.5, colour = "grey10") +
  scale_fill_manual(values = COMP_COLS, breaks = rev(COMP_LEVELS), name = "Completeness\n(% of 54 samples)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03)), labels = comma) +
  labs(x = NULL, y = "Protein groups", title = "Per-protein completeness distribution",
       subtitle = "Counts >= 40 labelled. Green = quantified in every sample.") + th
ggsave(file.path(FIG, "fig3_completeness_bins.png"), f3, width = 9, height = 5.5, dpi = DPI, bg = "white")
ggsave(file.path(FIG, "fig3_completeness_bins.pdf"), f3, width = 9, height = 5.5, bg = "white")

# Fig 4: where Try2 proteins lose completeness — by abundance, vs prep reference
st <- strat[comparison %in% c("Mix30 +imp vs Try2", "Try1 vs Try2 (reference)")]
st <- melt(st, id.vars = c("comparison","quintile"), measure.vars = c("frac_less_complete","frac_more_complete"),
           variable.name = "dir", value.name = "frac")
st[, dir := factor(dir, labels = c("less complete than Trypsin 2", "more complete than Trypsin 2"))]
st[, comparison := factor(comparison, levels = c("Mix30 +imp vs Try2","Try1 vs Try2 (reference)"),
                          labels = c("Mix30 deliverable", "Trypsin 1 (reference)"))]
f4 <- ggplot(st, aes(quintile, frac, fill = comparison)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72, colour = "white", linewidth = 0.4) +
  facet_wrap(~ dir) + scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("Mix30 deliverable" = COL[["Mix +imp"]], "Trypsin 1 (reference)" = COL[["Try1"]]), name = NULL) +
  labs(x = "Trypsin 2 abundance quintile (Q1 = lowest)", y = "Share of Trypsin 2 proteins",
       title = "Which Trypsin 2 proteins change completeness",
       subtitle = "Proteins quantified in both. The Trypsin 1 bars show normal prep-to-prep change.") +
  th + theme(legend.position = "bottom")
ggsave(file.path(FIG, "fig4_loss_by_abundance.png"), f4, width = 8, height = 4.8, dpi = DPI, bg = "white")
ggsave(file.path(FIG, "fig4_loss_by_abundance.pdf"), f4, width = 8, height = 4.8, bg = "white")

cat("\nDone. Tables ->", TBL, "\nFigures ->", FIG, "\n")
