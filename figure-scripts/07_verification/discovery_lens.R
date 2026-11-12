# =============================================================================
# Discovery lens: is the Mix30 COMBINED matrix (QuantAdvisor deliverable) better
# for discovery than trypsin-only Trypsin 2?   2026-09-28
#
# Discovery question = how many proteins are analysis-ready under the usual
# presence filters (50-80% of samples), and are the extra proteins usable?
#   1. Discovery-ready proteins at 50/60/70/80% filters, with paired bootstrap CI
#   2. Gain/loss waterfall: retained / lost / gained (shared) / gained (new)
#   3. Evidence depth of each class (median distinct peptides per run)
#   4. Robust discovery set = passes filter AND >= 2 peptides per run (median)
#   5. Quantifiability: inter-day CV of the same proteins in the reproducibility cohort
# Combined = one channel per protein (more complete protease) + strict
# cross-protease imputation, exactly as SI Fig S3 describes.
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales); library(duckdb)
})
set.seed(20260928)

.sd <- function() { a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f)) normalizePath(dirname(f)) else getwd() }
source(file.path(.sd(), "..", "config.R"))   # DATA_DIR, OUT_ROOT, step_dirs()
source(file.path(.sd(), "helpers.R"))
D <- step_dirs("07_verification")
FIG <- file.path(D$figures, "discovery"); TBL <- file.path(D$tables, "discovery")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE); dir.create(TBL, showWarnings = FALSE, recursive = TRUE)

FILTERS <- c(0.50, 0.60, 0.70, 0.80)
MIN_PEP <- 2                              # robust = median distinct peptides per detected run >= 2
C_TRY <- "#2C6FB3"; C_KRA <- "#E08214"; C_CMB <- "#1B7837"; C_LOSS <- "#B2474A"; C_GREY <- "#9A9A9A"

# ====================================================================== data ==
cat("== Discovery lens: Mix30 combined vs Trypsin 2 ==\n")
assert_panel_a()
TR <- load_disc(PG_TRY); KR <- load_disc(PG_KRA)
S <- sort(unique(TR[Cond == "Try2", SampleID]))
stopifnot(length(S) == 54L, setequal(S, unique(TR[Cond == "Mix30", SampleID])))
K <- length(S)
A  <- pres_mat(TR[Cond == "Try2"], S)
MX <- build_mix_arms(TR, KR, "Mix30", S)
D  <- MX$D; chosen <- MX$chosen
# guard: must match the main analysis
stopifnot(sum(rowSums(A) >= 27) == 1053L, sum(rowSums(D) >= 27) == 1278L, sum(rowSums(D) == 54) == 602L)
cat("Arms rebuilt; main-analysis counts reproduced (1053 / 1278 / 602).\n")

# ---- evidence depth from report.parquet (1% run-level precursor q and PG q) ----
pep_depth <- function(parq, where) {
  con <- dbConnect(duckdb()); on.exit(dbDisconnect(con, shutdown = TRUE))
  as.data.table(dbGetQuery(con, sprintf("
    SELECT \"Protein.Group\" AS FeatureID, Run, count(DISTINCT \"Stripped.Sequence\") AS n_pep
    FROM read_parquet('%s')
    WHERE (%s) AND Decoy = 0 AND \"Q.Value\" <= 0.01 AND \"PG.Q.Value\" <= 0.01
    GROUP BY 1, 2", parq, where)))
}
PQ_T <- file.path(SEARCH, "Trypsin_Library_Search", "report.parquet")
PQ_K <- file.path(SEARCH, "Krakatoa_Library_Search", "report.parquet")
W_TRY2  <- "Run LIKE '251117_DIA_Plasma_Try%'"
W_MIX30 <- "Run LIKE '251117_DIA_Plasma_Mix%' AND Run NOT LIKE '%20min%'"
pep_runs <- function(parq, where) {
  con <- dbConnect(duckdb()); on.exit(dbDisconnect(con, shutdown = TRUE))
  as.data.table(dbGetQuery(con, sprintf("
    SELECT \"Protein.Group\" AS FeatureID, \"Stripped.Sequence\" AS pep, count(DISTINCT Run) AS runs
    FROM read_parquet('%s')
    WHERE (%s) AND Decoy = 0 AND \"Q.Value\" <= 0.01 AND \"PG.Q.Value\" <= 0.01
    GROUP BY 1, 2", parq, where)))
}
# per protein: median distinct peptides per detected run; distinct peptides across the cohort;
# consistency = share of detected runs that contain the protein's most frequent peptide
depth <- function(parq, where) {
  a <- pep_depth(parq, where)[, .(med_pep = median(n_pep), det_runs = .N), by = FeatureID]
  b <- pep_runs(parq, where)[, .(tot_pep = .N, top_runs = max(runs)), by = FeatureID]
  merge(a, b, by = "FeatureID")[, consistency := top_runs / det_runs][]
}
dp_try2 <- depth(PQ_T, W_TRY2)
dp_mixT <- depth(PQ_T, W_MIX30)
dp_mixK <- depth(PQ_K, W_MIX30)

# ---- per-protein discovery table ----
prots <- union(rownames(A), rownames(D))
nA <- setNames(rep(0L, length(prots)), prots); nD <- nA
nA[rownames(A)] <- as.integer(rowSums(A)); nD[rownames(D)] <- as.integer(rowSums(D))
P <- data.table(FeatureID = prots, n_try2 = nA, n_mix = nD, channel = unname(chosen[prots]))
P <- merge(P, dp_try2[, .(FeatureID, pep_try2 = med_pep, tot_try2 = tot_pep, cons_try2 = consistency)], by = "FeatureID", all.x = TRUE)
P <- merge(P, dp_mixT[, .(FeatureID, pep_mixT = med_pep, tot_mixT = tot_pep, cons_mixT = consistency)], by = "FeatureID", all.x = TRUE)
P <- merge(P, dp_mixK[, .(FeatureID, pep_mixK = med_pep, tot_mixK = tot_pep, cons_mixK = consistency)], by = "FeatureID", all.x = TRUE)
P[, `:=`(pep_mix  = fifelse(channel == "Krakatoa", pep_mixK, pep_mixT),
         tot_mix  = fifelse(channel == "Krakatoa", tot_mixK, tot_mixT),
         cons_mix = fifelse(channel == "Krakatoa", cons_mixK, cons_mixT))]
# evidence tier in the matrix that carries the protein
tier <- function(med, tot) fcase(med >= MIN_PEP, ">= 2 peptides per run", tot >= MIN_PEP, ">= 2 peptides in cohort",
                                 !is.na(tot), "1 peptide only", default = NA_character_)
P[, `:=`(tier_try2 = tier(pep_try2, tot_try2), tier_mix = tier(pep_mix, tot_mix))]
genes <- unique(rbind(fread(PG_TRY, select = c("Protein.Group","Genes","First.Protein.Description")),
                      fread(PG_KRA, select = c("Protein.Group","Genes","First.Protein.Description"))), by = "Protein.Group")
P <- merge(P, genes[, .(FeatureID = Protein.Group, Genes, Description = First.Protein.Description)], by = "FeatureID", all.x = TRUE)

classify <- function(n0, n1, kk) fcase(n0 >= kk & n1 >= kk, "retained", n0 >= kk & n1 < kk, "lost",
                                       n0 > 0 & n0 < kk & n1 >= kk, "gained: shared",
                                       n0 == 0 & n1 >= kk, "gained: new", default = NA_character_)
for (f in FILTERS) P[, (paste0("class_", 100 * f)) := classify(n_try2, n_mix, ceiling(f * K))]
fwrite(P[order(-n_mix)], file.path(TBL, "discovery_per_protein.csv"))

# ======================================================= 1 + 4. filter summary ==
robust_rows <- function(M, pep) { r <- rownames(M); keep <- r[r %in% pep[!is.na(v) & v >= MIN_PEP, FeatureID]]
  M[keep, , drop = FALSE] }
A_rob <- robust_rows(A, P[, .(FeatureID, v = pep_try2)])          # strict: >= 2 peptides per run
D_rob <- robust_rows(D, P[, .(FeatureID, v = pep_mix)])
A_len <- robust_rows(A, P[, .(FeatureID, v = tot_try2)])          # lenient: >= 2 peptides across the cohort
D_len <- robust_rows(D, P[, .(FeatureID, v = tot_mix)])
boot_gain <- function(M1, M0, kk, B = 2000L) {
  g <- replicate(B, { i <- sample.int(K, K, TRUE)
    sum(rowSums(M1[, i, drop = FALSE]) >= kk) - sum(rowSums(M0[, i, drop = FALSE]) >= kk) })
  quantile(g, c(0.025, 0.975))
}
fs <- rbindlist(lapply(FILTERS, function(f) { kk <- ceiling(f * K)
  ci <- boot_gain(D, A, kk); cr <- boot_gain(D_rob, A_rob, kk); cl <- boot_gain(D_len, A_len, kk)
  data.table(filter_pct = 100 * f, min_samples = kk,
             try2 = sum(rowSums(A) >= kk), mix30 = sum(rowSums(D) >= kk),
             gain = sum(rowSums(D) >= kk) - sum(rowSums(A) >= kk), gain_lo = ci[[1]], gain_hi = ci[[2]],
             try2_robust = sum(rowSums(A_rob) >= kk), mix30_robust = sum(rowSums(D_rob) >= kk),
             gain_robust = sum(rowSums(D_rob) >= kk) - sum(rowSums(A_rob) >= kk),
             gain_robust_lo = cr[[1]], gain_robust_hi = cr[[2]],
             try2_lenient = sum(rowSums(A_len) >= kk), mix30_lenient = sum(rowSums(D_len) >= kk),
             gain_lenient = sum(rowSums(D_len) >= kk) - sum(rowSums(A_len) >= kk),
             gain_lenient_lo = cl[[1]], gain_lenient_hi = cl[[2]]) }))
fwrite(fs, file.path(TBL, "discovery_filter_summary.csv"))
cat("\n--- discovery-ready proteins (all / robust = >=2 peptides per run) ---\n"); print(fs)

# ============================================================ 2. waterfall ==
wf <- rbindlist(lapply(FILTERS, function(f) { cl <- P[[paste0("class_", 100 * f)]]
  data.table(filter_pct = 100 * f, try2 = sum(cl %in% c("retained","lost")), retained = sum(cl == "retained", na.rm = TRUE),
             lost = sum(cl == "lost", na.rm = TRUE), gained_shared = sum(cl == "gained: shared", na.rm = TRUE),
             gained_new_krakatoa = sum(cl == "gained: new" & P$channel == "Krakatoa", na.rm = TRUE),
             gained_new_trypsin = sum(cl == "gained: new" & P$channel == "Trypsin", na.rm = TRUE)) }))
wf[, mix30 := retained + gained_shared + gained_new_krakatoa + gained_new_trypsin]
stopifnot(all(wf$mix30 == fs$mix30), all(wf$try2 == fs$try2))
fwrite(wf, file.path(TBL, "discovery_waterfall.csv"))
cat("\n--- waterfall ---\n"); print(wf)

# ======================================================= 3. class evidence ==
ev <- rbindlist(lapply(FILTERS, function(f) {
  x <- P[!is.na(get(paste0("class_", 100 * f)))][, cls := get(paste0("class_", 100 * f))]
  x[, `:=`(pep = fifelse(cls == "lost", pep_try2, pep_mix), tot = fifelse(cls == "lost", tot_try2, tot_mix),
           cons = fifelse(cls == "lost", cons_try2, cons_mix))]     # depth in the matrix the protein comes from
  x[, .(proteins = .N, median_peptides = as.numeric(median(pep, na.rm = TRUE)),
        frac_single_peptide_per_run = round(mean(pep < MIN_PEP, na.rm = TRUE), 3),
        frac_one_peptide_in_cohort = round(mean(tot < MIN_PEP, na.rm = TRUE), 3),
        median_peptides_in_cohort = as.numeric(median(tot, na.rm = TRUE)),
        median_top_peptide_consistency = round(median(cons, na.rm = TRUE), 3),
        frac_krakatoa_channel = round(mean(channel == "Krakatoa"), 3),
        median_completeness_try2 = round(100 * median(n_try2) / K, 1),
        median_completeness_mix = round(100 * median(n_mix) / K, 1)), by = cls][, filter_pct := 100 * f][] }))
fwrite(ev, file.path(TBL, "discovery_class_evidence.csv"))
cat("\n--- evidence depth by class ---\n"); print(ev)
fwrite(P[class_50 %in% c("gained: new","gained: shared")][order(class_50, -n_mix),
         .(FeatureID, Genes, Description, class_50, channel, n_try2, n_mix, pep_mix, tot_mix, cons_mix, tier_mix)],
       file.path(TBL, "gained_proteins_at_50pct.csv"))
fwrite(P[class_50 == "lost"][order(n_mix), .(FeatureID, Genes, Description, n_try2, n_mix, pep_try2, tot_try2, cons_try2, pep_mixT)],
       file.path(TBL, "lost_proteins_at_50pct.csv"))

# ============================================== 5. reproducibility-cohort CV ==
# The published definition (Reproducibility_Data_Processing.R): per sample, mean intensity per day (linear), inter-day CV; well-quantified =
# measured on all days and CV <= 20%. Channel matched to the discovery channel.
repro_cv <- function(path) {
  L <- melt_pg(fread(path, sep = "\t", header = TRUE, showProgress = FALSE))
  b <- sub("\\.raw$", "", basename(gsub("\\\\", "/", L$Run)), ignore.case = TRUE)
  p <- tstrsplit(b, "_", fixed = TRUE); n <- length(p)
  L[, `:=`(Day = as.integer(p[[n - 2]]), Sample = p[[n - 1]])]
  ed <- L[, .(ED = uniqueN(Day)), by = Sample]
  dm <- L[is.finite(Intensity), .(DayMean = mean(Intensity)), by = .(FeatureID, Sample, Day)]
  cv <- dm[, .(CV = if (.N < 2) NA_real_ else 100 * sd(DayMean) / mean(DayMean), MD = .N), by = .(FeatureID, Sample)]
  cv <- merge(cv, ed, by = "Sample")[, wq := MD == ED & is.finite(CV) & CV <= 20]
  cv[, .(repro_samples = .N, median_interday_cv = median(CV[MD == ED], na.rm = TRUE),
         frac_samples_wq = mean(wq)), by = FeatureID]
}
rc <- rbind(repro_cv(RP_TRY)[, channel := "Trypsin"], repro_cv(RP_HTA)[, channel := "Krakatoa"])
qa <- merge(P[!is.na(class_50), .(FeatureID, class_50, channel)], rc, by = c("FeatureID","channel"), all.x = TRUE)
qa_sum <- qa[, .(proteins = .N, in_repro = sum(!is.na(repro_samples)),
                 frac_in_repro = round(mean(!is.na(repro_samples)), 3),
                 median_interday_cv = round(median(median_interday_cv, na.rm = TRUE), 1),
                 frac_wq_in_any_sample = round(mean(frac_samples_wq > 0, na.rm = TRUE), 3)), by = .(class_50)]
fwrite(qa, file.path(TBL, "repro_quantifiability_per_protein.csv")); fwrite(qa_sum, file.path(TBL, "repro_quantifiability_summary.csv"))
cat("\n--- reproducibility-cohort check (channel-matched; class at 50% filter) ---\n"); print(qa_sum)

# ============================================ 6. plasma plausibility (HPA) ==
# Human Protein Atlas table (proteinatlas.tsv, downloaded 2026-09-28, HPA release dated 2026-06-03). Not
# redistributed here (CC BY-SA 4.0): download https://www.proteinatlas.org/download/proteinatlas.tsv.zip,
# unzip, and set HPA_TSV to the file (default: <TWO_PROTEASE_DATA>/reference/proteinatlas.tsv).
#   ms_blood       = has an MS-based blood concentration (HPA plasma proteome; PeptideAtlas, i.e. TRYPTIC
#                    evidence -> biased against proteins only a non-tryptic protease can see)
#   secreted_blood = HPA secretome location "Secreted to blood" (prediction/annotation; digestion-independent)
#   any_blood      = ms_blood OR immunoassay blood concentration OR secreted_blood
HPA_TSV <- Sys.getenv("HPA_TSV", file.path(DATA_DIR, "reference", "proteinatlas.tsv"))
if (!file.exists(HPA_TSV)) stop("Human Protein Atlas table not found: ", HPA_TSV,
  "\n  Download https://www.proteinatlas.org/download/proteinatlas.tsv.zip, unzip, and set HPA_TSV.")
hpa <- fread(HPA_TSV, sep = "\t", quote = "\"",
             select = c("Gene", "Uniprot", "Secretome location", "Blood concentration - Conc. blood IM [pg/L]",
                        "Blood concentration - Conc. blood MS [pg/L]"), colClasses = "character")
setnames(hpa, c("Gene", "Uniprot", "secretome", "conc_im", "conc_ms"))
hpa <- hpa[Uniprot != "", .(acc = trimws(unlist(strsplit(Uniprot, ",")))), by = .(Gene, secretome, conc_im, conc_ms)]
hpa <- hpa[, .(ms_blood = any(!is.na(conc_ms) & conc_ms != ""), im_blood = any(!is.na(conc_im) & conc_im != ""),
               secreted_blood = any(grepl("Secreted to blood", secretome))), by = acc]
hpa[, any_blood := ms_blood | im_blood | secreted_blood]
P[, acc := sub(";.*", "", FeatureID)]
P <- merge(P, hpa, by = "acc", all.x = TRUE)
P[, in_hpa := !is.na(ms_blood)]
for (v in c("ms_blood","im_blood","secreted_blood","any_blood")) P[is.na(get(v)), (v) := FALSE]
cat(sprintf("\nHPA mapping: %d of %d proteins matched by UniProt accession\n", sum(P$in_hpa), nrow(P)))
pl_cls <- P[!is.na(class_50), .(proteins = .N, in_hpa = sum(in_hpa), frac_ms_blood = round(mean(ms_blood), 3),
                                frac_secreted_blood = round(mean(secreted_blood), 3), frac_any_blood = round(mean(any_blood), 3)),
            by = class_50][order(match(class_50, c("retained","lost","gained: shared","gained: new")))]
pl_tier <- P[class_50 == "gained: new", .(proteins = .N, frac_ms_blood = round(mean(ms_blood), 3),
                                          frac_any_blood = round(mean(any_blood), 3)), by = tier_mix]
pl_filt <- rbindlist(lapply(FILTERS, function(f) { kk <- ceiling(f * K)
  data.table(filter_pct = 100 * f,
             try2_ms_blood = P[n_try2 >= kk, sum(ms_blood)], mix30_ms_blood = P[n_mix >= kk, sum(ms_blood)],
             try2_any_blood = P[n_try2 >= kk, sum(any_blood)], mix30_any_blood = P[n_mix >= kk, sum(any_blood)]) }))
pl_filt[, `:=`(gain_ms_blood = mix30_ms_blood - try2_ms_blood, gain_any_blood = mix30_any_blood - try2_any_blood)]
fwrite(pl_cls, file.path(TBL, "plasma_plausibility_by_class.csv"))
fwrite(pl_tier, file.path(TBL, "plasma_plausibility_gained_new_by_evidence.csv"))
fwrite(pl_filt, file.path(TBL, "plasma_plausibility_by_filter.csv"))
fwrite(P[order(-n_mix)], file.path(TBL, "discovery_per_protein.csv"))    # re-write with HPA columns
cat("\n--- plasma plausibility by class (50% filter) ---\n"); print(pl_cls)
cat("\n--- gained-new proteins: plausibility by evidence tier ---\n"); print(pl_tier)
cat("\n--- plasma-supported discovery-ready proteins ---\n"); print(pl_filt)

# ================================================================= figures ==
th <- theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold", size = 12),
  plot.subtitle = element_text(size = 8.8, colour = "grey30"), axis.text = element_text(colour = "black"),
  panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
  strip.background = element_blank(), strip.text = element_text(face = "bold"), legend.position = "bottom")
save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 600, bg = "white")
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
}

# Fig D1: discovery-ready proteins at common filters, split into three evidence tiers
TIERS <- c(">= 2 peptides per run", ">= 2 peptides in cohort", "1 peptide only")
d1 <- rbindlist(lapply(FILTERS, function(f) { kk <- ceiling(f * K)
  rbind(P[n_try2 >= kk, .N, by = .(tier = tier_try2)][, arm := "Trypsin 2"],
        P[n_mix  >= kk, .N, by = .(tier = tier_mix)][, arm := "Mix30 combined"])[, filter_pct := 100 * f] }))
d1[is.na(tier), tier := "1 peptide only"]          # no parquet evidence at 1% run q (PG passes only globally)
d1 <- d1[, .(n = sum(N)), by = .(filter_pct, arm, tier)]
d1[, `:=`(arm = factor(arm, levels = c("Trypsin 2","Mix30 combined")), tier = factor(tier, levels = rev(TIERS)),
          fx = factor(paste0(">= ", filter_pct, "% of samples"), levels = paste0(">= ", 100 * FILTERS, "% of samples")))]
d1[, fillkey := factor(paste(arm, tier, sep = " | "),          # first level stacks on top -> strongest evidence at the bottom
      levels = as.vector(t(outer(c("Trypsin 2","Mix30 combined"), rev(TIERS), paste, sep = " | "))))]
fwrite(d1, file.path(TBL, "discovery_evidence_tiers.csv"))
tot <- d1[, .(n = sum(n)), by = .(fx, arm, filter_pct)]
lab <- merge(tot[arm == "Mix30 combined"], fs[, .(filter_pct, gain, gain_lenient, gain_robust)], by = "filter_pct")
FILLS <- c("Trypsin 2 | >= 2 peptides per run" = "#1F4E85", "Trypsin 2 | >= 2 peptides in cohort" = "#5C8FCB",
           "Trypsin 2 | 1 peptide only" = "#BCD3EC",
           "Mix30 combined | >= 2 peptides per run" = "#135A27", "Mix30 combined | >= 2 peptides in cohort" = "#3E9A55",
           "Mix30 combined | 1 peptide only" = "#B9DFC0")
p_d1 <- ggplot(d1, aes(arm, n, fill = fillkey)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.5) +
  geom_text(data = tot, aes(arm, n, label = comma(n)), inherit.aes = FALSE, vjust = -0.4, size = 3.2) +
  geom_text(data = lab, aes(arm, n, label = sprintf("%+d total\n%+d (>=2 in cohort)", gain, gain_lenient)),
            inherit.aes = FALSE, vjust = -1.0, size = 2.7, fontface = "bold", colour = C_CMB, lineheight = 0.9) +
  facet_wrap(~ fx, nrow = 1) +
  scale_fill_manual(values = FILLS, breaks = names(FILLS), labels = sub(" \\| ", ": ", names(FILLS)), name = NULL) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.24))) +
  labs(x = NULL, y = "Proteins passing the filter",
       title = "Discovery-ready proteins: Mix30 combined vs trypsin-only",
       subtitle = paste0("Same 54 plasma samples. Combined = more complete protease per protein + cross-protease imputation.\n",
                         "Shading = peptide evidence in the channel that carries the protein (dark = strongest). Green: gain vs Trypsin 2.")) +
  th + theme(axis.text.x = element_text(size = 8.5)) + guides(fill = guide_legend(nrow = 2, byrow = TRUE, reverse = FALSE))
save_fig(p_d1, "figD1_discovery_ready_proteins", 10, 5.6)

# Fig D2: waterfall at 50% and 70%
stopifnot(all(wf$gained_new_trypsin == 0))      # every new protein comes through the Krakatoa channel
mk_wf <- function(f) { w <- wf[filter_pct == 100 * f]
  steps <- data.table(step = c("Trypsin 2", "lost", "gained:\nshared", "gained: new\n(Krakatoa ch.)", "Mix30\ncombined"),
                      delta = c(w$try2, -w$lost, w$gained_shared, w$gained_new_krakatoa, w$mix30),
                      kind = c("total","loss","gain_t","gain_k","total_mix"))
  steps[, end := c(delta[1], delta[1] + cumsum(delta[2:4]), delta[5])]
  steps[, start := c(0, head(end, -1)[1:3], 0)]
  steps[, `:=`(i = .I, filter = paste0(">= ", 100 * f, "% of samples"))][] }
wfd <- rbind(mk_wf(0.5), mk_wf(0.7))
wfd[, step := factor(step, levels = unique(step))]
wfd[, lab := fifelse(kind %in% c("total","total_mix"), comma(delta), sprintf("%+d", delta))]
p_d2 <- ggplot(wfd) +
  geom_rect(aes(xmin = i - 0.38, xmax = i + 0.38, ymin = start, ymax = end, fill = kind), colour = "white", linewidth = 0.4) +
  geom_segment(data = wfd[i < 5], aes(x = i + 0.38, xend = i + 0.62, y = end, yend = end), colour = "grey55", linewidth = 0.3) +
  geom_text(aes(i, pmax(start, end), label = lab), vjust = -0.45, size = 3.1) +
  facet_wrap(~ filter, nrow = 1) +
  scale_x_continuous(breaks = 1:5, labels = levels(wfd$step)) +
  scale_fill_manual(values = c(total = C_TRY, loss = C_LOSS, gain_t = "#7FA9D6", gain_k = C_KRA, total_mix = C_CMB), guide = "none") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
  coord_cartesian(ylim = c(800, NA)) +
  labs(x = NULL, y = "Proteins passing the filter", title = "Where the discovery gain comes from",
       subtitle = paste0("Lost = passes in Trypsin 2 but not in Mix30; gained: shared = trypsin sees it but below the filter;\n",
                         "gained: new = not quantified by Trypsin 2 at all. y-axis starts at 800.")) +
  th + theme(axis.text.x = element_text(size = 8))
save_fig(p_d2, "figD2_gain_loss_waterfall", 10, 5)

# Fig D3: evidence tier and peptide consistency by class (50% filter)
e50 <- P[!is.na(class_50)][, `:=`(tier = fifelse(class_50 == "lost", tier_try2, tier_mix),
                                  cons = fifelse(class_50 == "lost", cons_try2, cons_mix))]
e50[is.na(tier), tier := "1 peptide only"]
CLS <- c("retained","lost","gained: shared","gained: new")
e50[, `:=`(class_50 = factor(class_50, levels = CLS), tier = factor(tier, levels = rev(TIERS)))]
cnt <- e50[, .N, by = class_50]
p3a <- ggplot(e50, aes(class_50, fill = tier)) +
  geom_bar(position = "fill", width = 0.7, colour = "white", linewidth = 0.4) +
  geom_text(data = cnt, aes(class_50, 1.04, label = paste0("n = ", N)), inherit.aes = FALSE, size = 3) +
  scale_fill_manual(values = c(">= 2 peptides per run" = "#3F007D", ">= 2 peptides in cohort" = "#807DBA",
                               "1 peptide only" = "#DADAEB"), breaks = TIERS, name = "Peptide evidence") +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.07))) +
  labs(x = "Protein class at the 50% filter", y = "Share of proteins") + th + theme(legend.position = "right")
p_d3 <- p3a + labs(title = NULL) + plot_annotation(
  title = "Evidence behind the proteins gained and lost at the 50% filter",
  subtitle = "Distinct peptides, measured in the matrix that carries the protein\n(lost: Trypsin 2; others: Mix30 chosen channel).",
  theme = theme(plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 8.8, colour = "grey30")))
save_fig(p_d3, "figD3_evidence_depth_by_class", 6.5, 5.6)

# Fig D4: reproducibility-cohort precision of the same proteins (channel-matched)
q4 <- qa[!is.na(median_interday_cv)][, class_50 := factor(class_50, levels = c("retained","lost","gained: shared","gained: new"))]
q4n <- merge(qa_sum, q4[, .(ymax = max(median_interday_cv)), by = class_50], by = "class_50", all.x = TRUE)
q4n[, class_50 := factor(class_50, levels = levels(q4$class_50))]
p_d4 <- ggplot(q4, aes(class_50, median_interday_cv, colour = channel)) +
  geom_hline(yintercept = 20, linetype = "22", colour = "grey50") +
  geom_boxplot(aes(group = class_50), outlier.shape = NA, width = 0.55, colour = "grey30", fill = NA) +
  geom_jitter(width = 0.15, height = 0, size = 1.3, alpha = 0.75) +
  geom_text(data = q4n, aes(class_50, -Inf, label = sprintf("%d of %d\nin repro", in_repro, proteins)), inherit.aes = FALSE,
            vjust = -0.3, size = 2.8, lineheight = 0.9) +
  scale_colour_manual(values = c(Trypsin = C_TRY, Krakatoa = C_KRA), name = "Channel") +
  scale_y_continuous(limits = c(-12, NA), breaks = seq(0, 200, 20)) +
  labs(x = "Protein class at the 50% discovery filter", y = "Median inter-day CV (%)",
       title = "Are the gained proteins quantifiable? Reproducibility-cohort check",
       subtitle = paste0("Same proteins, same protease channel, in the separate 8-sample x 3-day x 8-replicate cohort (Exploris 480, 15 min).\n",
                         "Only proteins that the shallower reproducibility run detects can be assessed. Dashed line = 20% CV.")) + th
save_fig(p_d4, "figD4_repro_quantifiability", 7.5, 5.2)

# Fig D5: plasma plausibility of each class (HPA)
d5 <- melt(pl_cls, id.vars = c("class_50","proteins"), measure.vars = c("frac_ms_blood","frac_secreted_blood","frac_any_blood"),
           variable.name = "evidence", value.name = "frac")
d5[, `:=`(class_50 = factor(class_50, levels = c("retained","lost","gained: shared","gained: new")),
          evidence = factor(evidence, levels = c("frac_ms_blood","frac_secreted_blood","frac_any_blood"),
                            labels = c("MS-detected in blood (PeptideAtlas; tryptic)", "Predicted secreted to blood",
                                       "Any HPA blood evidence")))]
p_d5 <- ggplot(d5, aes(class_50, frac, fill = evidence)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.74, colour = "white", linewidth = 0.4) +
  geom_text(aes(label = percent(frac, accuracy = 1)), position = position_dodge(width = 0.8), vjust = -0.4, size = 2.7) +
  geom_text(data = unique(d5[, .(class_50, proteins)]), aes(class_50, -0.04, label = paste0("n = ", proteins)),
            inherit.aes = FALSE, size = 2.9, colour = "grey30") +
  scale_fill_manual(values = c("#4C4C4C", "#9E9AC8", "#C7C7C7"), name = NULL) +
  scale_y_continuous(labels = percent, limits = c(-0.07, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Protein class at the 50% discovery filter", y = "Share of proteins",
       title = "Are the proteins plausible plasma proteins? (Human Protein Atlas)",
       subtitle = paste0("PeptideAtlas plasma evidence is built from tryptic digests, so it under-counts proteins\n",
                         "that only a non-tryptic protease can detect; the secretion annotation does not depend on digestion.")) +
  th + guides(fill = guide_legend(nrow = 1))
save_fig(p_d5, "figD5_plasma_plausibility", 8.5, 5.2)

# Summary composite (main-text candidate): D1 over D2
save_fig(p_d1 / p_d2 + plot_layout(heights = c(1.1, 1)) + plot_annotation(tag_levels = "A"),
         "figD0_discovery_summary", 10, 10.5)
cat("\nDone. Tables ->", TBL, "\nFigures ->", FIG, "\n")
