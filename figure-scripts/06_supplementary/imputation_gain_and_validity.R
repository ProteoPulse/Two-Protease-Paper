# =============================================================================
# imputation_gain_and_validity.R — one representative figure: what cross-protease
# imputation BUYS (given the app commits to one channel per protein), and that the
# filled values are VALID.
#
# The advisor delivers one value per protein per sample from a single CHOSEN
# channel. Imputation fills that channel's gaps from the other protease. So the
# gain that matters is per-channel completeness, not the union (which two proteases
# already make near-complete). This figure puts the two questions side by side:
#   (A) GAIN     — well-quantified counts per channel + union, measured vs filled.
#   (B) VALIDITY — masking recovery error under MNAR: the fill beats every standard
#                  imputer, so the recovered values can be trusted.
#
# Reuses the tables already written by discovery_imputation.R and masking_validation.R;
# recomputes nothing. Usage: Rscript figure-scripts/06_supplementary/imputation_gain_and_validity.R
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

.sd <- function() { a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f)) normalizePath(dirname(f)) else getwd() }
source(file.path(.sd(), "..", "config.R"))   # DATA_DIR, OUT_ROOT, step_dirs()
D <- step_dirs("06_supplementary"); FIG <- D$figures; TBL <- D$tables

TRY_COL <- "#3C7DA6"; HTA_COL <- "#E1812C"; CMB_COL <- "#4C4C4C"
PALE <- function(x) grDevices::adjustcolor(x, alpha.f = 0.42)
TAU <- 0.50
CONDS <- c("Mix30", "Mix20")

# ---------------------------------------------------------------- (A) gain ---
cc <- fread(file.path(TBL, "discovery_completeness_counts.csv"))[tau == TAU]
gain <- rbindlist(lapply(CONDS, function(cd) {
  r <- cc[Condition == cd]
  data.table(Condition = cd,
    channel = factor(c("Krakatoa","Trypsin","Combined\n(union)"),
                     levels = c("Krakatoa","Trypsin","Combined\n(union)")),
    before = c(r$Krakatoa_before, r$Trypsin_before, r$Combined_before),
    after  = c(r$Krakatoa_strict, r$Trypsin_strict, r$Combined_strict))
}))
gain[, delta := after - before]
gl <- melt(gain, id.vars = c("Condition","channel","delta"),
           measure.vars = c("before","after"), variable.name = "phase", value.name = "N")
gl[, phase := factor(phase, levels = c("before","after"), labels = c("Measured","+ imputation"))]
gl[, Condition := factor(Condition, levels = CONDS)]
gl[, fillkey := paste(channel, phase)]
chan_col <- c("Krakatoa" = HTA_COL, "Trypsin" = TRY_COL, "Combined\n(union)" = CMB_COL)
fillvals <- unlist(lapply(names(chan_col), function(ch)
  setNames(c(PALE(chan_col[[ch]]), chan_col[[ch]]),
           c(paste(ch, "Measured"), paste(ch, "+ imputation")))))
dlab <- gain[, .(Condition = factor(Condition, levels = CONDS), channel,
                 N = after, lab = sprintf("+%d", delta))]

pA <- ggplot(gl, aes(channel, N, fill = fillkey, group = phase)) +
  geom_col(position = position_dodge(0.72), width = 0.66, color = "grey25", linewidth = 0.2) +
  geom_text(data = dlab, aes(x = channel, y = N, label = lab), inherit.aes = FALSE,
            vjust = -0.5, size = 3.2, fontface = "bold", color = "grey20") +
  facet_wrap(~Condition, nrow = 1) +
  scale_fill_manual(values = fillvals, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Well-quantified proteins (>=50% of runs)",
       title = "A  Gain from imputation in the chosen channel",
       subtitle = "Pale = measured, solid = after cross-protease fill. The sparse channel gains; the union\nis already near-complete from having two proteases, so imputation barely moves it.") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8.4, color = "grey30"),
        axis.text = element_text(color = "black"),
        strip.background = element_rect(fill = "grey92", color = NA), strip.text = element_text(face = "bold"))

# ------------------------------------------------------------ (B) validity ---
ma <- fread(file.path(TBL, "masking_aggregate_discovery.csv"))
val <- ma[regime == "mnar" & Recipient == "Krakatoa" &
          method %in% c("cross_protease","protein_mean","protein_knn","sample_median","minprob","mindet")]
lab_map <- c(cross_protease = "cross-protease", protein_mean = "protein mean",
             protein_knn = "protein kNN", sample_median = "sample median",
             minprob = "MinProb", mindet = "MinDet")
val[, mlab := factor(lab_map[method], levels = lab_map[c("cross_protease","protein_knn","protein_mean",
                                                         "sample_median","minprob","mindet")])]
val[, Condition := factor(Condition, levels = CONDS)]
val[, is_cross := method == "cross_protease"]
pB <- ggplot(val, aes(mlab, median_abs_error, fill = is_cross)) +
  geom_col(width = 0.7, color = "grey25", linewidth = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4, color = "grey40") +
  geom_text(aes(label = sprintf("%.2f", median_abs_error)), vjust = -0.4, size = 2.8, color = "grey20") +
  facet_wrap(~Condition, nrow = 1) +
  scale_fill_manual(values = c(`TRUE` = "#1b7837", `FALSE` = "#B8B8B8"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "MNAR recovery error (median |log2|)",
       title = "B  Validity: masking recovery under realistic (MNAR) missingness",
       subtitle = "Hide the low-abundance cells, refill, compare to truth (Krakatoa recipient). Cross-protease is\nthe most accurate; the dashed line is 2-fold error. The downshifted draws (MinProb/MinDet) fail.") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8.4, color = "grey30"),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 25, hjust = 1),
        strip.background = element_rect(fill = "grey92", color = NA), strip.text = element_text(face = "bold"))

fig <- (pA / pB) +
  plot_annotation(
    title = "Cross-protease imputation: what it buys, and that it is valid",
    caption = paste0("Discovery cohort, threshold 50% completeness. Imputation = per-protein median log2 offset from the other ",
                     "protease's measurement of the\nsame run (>=5 paired runs); measured values unchanged; the two proteases are never averaged. ",
                     "Counts: output/06_supplementary/tables."),
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.caption = element_text(hjust = 0, size = 8, color = "grey35")))

ggsave(file.path(FIG, "imputation_gain_and_validity.png"), fig, width = 9.5, height = 8.6, dpi = 600, bg = "white")
ggsave(file.path(FIG, "imputation_gain_and_validity.pdf"), fig, width = 9.5, height = 8.6, bg = "white")
cat("Saved: figures/imputation_gain_and_validity.png + .pdf\n")
cat("\nGain (tau>=50%, strict):\n"); print(gain[, .(Condition, channel, before, after, delta)])
cat("\nValidity (MNAR, Krakatoa recipient) median |log2|:\n")
print(val[order(Condition, median_abs_error), .(Condition, method, median_abs_error = round(median_abs_error,3))])
