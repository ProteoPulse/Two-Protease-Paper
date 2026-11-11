# =============================================================================
# make_manuscript_figures.R — standalone, publication-ready CV figures.
#
# Produces, for BOTH definitions (choose per figure when assembling):
#
#   published/     CV = 100*sd(3 day-means)/mean      filter = all 3 days detected
#               (reproduces Panel A exactly)
#   corrected/  CV = injection noise removed:      filter = all 24 injections
#               sigma2_day = var(day means) - MS_within/8, floored at 0
#
#   1. cv_<protease>_baseline.(png|pdf)        measured data only
#   2. cv_<protease>_<tier>.(png|pdf)          before vs after imputation, x4 tiers
#   3. stepwise_strictness.(png|pdf)           counts across the tier ladder
#   4. cv_both_proteases_<tier>.(png|pdf)      stacked two-protease version
#
# Imputation: cross-protease gap fill into each protease's OWN matrix (never
# averaged); donor rescaled by the protein's median log2 offset; measured values
# never altered. Tiers control how many of the 3 day-means may be built entirely
# from filled values:  strict 0 | onefill <=1 | partial <=2 | ungated any.
#
# Usage: Rscript figure-scripts/05_cross_protease/make_manuscript_figures.R
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

.sd <- function() { a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f)) normalizePath(dirname(f)) else getwd() }
source(file.path(.sd(), "..", "config.R"))   # DATA_DIR, OUT_ROOT, step_dirs()
source(file.path(.sd(), "cv_common.R"))      # data, CV metrics, Panel A gate, counts table

# ------------------------------------------------------------- plotting -----
# One y-axis per definition, shared by every figure of that definition so the
# baseline, the four tiers and both proteases can be compared side by side.
YMAX <- sapply(c("published","corrected"), function(def) {
  cc <- cvcol(def)
  v <- unlist(lapply(c("Krakatoa","Trypsin"), function(pr)
        c(elig(M[[pr]]$b, def)[[cc]], elig(M[[pr]]$a, def)[[cc]])))
  ceiling(max(v, na.rm = TRUE) / 10) * 10
})
YLAB_AT <- YMAX * 1.045
YLIM    <- YMAX * 1.10
YBREAKS <- lapply(YMAX, function(m) seq(0, m, by = if (m > 120) 50 else 25))
cat(sprintf("\nY-axis: published 0-%.0f, corrected 0-%.0f\n", YMAX[["published"]], YMAX[["corrected"]]))

base_theme <- theme_classic(base_size = 12) +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold", size = 11),
        legend.position = "top", legend.text = element_text(size = 9),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8.5, color = "grey30"),
        plot.margin = margin(6, 8, 6, 8))
ylab_for <- function(pr, def) {
  if (pr == "Combined") "Combined inter-day %CV  (better protease)"
  else sprintf("%s inter-day %%CV", pr)
}

# --- 1. standalone baseline violin (one protease, measured data only) --------
fig_baseline <- function(pr, def) {
  e <- pdata(pr, def, "b")
  d <- e[, .(SampleID = factor(SampleID, levels = samples), CV)]
  lab <- data.table(SampleID = factor(samples, levels = samples),
                    n = vv(e[CV <= CV_THRESHOLD]))
  ggplot(d, aes(SampleID, CV)) +
    geom_violin(fill = COL[[pr]], trim = TRUE, scale = "width", color = "black",
                linewidth = 0.3, alpha = 0.9) +
    geom_boxplot(width = 0.11, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.25) +
    geom_hline(yintercept = CV_THRESHOLD, linetype = "dashed", linewidth = 0.5) +
    geom_label(data = lab, aes(SampleID, YLAB_AT[[def]], label = n), inherit.aes = FALSE,
               fill = COL[[pr]], color = "white", fontface = "bold", size = 3.4,
               label.size = 0, label.r = unit(2,"pt"), label.padding = unit(2.4,"pt")) +
    scale_y_continuous(limits = c(0, YLIM[[def]]), breaks = YBREAKS[[def]], expand = c(0,0)) +
    labs(x = "Sample", y = ylab_for(pr, def),
         subtitle = if (def == "corrected")
           sprintf("Proteins with CV <%d%% and all 24 injections detected (injection noise removed from the day-means)", CV_THRESHOLD)
         else sprintf("Proteins with CV <%d%% detected on all 3 days", CV_THRESHOLD)) +
    base_theme
}

# --- 2. standalone before/after imputation (one protease, one tier) ---------
fig_tier <- function(pr, def, tier) {
  eb <- pdata(pr, def, "b"); ea <- pdata(pr, def, "a", tier)
  d <- rbind(eb[, .(SampleID, CV, cond = "Before")],
             ea[, .(SampleID, CV, cond = "After")])
  d[, `:=`(SampleID = factor(SampleID, levels = samples),
           cond = factor(cond, levels = c("Before","After")))]
  lab <- rbind(data.table(SampleID = samples, cond = "Before", n = vv(eb[CV <= CV_THRESHOLD])),
               data.table(SampleID = samples, cond = "After",  n = vv(ea[CV <= CV_THRESHOLD])))
  lab[, `:=`(SampleID = factor(SampleID, levels = samples),
             cond = factor(cond, levels = c("Before","After")))]
  gain <- mean(lab[cond=="After", n] - lab[cond=="Before", n])
  ggplot(d, aes(SampleID, CV, fill = cond)) +
    geom_violin(trim = TRUE, scale = "width", color = "black", linewidth = 0.3,
                alpha = 0.92, position = position_dodge(0.85)) +
    geom_boxplot(width = 0.11, outlier.shape = NA, fill = "white", color = "black",
                 linewidth = 0.25, position = position_dodge(0.85)) +
    geom_hline(yintercept = CV_THRESHOLD, linetype = "dashed", linewidth = 0.5) +
    geom_label(data = lab, aes(SampleID, YLAB_AT[[def]], label = n, fill = cond), inherit.aes = FALSE,
               position = position_dodge(0.85), color = "white", fontface = "bold", size = 3.1,
               label.size = 0, label.r = unit(2,"pt"), label.padding = unit(2.2,"pt"),
               show.legend = FALSE) +
    scale_fill_manual(values = c(Before = LIGHT[[pr]], After = COL[[pr]]),
                      labels = c("Measured only", "After cross-protease imputation")) +
    scale_y_continuous(limits = c(0, YLIM[[def]]), breaks = YBREAKS[[def]], expand = c(0,0)) +
    labs(x = "Sample", y = ylab_for(pr, def), fill = NULL,
         subtitle = sprintf("%s  |  tier '%s': %s  |  mean gain %+.1f proteins/sample",
                            if (def == "corrected") sprintf("CV <%d%%, all 24 injections", CV_THRESHOLD)
                            else sprintf("CV <%d%%, all 3 days", CV_THRESHOLD),
                            tier, TIERS[[tier]], gain)) +
    base_theme
}

sv <- function(p, dir, name, w = 8.6, h = 5.2) {
  ggsave(file.path(OUT, "figures", dir, paste0(name, ".png")), p, width = w, height = h,
         dpi = PUB_DPI, bg = "white")
  ggsave(file.path(OUT, "figures", dir, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
}

cat("\n--- Writing figures ---\n")
for (def in c("published","corrected")) {
  for (pr in PROTEASES) {
    sv(fig_baseline(pr, def), def, sprintf("cv_%s_baseline", tolower(pr)))
    for (tr in names(TIERS)) sv(fig_tier(pr, def, tr), def, sprintf("cv_%s_%s", tolower(pr), tr))
  }
  # stacked two-protease version per tier
  for (tr in names(TIERS)) {
    p <- (fig_tier("Krakatoa", def, tr) + labs(x = NULL)) / fig_tier("Trypsin", def, tr)
    sv(p, def, sprintf("cv_both_proteases_%s", tr), w = 9.2, h = 9.4)
    p3 <- (fig_baseline("Krakatoa", def) + labs(x = NULL)) /
          (fig_tier("Combined", def, tr) + labs(x = NULL)) /
          (fig_baseline("Trypsin", def))
    sv(p3, def, sprintf("cv_panelA_layout_%s", tr), w = 9.2, h = 13.2)
    pAll <- (fig_tier("Krakatoa", def, tr) + labs(x = NULL)) /
            (fig_tier("Combined", def, tr) + labs(x = NULL)) /
            (fig_tier("Trypsin", def, tr))
    sv(pAll, def, sprintf("cv_all_three_%s", tr), w = 9.2, h = 13.2)
  }
  # --- 3. stepwise strictness ladder ---
  s <- allc[Definition == def, .(total = sum(n)), by = .(Protease, Tier)]
  s[, Tier := factor(Tier, levels = c("Measured only", names(TIERS)))]
  s[, Protease := factor(Protease, levels = PROTEASES)]
  base <- s[Tier == "Measured only", .(Protease, b = total)]
  s <- merge(s, base, by = "Protease")[, gain := total - b]
  p3 <- ggplot(s, aes(Tier, total, fill = Protease)) +
    geom_col(width = 0.62, position = position_dodge(0.7)) +
    geom_text(aes(label = ifelse(Tier == "Measured only", as.character(total),
                                 sprintf("%d\n(%+d)", total, gain))),
              position = position_dodge(0.7), vjust = -0.25, size = 3.1, lineheight = 0.9) +
    scale_fill_manual(values = COL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = NULL, y = "Well-quantified protein x samples (8 samples)",
         title = "Effect of imputation strictness",
         subtitle = paste0(if (def == "corrected")
                             sprintf("CV <%d%% (corrected) and all 24 injections detected. ", CV_THRESHOLD)
                           else sprintf("CV <%d%% (published definition) and all 3 days detected. ", CV_THRESHOLD),
                           "Tiers cap how many day-means may be built entirely from imputed values.")) +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 12, hjust = 1, color = "black"),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 8.5, color = "grey30"))
  sv(p3, def, "stepwise_strictness", w = 9.4, h = 6)
  cat(sprintf("  %s/ : %d files\n", def, length(list.files(file.path(OUT,"figures",def)))))
}
cat("\nDONE. Figures: ", file.path(OUT, "figures"), "\n", sep = "")
