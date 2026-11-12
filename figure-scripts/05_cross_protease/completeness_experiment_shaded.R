# =============================================================================
# completeness_experiment_shaded.R — Fig 4C (Oct 2026): the
# completeness_experiment figures with the "with imputation" bars HATCHED so they
# read as visually distinct from the measured bars. Same colours, bins and counts.
#
# Reads the counts written by completeness_experiment.R (tables/completeness_
# experiment_{mix,repro}.csv) — no recomputation; guards that the known totals
# (100% bin: Mix30 424->573, Mix20 344->471, Repro 40->73) are what is plotted.
# Writes NEW files (*_shaded.*); the originals are left untouched.
#
# Usage: Rscript figure-scripts/05_cross_protease/completeness_experiment_shaded.R
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(grid) })

.sd <- function() { a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f)) normalizePath(dirname(f)) else getwd() }
source(file.path(.sd(), "..", "config.R"))   # DATA_DIR, OUT_ROOT, step_dirs()
D <- step_dirs("05_cross_protease"); FIG <- D$figures; TBL <- D$tables
COMP_LEVELS <- c("0-20%","20-50%","50-80%","80-95%","95-100%","100%")
COMP_COLS <- setNames(colorRampPalette(c("red2","gold","green3"))(6), COMP_LEVELS)

mixb <- fread(file.path(TBL, "completeness_experiment_mix.csv"))
repb <- fread(file.path(TBL, "completeness_experiment_repro.csv"))
chk <- rbind(mixb, repb)[Bin == "100%"][order(Group, phase)]
stopifnot(identical(chk[, paste(Group, phase, N)],
  c("Mix20 Imputed 471", "Mix20 Measured 344", "Mix30 Imputed 573", "Mix30 Measured 424",
    "Repro Imputed 73", "Repro Measured 40")))
cat("Count guard: PASSED\n")

# Diagonal hatch tile (white stripes, semi-transparent) repeated across the bar
hatch <- function(spacing_mm = 2.2, lwd = 1.1, col = "white", alpha = 0.75, bg = NA) {
  # tile drawn in absolute mm (not npc) so stripes stay 45 degrees whatever the bar's aspect ratio
  d <- spacing_mm
  stripes <- segmentsGrob(x0 = unit(c(-d, 0, d), "mm"), y0 = unit(0, "mm"),
                          x1 = unit(c(0, d, 2 * d), "mm"), y1 = unit(d, "mm"),
                          gp = gpar(col = col, alpha = alpha, lwd = lwd, lineend = "square"))
  g <- if (is.na(bg)) stripes else
    grobTree(rectGrob(x = 0, y = 0, width = unit(d, "mm"), height = unit(d, "mm"),
                      just = c("left", "bottom"), gp = gpar(fill = bg, col = NA)), stripes)
  pattern(g, x = 0, y = 0, width = unit(d, "mm"), height = unit(d, "mm"),
          just = c("left", "bottom"), extend = "repeat")
}

plot_shaded <- function(bins, group_levels, xlabels, ylab, title, label_size = 4.2) {
  bins <- copy(bins)
  bins[, Bin := factor(Bin, levels = COMP_LEVELS)]
  bins[, xcol := factor(paste(Group, phase),
                        levels = as.vector(t(outer(group_levels, c("Measured","Imputed"), paste))))]
  bins <- bins[order(xcol, Bin)]
  imp <- bins[phase == "Imputed"]
  # label midpoints, stacked bottom-up from 100% (matches position_stack order below)
  lab <- bins[order(xcol, -as.integer(Bin))][, ymid := cumsum(N) - N / 2, by = xcol][N >= 20]
  ggplot(bins, aes(xcol, N)) +
    geom_col(aes(fill = Bin), width = 0.72, colour = "grey15", linewidth = 0.15,
             position = position_stack(reverse = FALSE)) +                      # green (100%) at bottom
    # hatch overlay on the imputed bars only (one block per bar, stacked total)
    geom_col(data = imp[, .(N = sum(N)), by = xcol], aes(colour = "With cross-protease\nimputation"),
             fill = hatch(), width = 0.72, linewidth = 0.4) +
    # counts on a soft white chip so the hatch never runs through the digits
    geom_label(data = lab, aes(xcol, ymid, label = N), size = label_size, fontface = "bold",
               colour = "black",
               fill = scales::alpha("white", 0.8), linewidth = 0, label.r = unit(0.8, "mm"),
               label.padding = unit(0.5, "mm")) +
    scale_fill_manual(values = COMP_COLS, breaks = rev(COMP_LEVELS), name = "Data\ncompleteness") +
    scale_colour_manual(values = "grey15", name = NULL,
                        guide = guide_legend(override.aes = list(fill = list(hatch(bg = "grey55", alpha = 0.9)), colour = "grey15"),
                                             order = 2)) +
    scale_x_discrete(labels = xlabels) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
    labs(x = NULL, y = ylab, title = title) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"), axis.text.x = element_text(size = 10),
          axis.title.y = element_text(face = "bold"), legend.position = "right")
}

save_both <- function(p, stem, w, h) {
  ggsave(file.path(FIG, paste0(stem, ".png")), p, width = w, height = h, dpi = 600, bg = "white",
         device = ragg::agg_png)
  ggsave(file.path(FIG, paste0(stem, ".pdf")), p, width = w, height = h, bg = "white",
         device = grDevices::pdf, useDingbats = FALSE)
  cat("Saved:", stem, "(.png, .pdf)\n")
}

pm <- plot_shaded(mixb, c("Mix30","Mix20"),
                  c("Mix30 Measured" = "Mix 30\nmeasured", "Mix30 Imputed" = "Mix 30\nwith imputation",
                    "Mix20 Measured" = "Mix 20\nmeasured", "Mix20 Imputed" = "Mix 20\nwith imputation"),
                  "Protein groups (all in experiment)", "Discovery (Mix 30 / Mix 20)")
save_both(pm, "fig_mix_completeness_experiment_shaded", 8, 6)

pr <- plot_shaded(repb, "Repro",
                  c("Repro Measured" = "Selected", "Repro Imputed" = "With imputation"),
                  "Protein groups (all in experiment)", "Reproducibility", label_size = 5.2)
save_both(pr, "fig_repro_completeness_experiment_shaded", 5.5, 6)
cat("Done.\n")
