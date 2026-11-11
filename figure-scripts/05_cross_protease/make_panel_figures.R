# =============================================================================
# make_panel_figures.R — Panel A in the manuscript's own layout, with updated
# numbers, plus the before/after-imputation counterpart.
#
# Layout reproduced from the draft figure:
#   [ Krakatoa count strip ]
#   [ Krakatoa violin      ]
#   [ "Combined proteins" title + count strip ]
#   [ Trypsin violin       ]
#   [ Trypsin count strip  ]
#
# Numbers come from cv_common.R, shared with make_manuscript_figures.R, so the
# panel can never disagree with the standalone figures.
#
# Two definitions are written to separate files:
#   published/     CV = sd of 3 day-means; filter = all 3 days detected
#               (reproduces the published Panel A exactly)
#   corrected/  CV = injection noise removed from the day-means;
#               filter = all 24 injections detected  (makes "100% completeness" true)
#
# NOTE ON THE COMBINED ROW.  The draft's combined row averaged the two
# proteases' raw intensities. That mixes two different intensity scales and is
# net-destructive (it loses 66 protein x samples that pass on a single protease
# and rescues only 7). Here "Combined" is the UNION: a protein counts if EITHER
# protease quantifies it well. Nothing is averaged.
#
# Usage: Rscript figure-scripts/05_cross_protease/make_panel_figures.R
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

.sd <- function() { a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f)) normalizePath(dirname(f)) else getwd() }
source(file.path(.sd(), "..", "config.R"))   # DATA_DIR, OUT_ROOT, step_dirs()
source(file.path(.sd(), "cv_common.R"))      # same data + counts as make_manuscript_figures.R

OUT <- file.path(step_dirs("05_cross_protease")$figures, "fig4a_panels")
for (d in c(file.path(OUT, "published"), file.path(OUT, "corrected")))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

TIER_MAIN <- "strict"

# y-axis per definition, shared by both panels of that definition
YMAX <- sapply(c("published","corrected"), function(def) {
  v <- unlist(lapply(c("Krakatoa","Trypsin"), function(pr)
        c(pdata(pr, def, "b")$CV, pdata(pr, def, "a", TIER_MAIN)$CV)))
  ceiling(max(v, na.rm = TRUE) / 10) * 10
})
YBREAKS <- lapply(YMAX, function(m) seq(0, m, by = if (m > 120) 50 else 25))
cat(sprintf("Y-axis: published 0-%.0f, corrected 0-%.0f\n", YMAX[["published"]], YMAX[["corrected"]]))

# ------------------------------------------------------------------ pieces ---
# a horizontal strip of count boxes, one per sample (the draft's "table row")
strip <- function(vals, fill, txt = "white", vals2 = NULL, fill2 = NULL) {
  if (is.null(vals2)) {
    d <- data.table(x = seq_along(samples), y = 1, lab = vals, f = fill)
  } else {
    d <- rbind(data.table(x = seq_along(samples) - 0.21, y = 1, lab = vals,  f = fill2),
               data.table(x = seq_along(samples) + 0.21, y = 1, lab = vals2, f = fill))
  }
  w <- if (is.null(vals2)) 0.94 else 0.40
  ggplot(d, aes(x, y)) +
    geom_tile(aes(fill = f), width = w, height = 0.86) +
    geom_text(aes(label = lab), color = txt, fontface = "bold",
              size = if (is.null(vals2)) 4.6 else 3.6) +
    scale_fill_identity() +
    scale_x_continuous(limits = c(0.5, length(samples) + 0.5), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void() + theme(plot.margin = margin(1, 2, 1, 2))
}

# violin row; `paired` adds the after-imputation distribution alongside
violin_row <- function(pr, def, ylab, paired = FALSE, tier = TIER_MAIN) {
  b <- pdata(pr, def, "b")[, cond := "Before"]
  d <- if (paired) rbind(b, pdata(pr, def, "a", tier)[, cond := "After"]) else b
  d[, `:=`(SampleID = factor(SampleID, levels = samples),
           cond = factor(cond, levels = c("Before", "After")))]
  p <- ggplot(d, aes(SampleID, CV))
  if (paired) {
    p <- p + geom_violin(aes(fill = cond), trim = TRUE, scale = "width", color = "black",
                         linewidth = 0.3, position = position_dodge(0.85)) +
             geom_boxplot(aes(group = interaction(SampleID, cond)), width = 0.10,
                          outlier.shape = NA, fill = "white", color = "black",
                          linewidth = 0.25, position = position_dodge(0.85)) +
             scale_fill_manual(values = c(Before = LIGHT[[pr]], After = COL[[pr]]), guide = "none")
  } else {
    p <- p + geom_violin(fill = COL[[pr]], trim = TRUE, scale = "width", color = "black",
                         linewidth = 0.3) +
             geom_boxplot(width = 0.11, outlier.shape = NA, fill = "white",
                          color = "black", linewidth = 0.25)
  }
  p + geom_hline(yintercept = CV_THRESHOLD, linetype = "dashed", linewidth = 0.5) +
    scale_y_continuous(limits = c(0, YMAX[[def]]), breaks = YBREAKS[[def]], expand = c(0, 0)) +
    labs(x = NULL, y = ylab) +
    theme_classic(base_size = 13) +
    theme(axis.text.y = element_text(color = "black", size = 11),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.title.y = element_text(face = "bold", size = 12),
          plot.margin = margin(2, 4, 2, 4))
}

cnt <- function(pr, def, which = "b", tier = TIER_MAIN) {
  d <- pdata(pr, def, which, tier)
  r <- d[CV <= CV_THRESHOLD, .N, by = SampleID]
  x <- r$N[match(samples, r$SampleID)]; x[is.na(x)] <- 0L; as.integer(x)
}

sample_axis <- function() {
  ggplot(data.table(x = seq_along(samples), lab = paste0("S", samples)), aes(x, 1)) +
    geom_text(aes(label = lab), size = 4) +
    scale_x_continuous(limits = c(0.5, length(samples) + 0.5), expand = c(0, 0)) +
    theme_void() + theme(plot.margin = margin(0, 2, 2, 2))
}
title_row <- function(txt, size = 11.5) {
  ggplot() + annotate("text", 0, 0, label = txt, fontface = "bold", size = size/2.6) +
    theme_void() + theme(plot.margin = margin(4, 2, 1, 2))
}

comb_label <- function(def)
  sprintf("Combined proteins (<%d%% CV, %s), quantified by either protease",
          CV_THRESHOLD, if (def == "corrected") "100% completeness" else "all 3 days")

# ------------------------------------------------------------------ panels ---
build <- function(def, paired, tier = TIER_MAIN, tag = NULL) {
  if (is.null(tag)) tag <- if (paired) sprintf("imputation_%s", tier) else "A"
  kr_b <- cnt("Krakatoa", def, "b"); tr_b <- cnt("Trypsin", def, "b"); cb_b <- cnt("Combined", def, "b")
  if (paired) {
    kr_a <- cnt("Krakatoa", def, "a", tier); tr_a <- cnt("Trypsin", def, "a", tier)
    cb_a <- cnt("Combined", def, "a", tier)
    s_kr <- strip(kr_b, COL[["Krakatoa"]], vals2 = kr_a, fill2 = LIGHT[["Krakatoa"]])
    s_cb <- strip(cb_b, COL[["Combined"]], vals2 = cb_a, fill2 = LIGHT[["Combined"]])
    s_tr <- strip(tr_b, COL[["Trypsin"]],  vals2 = tr_a, fill2 = LIGHT[["Trypsin"]])
  } else {
    s_kr <- strip(kr_b, COL[["Krakatoa"]]); s_cb <- strip(cb_b, COL[["Combined"]])
    s_tr <- strip(tr_b, COL[["Trypsin"]])
  }
  ttl <- if (paired) "Cross-protease imputation" else "Inter-day variability by protease"
  sub <- paste0(
    if (def == "corrected")
      sprintf("Inter-day CV with injection noise removed from the day-means; proteins detected in all 24 injections. Threshold %d%%.", CV_THRESHOLD)
    else
      sprintf("Inter-day CV as sd of the three day-means; proteins detected on all 3 days. Threshold %d%%.", CV_THRESHOLD),
    if (paired) sprintf("\nPale = measured only, solid = after cross-protease imputation (tier '%s': %s). Counts are proteins below the threshold.", tier, TIERS[[tier]])
    else "\nCounts are proteins below the threshold. Combined = union of the two proteases; the proteases are never averaged.")

  p <- title_row(ttl, 15) / s_kr /
       violin_row("Krakatoa", def, "Krakatoa inter-day %CV", paired, tier) /
       title_row(comb_label(def)) / s_cb /
       violin_row("Trypsin", def, "Trypsin inter-day %CV", paired, tier) / s_tr / sample_axis() +
       plot_layout(heights = c(0.30, 0.34, 3.1, 0.26, 0.34, 3.1, 0.34, 0.22)) +
       plot_annotation(caption = sub,
         theme = theme(plot.caption = element_text(hjust = 0, size = 8.6, color = "grey30",
                                                   margin = margin(6, 0, 0, 0))))
  nm <- file.path(OUT, def, sprintf("panel_%s", tag))
  ggsave(paste0(nm, ".png"), p, width = 8.4, height = 9.2, dpi = 600, bg = "white")
  ggsave(paste0(nm, ".pdf"), p, width = 8.4, height = 9.2, bg = "white")
  cat(sprintf("  %s/panel_%s.png + .pdf\n", def, tag))
  data.table(Definition = def, Figure = tag, Tier = if (paired) tier else NA_character_,
             SampleID = samples,
             Krakatoa_before = kr_b, Krakatoa_after = if (paired) kr_a else NA_integer_,
             Combined_before = cb_b, Combined_after = if (paired) cb_a else NA_integer_,
             Trypsin_before  = tr_b, Trypsin_after  = if (paired) tr_a else NA_integer_)
}

cat("\n--- Panel figures ---\n")
res <- rbindlist(lapply(c("published","corrected"), function(d) {
  out <- rbindlist(c(list(build(d, FALSE)),
                     lapply(names(TIERS), function(tr) build(d, TRUE, tr))))
  # also write the original filename, pointing at the strict tier; not tabulated
  invisible(build(d, TRUE, "strict", tag = "imputation"))
  out
}))
fwrite(res, file.path(TBL, "panel_figure_numbers.csv"))
cat("\nNumbers written to panel_figure_numbers.csv\n")
cat("\n=== totals by tier ===\n")
print(dcast(res[!is.na(Tier), .(Krakatoa = sum(Krakatoa_after), Combined = sum(Combined_after),
                                Trypsin = sum(Trypsin_after)), by = .(Definition, Tier)],
            Definition + Tier ~ ., value.var = c("Krakatoa","Combined","Trypsin")))
cat("\n=== ungated, per sample ===\n")
print(res[Tier == "ungated", .(Definition, SampleID, Krakatoa_before, Krakatoa_after,
                               Combined_before, Combined_after, Trypsin_before, Trypsin_after)])
