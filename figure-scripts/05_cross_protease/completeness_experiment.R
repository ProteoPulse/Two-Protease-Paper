# =============================================================================
# completeness_experiment.R — completeness of ALL proteins in the experiment,
# before / after cross-protease imputation, as a single combined distribution.
#
# Combined experiment matrix, one value per protein per run:
#   * trypsin-only proteins  -> trypsin value
#   * Krakatoa-only proteins -> Krakatoa value
#   * shared proteins        -> AVERAGE of the two channels. An averaged value
#     needs both proteases, so before imputation a shared protein is quantified in
#     a run only where BOTH searches saw it (intersection); strict cross-protease
#     imputation supplies the missing channel, so after imputation the average can
#     be formed wherever EITHER saw it (union) -> completeness rises. Trypsin-only
#     and Krakatoa-only proteins have no donor and do not change.
#
# Separate from the existing figure set. Panel A asserted first.
# Usage: Rscript figure-scripts/05_cross_protease/completeness_experiment.R
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

.sd <- function() { a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f)) normalizePath(dirname(f)) else getwd() }
source(file.path(.sd(), "..", "config.R"))   # DATA_DIR, OUT_ROOT, step_dirs()
D <- step_dirs("05_cross_protease"); FIG <- D$figures; TBL <- D$tables
PG_TRY <- file.path(SEARCH, "Trypsin_Library_Search",  "trypsin_library_search_report.pg_matrix.tsv")
PG_KRA <- file.path(SEARCH, "Krakatoa_Library_Search", "report.pg_matrix.tsv")
REPRO  <- file.path(SEARCH, "Reproducibility_Search")
RP_TRY <- file.path(REPRO, "Trypsin Library",  "report.pg_matrix.tsv")
RP_HTA <- file.path(REPRO, "Krakatoa Library", "report.pg_matrix.tsv")
ANN <- c("Genes","Protein.Group","Protein.Ids","Protein.Names",
         "First.Protein.Description","N.Sequences","N.Proteotypic.Sequences")
MIN_PAIRED <- 5L
COMP_LEVELS <- c("0-20%","20-50%","50-80%","80-95%","95-100%","100%")
COMP_COLS <- setNames(colorRampPalette(c("red2","gold","green3"))(6), COMP_LEVELS)
bin_c <- function(pct) { if (!is.finite(pct)) return(NA_character_)
  if (pct>=100) "100%" else if (pct>=95) "95-100%" else if (pct>=80) "80-95%" else
  if (pct>=50) "50-80%" else if (pct>=20) "20-50%" else if (pct>0) "0-20%" else NA_character_ }

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
load_disc <- function(path) {
  L <- melt_pg(fread(path, sep = "\t", header = TRUE, showProgress = FALSE))
  b <- sub("\\.raw$", "", basename(gsub("\\\\", "/", L$Sample)), ignore.case = TRUE)
  is20 <- grepl("20min", b)
  L[, Group := fifelse(grepl("Plasma_Mix", b) & is20, "Mix20",
               fifelse(grepl("Plasma_Mix", b) & !is20, "Mix30", "other"))][]
}
load_repro <- function(path) {
  L <- melt_pg(fread(path, sep = "\t", header = TRUE, showProgress = FALSE))
  b <- sub("\\.raw$", "", basename(gsub("\\\\", "/", L$Sample)), ignore.case = TRUE)
  p <- tstrsplit(b, "_", fixed = TRUE); n <- length(p)
  L[, Day := as.integer(p[[n-2]])][, SampleID := p[[n-1]]][]
}
assert_panel_a <- function() {
  cnt <- function(L) {
    dm <- L[is.finite(Intensity), .(DayMean = mean(Intensity), .N), by = .(FeatureID, SampleID, Day)]
    ed <- L[, .(ED = uniqueN(Day)), by = SampleID]
    cv <- dm[, .(CV = if (.N < 2) NA_real_ else 100*sd(DayMean)/mean(DayMean), MD = .N), by = .(FeatureID, SampleID)]
    cv <- merge(cv, ed, by = "SampleID"); cv[, comp := MD == ED]
    so <- as.character(sort(as.integer(unique(L$SampleID))))
    w <- cv[comp & MD >= 2 & is.finite(CV) & CV <= 20, .N, by = SampleID]
    v <- w$N[match(so, w$SampleID)]; v[is.na(v)] <- 0L; as.integer(v)
  }
  ok <- identical(cnt(load_repro(RP_TRY)), c(120L,107L,127L,123L,123L,123L,146L,131L)) &&
        identical(cnt(load_repro(RP_HTA)), c(92L,80L,88L,89L,108L,101L,91L,82L))
  cat("Panel A guardrail:", if (ok) "PASSED" else "FAILED", "\n"); if (!ok) stop("Panel A failed")
}

# offset-estimable proteins (>=MIN_PAIRED co-observed) — needed to know if a shared
# protein's missing channel can be imputed
offset_ok <- function(tr, kr) {
  m <- merge(tr[is.finite(Intensity), .(FeatureID, Sample)], kr[is.finite(Intensity), .(FeatureID, Sample)],
             by = c("FeatureID","Sample"))
  m[, .N, by = FeatureID][N >= MIN_PAIRED, FeatureID]
}

# combined experiment completeness (before/after) for one run set
experiment_bins <- function(tr, kr, label) {
  nr <- uniqueN(rbind(tr[, .(Sample)], kr[, .(Sample)]))
  ok <- offset_ok(tr, kr)
  t_meas <- tr[is.finite(Intensity), .(FeatureID, Sample)][, tm := TRUE]
  k_meas <- kr[is.finite(Intensity), .(FeatureID, Sample)][, km := TRUE]
  grid <- merge(t_meas, k_meas, by = c("FeatureID","Sample"), all = TRUE)
  grid[is.na(tm), tm := FALSE][is.na(km), km := FALSE]
  seen <- grid[, .(t_seen = any(tm), k_seen = any(km)), by = FeatureID]
  seen[, cat := fifelse(t_seen & !k_seen, "trypsin_only",
              fifelse(!t_seen & k_seen, "krakatoa_only", "shared"))]
  grid <- merge(grid, seen[, .(FeatureID, cat)], by = "FeatureID")
  grid[, has_off := FeatureID %in% ok]
  # present before / after per (protein, run)
  grid[, before := fifelse(cat == "trypsin_only", tm,
                   fifelse(cat == "krakatoa_only", km, tm & km))]
  # after: shared needs both channels; imputation supplies the missing one if offset exists
  grid[, t_after := tm | (has_off & km)]
  grid[, k_after := km | (has_off & tm)]
  grid[, after := fifelse(cat == "trypsin_only", tm,
                  fifelse(cat == "krakatoa_only", km, t_after & k_after))]
  comp <- grid[, .(before = sum(before)/nr, after = sum(after)/nr), by = FeatureID]
  binc <- function(v, phase) data.table(FeatureID = comp$FeatureID,
    Bin = vapply(100*v, bin_c, character(1)), phase = phase)
  b <- rbind(binc(comp$before, "Measured"), binc(comp$after, "Imputed"))
  b <- b[!is.na(Bin), .N, by = .(Bin, phase)][, Group := label]
  list(bins = b, ncat = seen[, .N, by = cat])
}

plot_exp <- function(bins, group_levels, title) {
  bins[, Bin := factor(Bin, levels = COMP_LEVELS)]
  bins[, phase := factor(phase, levels = c("Measured","Imputed"), labels = c("Measured","Imputed"))]
  bins[, xcol := factor(paste(Group, phase), levels = as.vector(t(outer(group_levels, c("Measured","Imputed"), paste))))]
  ggplot(bins, aes(xcol, N, fill = Bin)) +
    geom_col(width = 0.72, colour = "grey15", linewidth = 0.15, position = position_stack(reverse = FALSE)) +   # green (100%) at bottom
    geom_text(aes(label = fifelse(N >= 20, as.character(N), "")),
              position = position_stack(vjust = 0.5, reverse = FALSE), size = 2.6, colour = "grey10") +
    scale_fill_manual(values = COMP_COLS, breaks = rev(COMP_LEVELS), name = "Data\ncompleteness") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
    labs(x = NULL, y = "Protein groups (all in experiment)", title = title,
         subtitle = "Combined = trypsin-only + Krakatoa-only + shared (averaged; needs both proteases). Cross-protease imputation\nsupplies the missing channel of shared proteins, so their averaged value can be formed in more runs. Counts >=20 shown.") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 8.3, colour = "grey30"),
          axis.text = element_text(color = "black"), axis.text.x = element_text(size = 9),
          axis.title.y = element_text(face = "bold"), legend.position = "right")
}

cat("Loading + Panel A...\n"); assert_panel_a()

# --- Mix (discovery): per condition ---
TRd <- load_disc(PG_TRY); KRd <- load_disc(PG_KRA)
mixres <- lapply(c("Mix30","Mix20"), function(cd) experiment_bins(TRd[Group==cd], KRd[Group==cd], cd))
mixb <- rbindlist(lapply(mixres, `[[`, "bins"))
fwrite(mixb, file.path(TBL, "completeness_experiment_mix.csv"))
cat("Mix protein categories:\n"); for (i in 1:2) { cat(" ", c("Mix30","Mix20")[i], ": "); print(mixres[[i]]$ncat) }
pm <- plot_exp(mixb, c("Mix30","Mix20"), "Mix (discovery): all proteins, completeness before / after imputation")
ggsave(file.path(FIG, "fig_mix_completeness_experiment.png"), pm, width = 8, height = 6, dpi = 600, bg = "white")
ggsave(file.path(FIG, "fig_mix_completeness_experiment.pdf"), pm, width = 8, height = 6, bg = "white")
cat("Saved: fig_mix_completeness_experiment.png\n")

# --- Reproducibility: whole experiment (all 192 runs) ---
TRr <- load_repro(RP_TRY); KRr <- load_repro(RP_HTA)
rres <- experiment_bins(TRr, KRr, "Repro")
fwrite(rres$bins, file.path(TBL, "completeness_experiment_repro.csv"))
cat("Repro protein categories:\n"); print(rres$ncat)
pr <- plot_exp(rres$bins, "Repro", "Reproducibility: all proteins, completeness before / after imputation")
ggsave(file.path(FIG, "fig_repro_completeness_experiment.png"), pr, width = 7, height = 6, dpi = 600, bg = "white")
ggsave(file.path(FIG, "fig_repro_completeness_experiment.pdf"), pr, width = 7, height = 6, bg = "white")
cat("Saved: fig_repro_completeness_experiment.png\n\nDone.\n")
