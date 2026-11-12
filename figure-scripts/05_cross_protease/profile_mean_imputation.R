#!/usr/bin/env Rscript
# profile_mean_imputation.R — Figure 2 (paper): per-sample abundance profile for a
# representative protein, showing the cross-protease fill against the alternatives.
#
# This is the Aug-27 "simplified advisor" profile figure (fig_simplified_advisor.R),
# reproduced for AGRN and CEP95 with ONE change requested for the paper: the flat
# "Median imputed" reference (red squares) is replaced by "Mean imputed" — i.e. the
# recipient protein's MEAN across observed samples, which is the protein_mean baseline
# benchmarked in the masking validation. Everything else (channel selection by
# precursor completeness, per-protein median-offset cross-protease fill, MinProb) is
# unchanged from the original figure.
#
# Channel-selection verification against the FULL Protease Quant Advisor is printed.
#
# Usage: Rscript figure-scripts/05_cross_protease/profile_mean_imputation.R
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
set.seed(20260827)

.sd <- function() { a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f)) normalizePath(dirname(f)) else getwd() }
source(file.path(.sd(), "..", "config.R"))   # DATA_DIR, OUT_ROOT, step_dirs()
OUT  <- step_dirs("05_cross_protease")$figures
DATA <- DATA_DIR

PG_HTA <- file.path(DATA, "Krakatoa_Library_Search", "report.pg_matrix.tsv")
PG_TRY <- file.path(DATA, "Trypsin_Library_Search", "trypsin_library_search_report.pg_matrix.tsv")
PR_HTA <- file.path(DATA, "Krakatoa_Library_Search", "report.pr_matrix.tsv")
PR_TRY <- file.path(DATA, "Trypsin_Library_Search", "report.pr_matrix.tsv")

# palette + theme (verbatim from fig_simplified_advisor.R)
TRY_COL <- "#2C6FB3"; HTA_COL <- "#E08214"
CONC_COL <- "#1B7837"; MP_COL <- "#7B3294"; MEAN_COL <- "#E31A1C"   # red square, now "Mean imputed"
TAU_THRESHOLD <- 0.284
FOCUS <- c(AGRN = "O00468", CEP95 = "Q96GE4")

theme_pub <- function(base = 13)
  theme_classic(base_size = base) +
  theme(plot.title = element_text(face = "bold", size = base + 1),
        plot.subtitle = element_text(colour = "grey30", size = base - 2),
        legend.position = "right", legend.title = element_blank(),
        legend.text = element_text(size = 10), axis.title = element_text(size = base - 1))

canonicalize_pg <- function(pg) vapply(pg, function(x) {
  if (is.na(x) || !nzchar(x)) return(NA_character_)
  parts <- trimws(unlist(strsplit(x, "[;,]"))); paste(sort(unique(parts[nzchar(parts)])), collapse = ";")
}, character(1), USE.NAMES = FALSE)
normalize_run <- function(x) sub("\\.(raw|d|mzML|mzml|wiff|dia)$", "", sub("^.*[\\\\/]", "", as.character(x)), ignore.case = TRUE)
PG_ANNOT <- c("Protein.Group","Protein.Ids","Protein.Names","Genes","First.Protein.Description",
              "N.Sequences","N.Proteotypic.Sequences","Molecular.Weight")
PR_ANNOT <- c("Protein.Group","Protein.Ids","Protein.Names","Genes","First.Protein.Description",
              "Proteotypic","Stripped.Sequence","Modified.Sequence","Precursor.Charge","Precursor.Id",
              "Precursor.Lib.Index","Lib.Q.Value","Lib.PG.Q.Value")
select_mix30_cols <- function(hdr, annot) {
  rc <- setdiff(hdr, annot); nm <- normalize_run(rc)
  keep <- grepl("Plasma_Mix", nm) & !grepl("20min", nm)
  data.table(col = rc[keep], sample_id = nm[keep])
}
load_pg <- function(path, channel) {
  hdr <- names(fread(path, nrows = 0L, sep = "\t", check.names = FALSE)); sel <- select_mix30_cols(hdr, PG_ANNOT)
  d <- fread(path, sep = "\t", select = c(intersect(PG_ANNOT, hdr), sel$col), check.names = FALSE,
             na.strings = c("","NA","NaN","0"))
  L <- melt(d, id.vars = intersect(c("Protein.Group","Genes"), names(d)), measure.vars = sel$col,
            variable.name = "run_col", value.name = "intensity", variable.factor = FALSE)
  L[sel, sample_id := i.sample_id, on = .(run_col = col)]
  L[, pg_canonical := canonicalize_pg(Protein.Group)][, gene := sub(";.*", "", Genes)]
  L[, log2q := log2(as.numeric(intensity))][, channel := channel]
  L[is.finite(log2q), .(pg_canonical, gene, sample_id, log2q, channel)]
}
load_pr <- function(path, channel) {
  hdr <- names(fread(path, nrows = 0L, sep = "\t", check.names = FALSE)); sel <- select_mix30_cols(hdr, PR_ANNOT)
  d <- fread(path, sep = "\t", select = c(intersect(PR_ANNOT, hdr), sel$col), check.names = FALSE,
             na.strings = c("","NA","NaN","0"))
  L <- melt(d, id.vars = intersect(c("Protein.Group","Precursor.Id"), names(d)), measure.vars = sel$col,
            variable.name = "run_col", value.name = "intensity", variable.factor = FALSE)
  L[sel, sample_id := i.sample_id, on = .(run_col = col)]
  L[, pg_canonical := canonicalize_pg(Protein.Group)][, precursor_id := as.character(Precursor.Id)]
  L[, detected := is.finite(as.numeric(intensity)) & as.numeric(intensity) > 0][, channel := channel]
  L[, .(pg_canonical, precursor_id, sample_id, detected, channel)]
}

cat("Loading Mix30 matrices...\n")
pg <- rbind(load_pg(PG_TRY, "trypsin"), load_pg(PG_HTA, "hta"))
samples <- sort(unique(pg$sample_id)); n_samples <- length(samples)
pr <- rbind(load_pr(PR_TRY, "trypsin"), load_pr(PR_HTA, "hta"))

# channel selection by precursor completeness (the simplified-advisor rule)
prec_det  <- pr[, .(detection_rate = sum(detected)/n_samples), by = .(pg_canonical, precursor_id, channel)]
prec_comp <- prec_det[, .(mean_prec_detection = mean(detection_rate)), by = .(pg_canonical, channel)]
comp_wide <- dcast(prec_comp, pg_canonical ~ channel, value.var = "mean_prec_detection", fill = 0)
shared <- comp_wide[trypsin > 0 & hta > 0][, winner := fifelse(trypsin >= hta, "trypsin", "hta")]
gene_lu <- unique(pg[, .(pg_canonical, gene)])[!is.na(gene) & gene != ""]
shared <- merge(shared, gene_lu, by = "pg_canonical", all.x = TRUE)

wide <- dcast(pg[pg_canonical %in% shared$pg_canonical], pg_canonical + sample_id ~ channel, value.var = "log2q")
tau_dt <- wide[, { ok <- is.finite(trypsin) & is.finite(hta)
  .(kendall_tau = if (sum(ok) >= 5L) suppressWarnings(cor(trypsin[ok], hta[ok], method = "kendall")) else NA_real_,
    n_paired = sum(ok)) }, by = pg_canonical]
shared <- merge(shared, tau_dt, by = "pg_canonical", all.x = TRUE)

full_grid <- CJ(pg_canonical = shared$pg_canonical, sample_id = samples, unique = TRUE)
full_wide <- merge(full_grid, wide, by = c("pg_canonical","sample_id"), all.x = TRUE)
full_wide <- merge(full_wide, shared[, .(pg_canonical, winner)], by = "pg_canonical")

# ---- imputation values: cross-protease (median offset), MEAN, MinProb --------
impute_one <- function(prot) {
  info <- shared[pg_canonical == prot]; win <- info$winner
  w <- copy(full_wide[pg_canonical == prot])
  if (win == "trypsin") { w[, rec_val := trypsin]; w[, don_val := hta] }
  else                  { w[, rec_val := hta];     w[, don_val := trypsin] }
  paired <- w[is.finite(rec_val) & is.finite(don_val)]
  if (nrow(paired) < 3L) return(NULL)
  median_offset <- median(paired$rec_val - paired$don_val, na.rm = TRUE)
  measured <- w[is.finite(rec_val), rec_val]
  mu_obs <- mean(measured); sigma_obs <- sd(measured)
  fill <- w[!is.finite(rec_val) & is.finite(don_val)]
  if (nrow(fill) == 0L) return(NULL)
  fill[, cross_protease := don_val + median_offset]
  fill[, minprob := rnorm(.N, mu_obs - 1.8 * sigma_obs, 0.3 * sigma_obs)]
  fill[, mean_imp := mean(measured)]      # <-- MEAN instead of median(measured)
  fill[, .(pg_canonical = prot, sample_id, don_val, cross_protease, minprob, mean_imp)]
}

# ---- profile plot (mean square blocks replace median) ------------------------
make_profile <- function(prot_id, gene_name) {
  info <- shared[pg_canonical == prot_id]; win <- info$winner; tau <- info$kendall_tau
  w <- copy(full_wide[pg_canonical == prot_id])
  if (win == "trypsin") { w[, rec_val := trypsin]; w[, don_val := hta]
    rec_label <- "Trypsin (chosen)"; don_label <- "HTA"; rec_col <- TRY_COL; don_col <- HTA_COL
  } else { w[, rec_val := hta]; w[, don_val := trypsin]
    rec_label <- "HTA (chosen)"; don_label <- "Trypsin"; rec_col <- HTA_COL; don_col <- TRY_COL }
  imp <- impute_one(prot_id)
  w <- merge(w, imp[, .(sample_id, cross_protease, minprob, mean_imp)], by = "sample_id", all.x = TRUE)
  w[, sort_val := fifelse(is.finite(rec_val), rec_val, cross_protease)]
  setorder(w, sort_val, na.last = TRUE); w[, rank := seq_len(.N)]

  series <- rbindlist(list(
    w[is.finite(rec_val),        .(rank, value = rec_val,        series = rec_label)],
    w[is.finite(don_val),        .(rank, value = don_val,        series = don_label)],
    w[is.finite(cross_protease), .(rank, value = cross_protease, series = "Cross-protease imputed")],
    w[is.finite(mean_imp),       .(rank, value = mean_imp,       series = "Mean imputed")],
    w[is.finite(minprob),        .(rank, value = minprob,        series = "MinProb imputed")]))
  gl <- w[is.finite(cross_protease)]
  if (nrow(gl)) { gl[, ymin := pmin(cross_protease, minprob, mean_imp, don_val, na.rm = TRUE)]
                  gl[, ymax := pmax(cross_protease, minprob, mean_imp, don_val, na.rm = TRUE)] }

  lv <- c(rec_label, don_label, "Cross-protease imputed", "Mean imputed", "MinProb imputed")
  col_map   <- setNames(c(rec_col, don_col, CONC_COL, MEAN_COL, MP_COL), lv)
  shape_map <- setNames(c(16, 16, 17, 15, 18), lv)     # 15 = square for Mean
  size_map  <- setNames(c(2.5, 2.5, 3, 3, 3), lv)
  series[, series := factor(series, levels = lv)]

  p <- ggplot()
  if (nrow(gl)) p <- p + geom_segment(data = gl, aes(x = rank, xend = rank, y = ymin, yend = ymax),
                                      colour = "grey80", linewidth = 0.4)
  p + geom_point(data = series, aes(rank, value, colour = series, shape = series, size = series)) +
    scale_colour_manual(values = col_map) + scale_shape_manual(values = shape_map) +
    scale_size_manual(values = size_map) +
    labs(title = sprintf("%s (%s chosen, tau = %.2f)", gene_name, win, tau),
         x = sprintf("sample (ranked by %s abundance)", win), y = expression(log[2]~abundance)) +
    theme_pub(13)
}

# ---- verification vs full advisor -------------------------------------------
cat("\n--- Channel-selection verification ---\n")
va <- tryCatch({
  # optional cross-check against a full advisor run on Mix30 (protease-quant-advisor/run_all_conditions.sh)
  res <- readRDS(Sys.getenv("ADVISOR_STATE", file.path(ADVISOR_DIR, "results", "mix30", "advisor_state.rds")))
  unique(as.data.table(res$final_matrix)[gene %in% names(FOCUS),
    .(gene, full_advisor = recommended_channel, confidence = recommendation_confidence)])
}, error = function(e) NULL)
for (g in names(FOCUS)) {
  s <- shared[pg_canonical == FOCUS[[g]]]
  fa <- if (!is.null(va)) va[gene == g] else NULL
  cat(sprintf("  %-6s simplified: winner=%s (tryp_prec=%.2f hta_prec=%.2f, tau=%.2f) | full advisor: %s\n",
              g, s$winner, s$trypsin, s$hta, s$kendall_tau,
              if (!is.null(fa) && nrow(fa)) sprintf("%s (%s)", fa$full_advisor, fa$confidence) else "n/a"))
}

# ---- render -----------------------------------------------------------------
cat("\n--- Figures ---\n")
for (g in names(FOCUS)) {
  p <- make_profile(FOCUS[[g]], g)
  fp <- file.path(OUT, sprintf("fig_profile_%s_mean.png", g))
  ggsave(fp, p, width = 8, height = 5.5, dpi = 600, bg = "white")
  ggsave(sub("\\.png$", ".pdf", fp), p, width = 8, height = 5.5, bg = "white")
  cat(sprintf("  Saved: %s (+ .pdf)\n", basename(fp)))
}
cat("Done.\n")
