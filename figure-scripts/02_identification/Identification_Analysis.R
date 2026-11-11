# --------------------------------------------------------------------------------------------------
# Two-Protease-Paper | Identification analysis across the discovery searches
# Manuscript: Fig 2A-B (identification matrices, completeness, peptides per protein). Its overlap
# heatmaps and per-sample R^2 scatters are the analyses behind the Fig 1B overlap counts and Fig S5.
# Authors: Simion Kreimer (original); Mario M. Alba (manuscript group labels and comparison plan).
# Inputs: TWO_PROTEASE_DATA/{Krakatoa_Library_Search,Trypsin_Library_Search}/report.{pr,pg}_matrix.tsv
# Outputs: TWO_PROTEASE_OUT/02_identification/DIANN_Compare_Report_<timestamp>/
# Run: Rscript figure-scripts/02_identification/Identification_Analysis.R (prompts take defaults).
# Packages: data.table, ggplot2, grid, patchwork, cowplot.
# Released as used for the manuscript: only file paths, line endings and this header were changed;
# the analysis code is untouched. Paths come from environment variables, see ../../README.md.
# --------------------------------------------------------------------------------------------------

# ============================================================
# DIA-NN multi-search comparison (up to 10 searches)
# - Peptide + Protein level
# - Clean outputs (few useful plots/tables)
# - Sample grouping: auto (path-aware) + scheme-learning for Unknowns (bulk)
#
# FIXES IN THIS REVISION (per your notes):
# 1) Peptides-per-protein plot: ALL searches share the SAME y-axis (no free_y).
# 2) Scatterplots: were taking forever due to insane binning; now fast + visible:
#    - reasonable bin count
#    - log-scaled density
#    - optional ID sampling for speed (keeps plots readable)
#    - legend removed (only label box remains: equation, N, R²)
# 3) Interface: YOU select which search(es) to generate scatters + R² violin for
#    (no longer forced to Union/All Searches).
# 4) Overlap heatmaps: proteins on ONE side of diagonal, peptides on the OTHER,
#    with y-axis reversed so the triangle meaning is visually consistent.
# 5) Bugfix: scatter title used undefined variable; fixed.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(patchwork)
})

# ==============================
# USER SETTINGS
# ==============================
pg_file <- "report.pg_matrix.tsv"
pr_file <- "report.pr_matrix.tsv"

peptide_collapse_method <- "sum"  # "sum" or "mean"
zero_as_na <- TRUE

max_searches <- 100
cell_ratio <- 0.70  # <— shorter cells (wider relative to height)

min_points_for_cor <- 50
ALL_SEARCHES_LABEL <- "All Searches"

enable_group_plots <- TRUE
top_group_pairs_scatter <- 12

# Scatter/violin speed controls (IMPORTANT for huge matrices)
max_ids_for_scatter_calc <- 500000   # if > this many IDs, sample IDs before building scatter points
max_ids_for_r2_calc      <- 500000   # if > this many IDs, sample IDs before computing R² per SampleKey

max_scatter_points_plot  <- 250000   # plot-only cap (downsample points after collecting)
use_log2_for_cor <- TRUE

# regression_method:
#   "lm"         = ordinary least squares y~x
#   "orthogonal" = orthogonal regression (total least squares)
regression_method <- "orthogonal"

PUB_DPI <- 600
pub_theme <- function(base_size = 16, strip_size = 16) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", color = "black", size = base_size + 2),
      axis.title = element_text(face = "bold", color = "black", size = base_size + 1),
      axis.text  = element_text(color = "black", face = "plain", size = base_size),
      legend.title = element_text(face = "bold", color = "black", size = base_size),
      legend.text  = element_text(color = "black", size = base_size - 1),
      strip.text   = element_text(color = "black", face = "bold", size = strip_size),
      axis.text.x  = element_text(angle = 45, hjust = 1),
      panel.grid   = element_blank(),
      strip.background = element_blank()
    )
}


# ==============================
# FOLDER PICKER + OUTPUT
# ==============================
pick_parent_folder <- function() {
  if (interactive() && .Platform$OS.type == "windows") {
    p <- utils::choose.dir(caption = "Select parent folder containing DIA-NN search result folders")
    if (is.na(p) || !nzchar(p)) stop("No folder selected.")
    return(p)
  }
  stop("Folder picker is Windows interactive only. Set base_dir manually.")
}

# base_dir <- pick_parent_folder()                      # interactive picker (disabled)
# Paths come from environment variables (see the repository README): TWO_PROTEASE_DATA = the DIA-NN
# search results (MassIVE MSV000103457), TWO_PROTEASE_OUT = where reports are written.
base_dir <- Sys.getenv("TWO_PROTEASE_DATA", "data")   # parent folder of the DIA-NN search folders

# ---- Paste/Rscript-safe prompts (2-protease): take all default answers --------
# Every readline() returns "" so each prompt uses its DEFAULT branch and can never
# consume a following pasted line (defaults: folder short-names, auto group map,
# auto-pick top group pairs, "All Searches"). Set ASSUME_DEFAULTS <- FALSE to answer
# prompts yourself when running via source() in an interactive console.
ASSUME_DEFAULTS <- TRUE
if (ASSUME_DEFAULTS) readline <- function(prompt = "") { if (nzchar(prompt)) cat(prompt); "" }

timestamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
out_dir   <- file.path(Sys.getenv("TWO_PROTEASE_OUT", "output"), "02_identification",
                       paste0("DIANN_Compare_Report_", timestamp))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================
# LOGGING
# ==============================
log_file <- file.path(out_dir, "run_log.txt")
sink(log_file, split = TRUE)
cat("DIA-NN compare run started: ", format(Sys.time()), "\n")
cat("Base dir: ", base_dir, "\n")
cat("Output dir: ", out_dir, "\n\n")

on.exit({
  cat("\nRun finished: ", format(Sys.time()), "\n")
  sink()
}, add = TRUE)

# ==============================
# EXPORT + THEME
# ==============================
save_pub <- function(plot_obj, out_png, width, height, dpi = PUB_DPI) {
  ggsave(out_png, plot_obj, width = width, height = height, dpi = dpi, units = "in")
  out_pdf <- sub("\\.png$", ".pdf", out_png, ignore.case = TRUE)
  ggsave(out_pdf, plot_obj, width = width, height = height, device = cairo_pdf, units = "in")
}


# ==============================
# COLOR HELPERS
# ==============================
# red -> green
rg_fill_scale <- function(limits = c(0, 100), na.value = "white") {
  scale_fill_gradient(low = "red2", high = "green3", limits = limits, na.value = na.value)
}

# Peptides-per-protein: 1 peptide = red; 2+ peptides shift toward blue
ppp_palette <- function() {
  c(
    "1"    = "red2",
    "2"    = "#8FBCE6",
    "3"    = "#6EA6DC",
    "4"    = "#4D90D2",
    "5"    = "#2C7AC8",
    "6-10" = "#1764B7",
    ">10"  = "#084594"
  )
}

# Completeness: low completeness red -> high completeness green
comp_palette <- function() {
  lv <- c("0-20%", "20-50%", "50-80%", "80-95%", "95-100%", "100%")
  cols <- colorRampPalette(c("red2", "gold", "green3"))(length(lv))
  setNames(cols, lv)
}
shared_fill <- rg_fill_scale(limits = c(0, 100), na.value = "white")
# Scatter density binning (FAST + visible)
scatter_bins <- 500                 # reasonable; do NOT set huge
scatter_fill_low <- "white"
scatter_fill_high <- "dodgerblue4"
scatter_density_trans <- "log10"    # makes low-density visible
# ==============================
# HELPERS
# ==============================
stop_if_missing <- function(path) if (!file.exists(path)) stop("Missing file: ", path)

build_precursor_matrix_from_pr <- function(pr_dt,
                                           zero_as_na = TRUE,
                                           diag_tag = "pr_precursor",
                                           diag_dir = NULL) {
  pr_dt <- as.data.table(pr_dt)
  
  # Prefer Precursor.Id if available (already unique precursor key)
  if ("Precursor.Id" %in% names(pr_dt)) {
    id_col <- "Precursor.Id"
    pr_dt[[id_col]] <- sanitize_ids(pr_dt[[id_col]])
  } else {
    # Otherwise build precursor key from Modified.Sequence + PrecursorCharge
    if (!("Modified.Sequence" %in% names(pr_dt))) stop("pr_matrix missing Modified.Sequence")
    if (!("PrecursorCharge" %in% names(pr_dt)) && !("Precursor.Charge" %in% names(pr_dt))) {
      stop("pr_matrix missing PrecursorCharge (or Precursor.Charge)")
    }
    ch_col <- if ("PrecursorCharge" %in% names(pr_dt)) "PrecursorCharge" else "Precursor.Charge"
    
    pr_dt[["Modified.Sequence"]] <- sanitize_ids(pr_dt[["Modified.Sequence"]])
    pr_dt[[ch_col]] <- suppressWarnings(as.integer(pr_dt[[ch_col]]))
    
    pr_dt[, PrecursorKey := paste0(Modified.Sequence, pr_dt[[ch_col]])]
    id_col <- "PrecursorKey"
  }
  
  pr_dt <- pr_dt[!is.na(pr_dt[[id_col]]) & nzchar(pr_dt[[id_col]])]
  
  num_cols <- numeric_cols_only(pr_dt, exclude = c(id_col, "Modified.Sequence", "PrecursorKey"),
                                diag_tag = diag_tag, diag_dir = diag_dir)
  if (length(num_cols) == 0) stop("No numeric/sample-like columns detected in pr_matrix.")
  
  for (cname in num_cols) pr_dt[[cname]] <- suppressWarnings(as.numeric(pr_dt[[cname]]))
  if (zero_as_na) for (cname in num_cols) pr_dt[[cname]][pr_dt[[cname]] == 0] <- NA_real_
  
  # Use your existing matrix_from_dt helper
  m <- matrix_from_dt(pr_dt, id_col = id_col, num_cols = num_cols, zero_as_na = FALSE)
  # matrix_from_dt already turned 0->NA if you asked; we already did that above
  m
}

make_overlap_heatmap_plot_nolegend <- function(df,
                                               fill_limits = c(0, 100),
                                               text_size = 4.0,
                                               cell_ratio = 0.45) {
  df <- as.data.table(df)
  if (nrow(df) == 0) return(ggplot() + theme_void())
  
  ggplot(df, aes(x = Col, y = Row, fill = Value)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(
      aes(label = Label),
      size = text_size,
      fontface = "bold",
      color = "black",
      lineheight = 0.88
    ) +
    rg_fill_scale(limits = fill_limits, na.value = "white") +
    coord_fixed(ratio = cell_ratio, clip = "off") +
    labs(x = NULL, y = NULL, fill = "% overlap") +
    pub_theme(base_size = 13, strip_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.margin = margin(2, 2, 2, 2),
      legend.position = "none"
    )
}

make_overlap_heatmap_plot_from_df <- function(df, text_size = 4.0, fill_limits = c(0, 100), cell_ratio = 0.38) {
  df <- as.data.table(df)
  
  ggplot(df, aes(x = Col, y = Row, fill = Value)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = Label), size = text_size, fontface = "bold", color = "black", lineheight = 0.92) +
    rg_fill_scale(limits = fill_limits, na.value = "white") +
    coord_fixed(ratio = cell_ratio, clip = "off") +
    labs(x = NULL, y = NULL, fill = "% overlap") +
    pub_theme(base_size = 13, strip_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.margin = margin(2, 2, 2, 2)
    )
}

plot_two_heatmaps_side_by_side_legend_over_left <- function(df_left, df_right, out_png,
                                                            width = 16.5, height = 5.0,
                                                            text_size = 4.0, cell_ratio = 0.38,
                                                            left_height_rel = 0.88) {
  # Build plots FROM DF (so no duplicate scales/coord)
  p_left  <- make_overlap_heatmap_plot_from_df(df_left,  text_size = text_size, cell_ratio = cell_ratio)
  p_right <- make_overlap_heatmap_plot_from_df(df_right, text_size = text_size, cell_ratio = cell_ratio)
  
  # Legend: take from LEFT, horizontal, then remove legends from both
  p_left_leg <- p_left + theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal"
  )
  p_left_noleg  <- p_left  + theme(legend.position = "none")
  p_right_noleg <- p_right + theme(legend.position = "none")
  
  # Extract legend using cowplot (this is OK; no plot_spacer needed from cowplot)
  # If you don't have cowplot loaded, add: library(cowplot)
  leg <- cowplot::get_legend(p_left_leg)
  
  # Make left/right widths proportional to number of columns (more cells => wider)
  ncol_left  <- length(levels(as.factor(df_left$Col)))
  ncol_right <- length(levels(as.factor(df_right$Col)))
  w_left  <- max(1, ncol_left)
  w_right <- max(1, ncol_right)
  
  # Left panel: legend above + shorter plot area
  left_panel <- (
    cowplot::ggdraw(leg) /
      p_left_noleg +
      patchwork::plot_layout(heights = c(0.18, left_height_rel))
  )
  
  # Right panel: full height plot
  combined <- left_panel | p_right_noleg
  combined <- combined + patchwork::plot_layout(widths = c(w_left, w_right))
  
  ggsave(out_png, combined, width = width, height = height, dpi = PUB_DPI, units = "in")
  ggsave(sub("\\.png$", ".pdf", out_png, ignore.case = TRUE),
         combined, width = width, height = height, device = cairo_pdf, units = "in")
  
  invisible(combined)
}


make_df_search_heatmap <- function(peptide_pres, protein_pres) {
  searches <- colnames(peptide_pres)
  if (!identical(searches, colnames(protein_pres))) {
    stop("Peptide/protein presence matrices have different search columns.")
  }
  
  pep <- pairwise_overlap_pct(peptide_pres)
  pro <- pairwise_overlap_pct(protein_pres)
  
  n <- length(searches)
  df <- CJ(Row = searches, Col = searches)
  df[, `:=`(Value = NA_real_, Label = "")]
  
  for (i in seq_len(n)) for (j in seq_len(n)) {
    r <- searches[i]; c <- searches[j]
    
    if (i == j) {
      df[Row == r & Col == c, `:=`(
        Value = NA_real_,
        Label = paste0(fmt_int(pro$diag_n[r]), " proteins\n", fmt_int(pep$diag_n[r]), " peptides")
      )]
    } else if (i < j) {
      # PROTEINS
      pct <- pro$pct[r, c]
      nn  <- pro$n[r, c]
      df[Row == r & Col == c, `:=`(
        Value = pct,
        Label = paste0(fmt_int(nn), " proteins")
      )]
    } else {
      # PEPTIDES
      pct <- pep$pct[r, c]
      nn  <- pep$n[r, c]
      df[Row == r & Col == c, `:=`(
        Value = pct,
        Label = paste0(fmt_int(nn), " peptides")
      )]
    }
  }
  
  df[, Row := factor(Row, levels = rev(searches))]
  df[, Col := factor(Col, levels = searches)]
  df
}

make_df_group_heatmap_from_union <- function(m_union_pep, m_union_prot, group_map) {
  gm <- as.data.table(group_map)
  groups <- sort(unique(gm$Group))
  if (length(groups) < 2) return(data.table())
  
  sets_from_presence <- function(pres_mat) {
    out <- setNames(vector("list", length(groups)), groups)
    diag_counts <- integer(length(groups)); names(diag_counts) <- groups
    for (g in groups) {
      smps <- gm[Group == g, Sample]
      smps <- intersect(smps, colnames(pres_mat))
      if (length(smps) == 0) {
        out[[g]] <- character(0); diag_counts[g] <- 0L; next
      }
      present <- rowSums(pres_mat[, smps, drop = FALSE]) > 0
      out[[g]] <- rownames(pres_mat)[present]
      diag_counts[g] <- sum(present)
    }
    list(sets = out, diag = diag_counts)
  }
  
  jaccard_from_sets <- function(sets) {
    g <- names(sets)
    n <- length(g)
    pct <- matrix(NA_real_, n, n, dimnames = list(g, g))
    nn  <- matrix(NA_integer_, n, n, dimnames = list(g, g))
    for (i in seq_len(n)) for (j in seq_len(n)) {
      a <- sets[[i]]; b <- sets[[j]]
      inter <- length(intersect(a, b))
      uni <- length(union(a, b))
      nn[i, j] <- inter
      pct[i, j] <- if (uni == 0) NA_real_ else 100 * inter / uni
    }
    list(pct = pct, n = nn)
  }
  
  pep_pres <- (!is.na(m_union_pep))
  pro_pres <- (!is.na(m_union_prot))
  rownames(pep_pres) <- rownames(m_union_pep); colnames(pep_pres) <- colnames(m_union_pep)
  rownames(pro_pres) <- rownames(m_union_prot); colnames(pro_pres) <- colnames(m_union_prot)
  
  pep_sets <- sets_from_presence(pep_pres)
  pro_sets <- sets_from_presence(pro_pres)
  
  pep <- jaccard_from_sets(pep_sets$sets)
  pro <- jaccard_from_sets(pro_sets$sets)
  
  df <- CJ(Row = groups, Col = groups)
  df[, `:=`(Value = NA_real_, Label = "")]
  
  nG <- length(groups)
  for (i in seq_len(nG)) for (j in seq_len(nG)) {
    r <- groups[i]; c <- groups[j]
    if (i == j) {
      df[Row == r & Col == c, `:=`(
        Value = NA_real_,
        Label = paste0(fmt_int(pro_sets$diag[r]), " proteins\n", fmt_int(pep_sets$diag[r]), " peptides")
      )]
    } else if (i < j) {
      pct <- pro$pct[r, c]; nn <- pro$n[r, c]
      df[Row == r & Col == c, `:=`(Value = pct, Label = paste0(fmt_int(nn), " proteins"))]
    } else {
      pct <- pep$pct[r, c]; nn <- pep$n[r, c]
      df[Row == r & Col == c, `:=`(Value = pct, Label = paste0(fmt_int(nn), " peptides"))]
    }
  }
  
  df[, Row := factor(Row, levels = rev(groups))]
  df[, Col := factor(Col, levels = groups)]
  df
}

add_stack_percent_labels <- function(dt, group_cols = c("Level","Search","Group"),
                                     min_pct = 6, min_count = 50) {
  dt <- copy(as.data.table(dt))
  dt[, Total := sum(Count, na.rm = TRUE), by = group_cols]
  dt[, Percent := ifelse(Total > 0, 100 * Count / Total, NA_real_)]
  
  # Only show labels if they are readable:
  # - at least min_pct of the bar AND
  # - at least min_count IDs in that segment
  dt[, LabelPct := ifelse(is.finite(Percent) & Percent >= min_pct & Count >= min_count,
                          sprintf("%.0f%%", Percent), "")]
  dt[]
}

sanitize_ids <- function(ids) {
  ids <- as.character(ids)
  ids <- trimws(ids)
  ids[ids %in% c("", "NA", "NaN")] <- NA_character_
  ids
}

discover_search_folders <- function(parent_dir, pg_file, pr_file) {
  subs <- list.dirs(parent_dir, full.names = TRUE, recursive = FALSE)
  ok <- vapply(subs, function(d) {
    file.exists(file.path(d, pg_file)) && file.exists(file.path(d, pr_file))
  }, logical(1))
  found <- subs[ok]
  names(found) <- basename(found)
  found <- found[order(names(found))]
  found
}

make_unique_names <- function(x) {
  nm <- names(x)
  if (anyDuplicated(nm)) names(x) <- make.unique(nm, sep = "_")
  x
}

prompt_search_short_names <- function(search_paths_named) {
  cat("\nShort names for search strategies (press Enter to accept default):\n")
  out <- search_paths_named
  new_names <- character(length(out))
  for (i in seq_along(out)) {
    default <- names(out)[i]
    ans <- readline(sprintf("  %s  -> short name [%s]: ", default, default))
    ans <- trimws(ans)
    if (!nzchar(ans)) ans <- default
    new_names[i] <- ans
  }
  names(out) <- make.unique(new_names, sep = "_")
  out
}

prompt_group_renames <- function(group_map) {
  gm <- as.data.table(group_map)
  groups <- sort(unique(gm$Group))
  if (length(groups) == 0) return(gm)
  
  cat("\nRename groups? (press Enter to keep the same name)\n")
  ren <- data.table(Old = groups, New = groups)
  for (i in seq_along(groups)) {
    g <- groups[i]
    ans <- readline(sprintf("  Group '%s' -> new name [%s]: ", g, g))
    ans <- trimws(ans)
    if (!nzchar(ans)) ans <- g
    ren$New[i] <- ans
  }
  ren[, New := make.unique(New, sep = "_")]
  gm <- merge(gm, ren, by.x = "Group", by.y = "Old", all.x = TRUE)
  gm[, Group := New][, New := NULL]
  gm[]
}

# ==============================
# Hardcoded manuscript group labels (2-protease paper) -- PI decision 2026-06
# ------------------------------------------------------------------------------
# Replaces the interactive prompt_group_renames() so every figure carries the
# manuscript nomenclature instead of acquisition-date tokens:
#     Try_<date>  -> "Trypsin 1" / "Trypsin 2"   (earliest date = set 1)
#     Kra_<date>  -> "Krakatoa 1" / "Krakatoa 2"
#     Mix_20min   -> "Mix 20min"
#     Mix_30min   -> "Mix 30min"
# Set rank is derived from the sorted distinct acquisition dates, so it is robust
# if the exact dates change (it does not hardcode 250807/251117). Unknown/other
# labels pass through unchanged. classify_group_2p() keys batch-matching on the
# trailing token, which agrees across proteases in BOTH schemes, so the
# hardcoded comparison plan is unaffected by this relabel.
# ==============================
apply_manuscript_group_labels <- function(group_map) {
  gm <- as.data.table(group_map)
  if (!"Group" %in% names(gm) || nrow(gm) == 0) return(gm)

  prot_dates <- unique(sub("^(Try|Kra)_", "",
                           gm[grepl("^(Try|Kra)_[0-9]{6,}$", Group), Group]))
  prot_dates <- sort(prot_dates)
  date_rank  <- stats::setNames(seq_along(prot_dates), prot_dates)

  relabel_one <- function(g) {
    if (grepl("^Try_[0-9]{6,}$", g)) {
      d <- sub("^Try_", "", g); return(paste0("Trypsin ",  date_rank[[d]]))
    }
    if (grepl("^Kra_[0-9]{6,}$", g)) {
      d <- sub("^Kra_", "", g); return(paste0("Krakatoa ", date_rank[[d]]))
    }
    if (grepl("^Mix_20min$", g, ignore.case = TRUE)) return("Mix 20min")
    if (grepl("^Mix_30min$", g, ignore.case = TRUE)) return("Mix 30min")
    g
  }

  gm[, Group := vapply(as.character(Group), relabel_one, character(1))]
  cat("\n[group-labels] Applied manuscript labels:\n")
  print(unique(gm[, .(Group)])[order(Group)])
  gm[]
}

make_unique_samples <- function(samples, level_name, search_name, out_dir) {
  orig <- as.character(samples)
  uniq <- make.unique(orig, sep = "__dup")
  if (any(orig != uniq)) {
    fwrite(
      data.table(Level = level_name, Search = search_name, OriginalColumn = orig, UniqueUsed = uniq),
      file.path(out_dir, paste0("sample_name_mapping_", level_name, "_", search_name, ".csv"))
    )
  }
  uniq
}

sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

ppp_levels <- c(">10","6-10","5","4","3","2","1")

# Build percent labels per bar segment
add_ppp_percent_labels <- function(dt_bins) {
  dt <- copy(as.data.table(dt_bins))
  if (nrow(dt) == 0) return(dt)
  dt[, Total := sum(Count, na.rm = TRUE), by = .(Search, Group)]
  dt[, Percent := ifelse(Total > 0, 100 * Count / Total, NA_real_)]
  dt[, LabelPct := ifelse(is.finite(Percent) & Percent >= 3, sprintf("%.0f%%", Percent), "")]
  dt[]
}


heatmap_theme <- function(base = 13) {
  pub_theme(base_size = base, strip_size = base + 2) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.position = "bottom"
    )
}
build_group_overlap_heatmap_df <- function(m_union_pep, m_union_prot, group_map) {
  gm <- as.data.table(group_map)
  groups <- sort(unique(gm$Group))
  if (length(groups) < 2) return(data.table())
  
  sets_from_presence <- function(pres_mat) {
    out <- setNames(vector("list", length(groups)), groups)
    diag_counts <- integer(length(groups)); names(diag_counts) <- groups
    
    for (g in groups) {
      smps <- gm[Group == g, Sample]
      smps <- intersect(smps, colnames(pres_mat))
      if (length(smps) == 0) {
        out[[g]] <- character(0); diag_counts[g] <- 0L; next
      }
      present <- rowSums(pres_mat[, smps, drop = FALSE]) > 0
      out[[g]] <- rownames(pres_mat)[present]
      diag_counts[g] <- sum(present)
    }
    list(sets = out, diag = diag_counts)
  }
  
  jaccard_from_sets <- function(sets) {
    g <- names(sets)
    n <- length(g)
    pct <- matrix(NA_real_, n, n, dimnames = list(g, g))
    nn  <- matrix(NA_integer_, n, n, dimnames = list(g, g))
    for (i in seq_len(n)) for (j in seq_len(n)) {
      a <- sets[[i]]; b <- sets[[j]]
      inter <- length(intersect(a, b))
      uni <- length(union(a, b))
      nn[i, j] <- inter
      pct[i, j] <- if (uni == 0) NA_real_ else 100 * inter / uni
    }
    list(pct = pct, n = nn)
  }
  
  pep_pres <- (!is.na(m_union_pep))
  pro_pres <- (!is.na(m_union_prot))
  rownames(pep_pres) <- rownames(m_union_pep); colnames(pep_pres) <- colnames(m_union_pep)
  rownames(pro_pres) <- rownames(m_union_prot); colnames(pro_pres) <- colnames(m_union_prot)
  
  pep_sets <- sets_from_presence(pep_pres)
  pro_sets <- sets_from_presence(pro_pres)
  
  pep <- jaccard_from_sets(pep_sets$sets)
  pro <- jaccard_from_sets(pro_sets$sets)
  
  df <- CJ(Row = groups, Col = groups)
  df[, `:=`(Value = NA_real_, Label = "")]
  
  nG <- length(groups)
  for (i in seq_len(nG)) for (j in seq_len(nG)) {
    r <- groups[i]; c <- groups[j]
    if (i == j) {
      df[Row == r & Col == c, `:=`(
        Value = NA_real_,
        Label = paste0(fmt_int(pro_sets$diag[r]), " proteins\n", fmt_int(pep_sets$diag[r]), " peptides")
      )]
    } else if (i < j) {
      pct <- pro$pct[r, c]; nn <- pro$n[r, c]
      df[Row == r & Col == c, `:=`(Value = pct, Label = paste0(fmt_int(nn), " proteins"))]
    } else {
      pct <- pep$pct[r, c]; nn <- pep$n[r, c]
      df[Row == r & Col == c, `:=`(Value = pct, Label = paste0(fmt_int(nn), " peptides"))]
    }
  }
  
  df[, Row := factor(Row, levels = rev(groups))]
  df[, Col := factor(Col, levels = groups)]
  df[]
}

build_search_overlap_heatmap_df <- function(peptide_pres, protein_pres) {
  searches <- colnames(peptide_pres)
  if (!identical(searches, colnames(protein_pres))) {
    stop("Peptide/protein presence matrices have different search columns.")
  }
  
  pep <- pairwise_overlap_pct(peptide_pres)
  pro <- pairwise_overlap_pct(protein_pres)
  
  n <- length(searches)
  df <- CJ(Row = searches, Col = searches)
  df[, `:=`(Value = NA_real_, Label = "")]
  
  for (i in seq_len(n)) for (j in seq_len(n)) {
    r <- searches[i]; c <- searches[j]
    
    if (i == j) {
      df[Row == r & Col == c, `:=`(
        Value = NA_real_,
        Label = paste0(fmt_int(pro$diag_n[r]), " proteins\n", fmt_int(pep$diag_n[r]), " peptides")
      )]
    } else if (i < j) {
      pct <- pro$pct[r, c]; nn <- pro$n[r, c]
      df[Row == r & Col == c, `:=`(
        Value = pct,
        Label = paste0(fmt_int(nn), " proteins")
      )]
    } else {
      pct <- pep$pct[r, c]; nn <- pep$n[r, c]
      df[Row == r & Col == c, `:=`(
        Value = pct,
        Label = paste0(fmt_int(nn), " peptides")
      )]
    }
  }
  
  # keep your visual convention (upper triangle looks upper)
  df[, Row := factor(Row, levels = rev(searches))]
  df[, Col := factor(Col, levels = searches)]
  df[]
}

make_overlap_heatmap_plot <- function(df,
                                      fill_limits = c(0, 100),
                                      text_size = 4.2,
                                      cell_ratio = 0.45,      # < 1 => wider-than-tall cells
                                      show_legend = TRUE) {
  df <- as.data.table(df)
  if (nrow(df) == 0) return(ggplot() + theme_void())
  
  p <- ggplot(df, aes(x = Col, y = Row, fill = Value)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(
      aes(label = Label),
      size = text_size,
      fontface = "bold",
      color = "black",
      lineheight = 0.88
    ) +
    rg_fill_scale(limits = fill_limits, na.value = "white") +
    coord_fixed(ratio = cell_ratio, clip = "off") +
    labs(x = NULL, y = NULL, fill = "% overlap") +
    pub_theme(base_size = 13, strip_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.margin = margin(2, 2, 2, 2),
      legend.position = if (show_legend) "top" else "none",
      legend.direction = "horizontal",
      legend.title = element_text(face = "bold"),
      legend.key.height = unit(0.35, "cm"),
      legend.key.width  = unit(1.4, "cm")
    )
  
  p
}


# ==============================
# Two-heatmap panel helper
# ==============================
plot_two_heatmaps_side_by_side_legend_left <- function(df_left, df_right, out_png,
                                                       width = 15.5, height = 5.2,
                                                       cell_ratio = 0.45,
                                                       text_size_left = 4.4,
                                                       text_size_right = 3.8) {
  # Determine widths by number of columns (real cells only, NO padding)
  nC_left  <- length(levels(df_left$Col))
  nC_right <- length(levels(df_right$Col))
  widths <- c(nC_left, nC_right)
  
  # Build plots
  p_left  <- make_overlap_heatmap_plot(df_left,  text_size = text_size_left,  cell_ratio = cell_ratio, show_legend = TRUE)
  p_right <- make_overlap_heatmap_plot(df_right, text_size = text_size_right, cell_ratio = cell_ratio, show_legend = FALSE)
  
  # Layout:
  # Row 1: legend area over LEFT only
  # Row 2: the two plots side-by-side
  p_combined <-
    (
      patchwork::guide_area() | patchwork::plot_spacer()
    ) /
    (
      p_left | p_right
    ) +
    patchwork::plot_layout(
      guides = "collect",
      widths = widths,
      heights = c(0.12, 1)   # legend row is short
    ) &
    theme(legend.position = "top", legend.direction = "horizontal")
  
  ggsave(out_png, p_combined, width = width, height = height, dpi = PUB_DPI, units = "in")
  ggsave(sub("\\.png$", ".pdf", out_png, ignore.case = TRUE),
         p_combined, width = width, height = height, device = cairo_pdf, units = "in")
  
  p_combined
}



# ==============================
# HARD EXCLUDE NON-SAMPLE COLUMNS
# ==============================
normalize_colname <- function(x) tolower(gsub("[^a-z0-9]+", "", x))

is_non_sample_col <- function(colname) {
  raw <- tolower(trimws(colname))
  nrm <- normalize_colname(colname)
  
  if (nrm %in% c(
    "proteotypic",
    "nsequences", "nsequence", "nseq", "nseqs",
    "precursorcharge",
    "qvalue", "pgqvalue",
    "pep", "pgpep",
    "libqvalue",
    "decoy", "decoyevidence",
    "rt", "rtstart", "rtstop",
    "ms1area", "ms1profilecorr"
  )) return(TRUE)
  
  if (grepl("^precursor\\.?\\s*charge$", raw)) return(TRUE)
  if (grepl("^proteotypic$", raw)) return(TRUE)
  if (grepl("^n\\.?\\s*sequences$", raw)) return(TRUE)
  if (grepl("^n\\s*seq(uences)?$", raw)) return(TRUE)
  
  FALSE
}

looks_like_sample_colname <- function(colname) {
  grepl("\\.(raw|d|mzml)$", colname, ignore.case = TRUE) ||
    grepl("^[A-Za-z]:\\\\", colname) || grepl("^[A-Za-z]:/", colname)
}

numeric_cols_only <- function(dt, exclude = character(), diag_tag = "unknown", diag_dir = NULL) {
  cols_all <- setdiff(names(dt), exclude)
  
  dropped_non_sample <- cols_all[vapply(cols_all, is_non_sample_col, logical(1))]
  cols <- setdiff(cols_all, dropped_non_sample)
  
  keep_numeric <- character()
  keep_sample_like_but_non_numeric <- character()
  dropped_non_numeric <- character()
  
  for (cname in cols) {
    v <- suppressWarnings(as.numeric(dt[[cname]]))
    nn <- sum(!is.na(v))
    if (nn > 0) {
      keep_numeric <- c(keep_numeric, cname)
    } else {
      if (looks_like_sample_colname(cname)) {
        keep_sample_like_but_non_numeric <- c(keep_sample_like_but_non_numeric, cname)
      } else {
        dropped_non_numeric <- c(dropped_non_numeric, cname)
      }
    }
  }
  
  if (!is.null(diag_dir)) {
    fwrite(data.table(Column = dropped_non_sample),
           file.path(diag_dir, paste0("excluded_non_sample_cols_", diag_tag, ".csv")))
    fwrite(data.table(Column = dropped_non_numeric),
           file.path(diag_dir, paste0("excluded_non_numeric_cols_", diag_tag, ".csv")))
    fwrite(data.table(Column = keep_sample_like_but_non_numeric),
           file.path(diag_dir, paste0("kept_sample_like_but_non_numeric_", diag_tag, ".csv")))
  }
  
  unique(c(keep_numeric, keep_sample_like_but_non_numeric))
}

read_diann_matrix <- function(path, id_col_guess = NULL, diag_tag = "matrix", diag_dir = NULL) {
  dt <- fread(path, sep = "\t", header = TRUE, data.table = TRUE, showProgress = FALSE)
  if (nrow(dt) == 0) stop("Empty file: ", path)
  
  if (!is.null(id_col_guess) && id_col_guess %in% names(dt)) {
    id_col <- id_col_guess
  } else {
    candidates <- c("Protein.Group", "Genes", "ModifiedPeptideSequence", "Modified.Sequence",
                    "Stripped.Sequence", names(dt)[1])
    candidates <- candidates[candidates %in% names(dt)]
    if (length(candidates) == 0) stop("Could not infer ID column in: ", path)
    id_col <- candidates[1]
  }
  
  num_cols <- numeric_cols_only(dt, exclude = id_col, diag_tag = diag_tag, diag_dir = diag_dir)
  for (cname in num_cols) dt[[cname]] <- suppressWarnings(as.numeric(dt[[cname]]))
  list(dt = dt, id_col = id_col, num_cols = num_cols)
}

matrix_from_dt <- function(dt, id_col, num_cols, zero_as_na = TRUE) {
  ids <- sanitize_ids(dt[[id_col]])
  keep <- !is.na(ids) & nzchar(ids)
  if (!all(keep)) {
    dt <- dt[keep]
    ids <- ids[keep]
  }
  
  m <- as.matrix(dt[, ..num_cols])
  mode(m) <- "numeric"
  if (zero_as_na) m[m == 0] <- NA_real_
  
  if (anyDuplicated(ids)) {
    agg <- data.table(ID = ids)
    for (j in seq_along(num_cols)) agg[[num_cols[j]]] <- m[, j]
    agg <- agg[, lapply(.SD, function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)),
               by = ID, .SDcols = num_cols]
    ids2 <- agg[["ID"]]
    m2 <- as.matrix(agg[, ..num_cols])
    rownames(m2) <- ids2
    colnames(m2) <- num_cols
    return(m2)
  } else {
    rownames(m) <- ids
    colnames(m) <- num_cols
    return(m)
  }
}

collapse_to_peptides <- function(pr_dt,
                                 collapse_method = "sum",
                                 zero_as_na = TRUE,
                                 peptide_id_candidates = c(
                                   "ModifiedPeptideSequence",
                                   "Modified.Sequence",
                                   "Modified sequence",
                                   "Modified.Peptide.Sequence",
                                   "ModifiedPeptide.Sequence",
                                   "Stripped.Sequence"
                                 ),
                                 diag_tag = "pr",
                                 diag_dir = NULL) {
  pep_id <- peptide_id_candidates[peptide_id_candidates %in% names(pr_dt)]
  if (length(pep_id) == 0) stop("pr_matrix missing peptide sequence column.")
  id_col <- pep_id[1]
  
  pr_dt[[id_col]] <- sanitize_ids(pr_dt[[id_col]])
  pr_dt <- pr_dt[!is.na(pr_dt[[id_col]]) & nzchar(pr_dt[[id_col]])]
  
  num_cols <- numeric_cols_only(pr_dt, exclude = id_col, diag_tag = diag_tag, diag_dir = diag_dir)
  if (length(num_cols) == 0) stop("No numeric/sample-like columns detected in pr_matrix after excluding ID column: ", id_col)
  
  for (cname in num_cols) pr_dt[[cname]] <- suppressWarnings(as.numeric(pr_dt[[cname]]))
  if (zero_as_na) for (cname in num_cols) pr_dt[[cname]][pr_dt[[cname]] == 0] <- NA_real_
  
  if (collapse_method == "sum") {
    pep <- pr_dt[, lapply(.SD, function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)),
                 by = id_col, .SDcols = num_cols]
  } else if (collapse_method == "mean") {
    pep <- pr_dt[, lapply(.SD, function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)),
                 by = id_col, .SDcols = num_cols]
  } else stop("Unknown collapse_method: ", collapse_method)
  
  pep_ids <- sanitize_ids(pep[[id_col]])
  keep <- !is.na(pep_ids) & nzchar(pep_ids)
  pep <- pep[keep]
  pep_ids <- pep_ids[keep]
  
  pep[[id_col]] <- NULL
  m <- as.matrix(pep)
  mode(m) <- "numeric"
  rownames(m) <- pep_ids
  colnames(m) <- num_cols
  m
}

align_samples_union <- function(mats) {
  all_cols <- lapply(mats, colnames)
  all_union <- sort(unique(unlist(all_cols)))
  if (length(all_union) == 0) stop("No sample columns found.")
  
  mats2 <- lapply(mats, function(m) {
    missing <- setdiff(all_union, colnames(m))
    if (length(missing) > 0) {
      add <- matrix(NA_real_, nrow = nrow(m), ncol = length(missing))
      colnames(add) <- missing
      rownames(add) <- rownames(m)
      m <- cbind(m, add)
    }
    m[, all_union, drop = FALSE]
  })
  
  list(mats = mats2, samples = all_union)
}

align_ids_union <- function(mats_named) {
  all_ids <- sort(unique(unlist(lapply(mats_named, rownames))))
  mats2 <- lapply(mats_named, function(m) {
    miss <- setdiff(all_ids, rownames(m))
    if (length(miss) > 0) {
      add <- matrix(NA_real_, nrow = length(miss), ncol = ncol(m))
      rownames(add) <- miss
      colnames(add) <- colnames(m)
      m <- rbind(m, add)
    }
    m[all_ids, , drop = FALSE]
  })
  mats2
}

build_union_matrix <- function(mats_named) {
  mats_named <- align_ids_union(mats_named)
  sample_ref <- colnames(mats_named[[1]])
  
  sum_m <- matrix(0, nrow = nrow(mats_named[[1]]), ncol = length(sample_ref))
  cnt_m <- matrix(0, nrow = nrow(mats_named[[1]]), ncol = length(sample_ref))
  rownames(sum_m) <- rownames(mats_named[[1]])
  colnames(sum_m) <- sample_ref
  rownames(cnt_m) <- rownames(sum_m)
  colnames(cnt_m) <- sample_ref
  
  for (s in names(mats_named)) {
    m <- mats_named[[s]][, sample_ref, drop = FALSE]
    ok <- !is.na(m)
    sum_m[ok] <- sum_m[ok] + m[ok]
    cnt_m[ok] <- cnt_m[ok] + 1
  }
  
  out <- sum_m
  out[cnt_m == 0] <- NA_real_
  out[cnt_m > 0] <- sum_m[cnt_m > 0] / cnt_m[cnt_m > 0]
  out
}

build_union_presence <- function(mats_named) {
  mats_named <- align_ids_union(mats_named)
  sample_ref <- colnames(mats_named[[1]])
  out <- matrix(FALSE, nrow = nrow(mats_named[[1]]), ncol = length(sample_ref))
  rownames(out) <- rownames(mats_named[[1]])
  colnames(out) <- sample_ref
  for (s in names(mats_named)) {
    m <- mats_named[[s]][, sample_ref, drop = FALSE]
    out <- out | (!is.na(m))
  }
  out
}

id_presence <- function(mats_named) {
  ids_all <- sort(unique(unlist(lapply(mats_named, function(m) {
    rn <- sanitize_ids(rownames(m))
    rn <- rn[!is.na(rn) & nzchar(rn)]
    rn
  }))))
  pres <- matrix(FALSE, nrow = length(ids_all), ncol = length(mats_named))
  rownames(pres) <- ids_all
  colnames(pres) <- names(mats_named)
  for (s in names(mats_named)) {
    rn <- sanitize_ids(rownames(mats_named[[s]]))
    rn <- rn[!is.na(rn) & nzchar(rn)]
    rn <- intersect(rn, ids_all)
    if (length(rn) > 0) pres[rn, s] <- TRUE
  }
  pres
}

pairwise_overlap_pct <- function(pres) {
  s <- colnames(pres); n <- length(s)
  overlap_pct <- matrix(NA_real_, n, n, dimnames = list(s, s))
  overlap_n   <- matrix(NA_integer_, n, n, dimnames = list(s, s))
  diag_n      <- integer(n); names(diag_n) <- s
  
  for (i in seq_len(n)) diag_n[i] <- sum(pres[, i])
  for (i in seq_len(n)) for (j in seq_len(n)) {
    inter <- sum(pres[, i] & pres[, j])
    uni   <- sum(pres[, i] | pres[, j])
    overlap_n[i, j] <- inter
    overlap_pct[i, j] <- if (uni == 0) NA_real_ else 100 * inter / uni
  }
  list(pct = overlap_pct, n = overlap_n, diag_n = diag_n)
}

build_idcount_df_with_all_searches <- function(mats_named, level_name) {
  searches <- names(mats_named)
  samples <- colnames(mats_named[[1]])
  
  per_search <- rbindlist(lapply(searches, function(s) {
    m <- mats_named[[s]]
    data.table(Level = level_name, Sample = colnames(m), Search = s, IDs = colSums(!is.na(m)))
  }))
  
  all_ids <- sort(unique(unlist(lapply(mats_named, rownames))))
  all_counts <- integer(length(samples)); names(all_counts) <- samples
  
  for (k in seq_along(samples)) {
    samp <- samples[k]
    present_any <- rep(FALSE, length(all_ids)); names(present_any) <- all_ids
    for (s in searches) {
      m <- mats_named[[s]]
      idx <- rownames(m)[!is.na(m[, samp])]
      if (length(idx) > 0) present_any[idx] <- TRUE
    }
    all_counts[k] <- sum(present_any)
  }
  
  all_dt <- data.table(Level = level_name, Sample = samples, Search = ALL_SEARCHES_LABEL, IDs = all_counts)
  rbind(per_search, all_dt)
}

# ==============================
# GROUPING (path-aware) + scheme learning
# ==============================
get_basename_noext <- function(x) {
  x2 <- gsub("\\\\", "/", x)
  b <- basename(x2)
  b <- sub("\\.raw$", "", b, ignore.case = TRUE)
  b <- sub("\\.d$",   "", b, ignore.case = TRUE)
  b <- sub("\\.mzML$", "", b, ignore.case = TRUE)
  b
}

get_parent_folder <- function(x) {
  x2 <- gsub("\\\\", "/", x)
  d <- dirname(x2)
  if (is.na(d) || d == "." || d == "") return(NA_character_)
  basename(d)
}

extract_date_token <- function(x) {
  nm <- get_basename_noext(x)
  m <- regexpr("^\\d{6}(?=_)", nm, perl = TRUE)
  if (m[1] != -1) return(regmatches(nm, m))
  NA_character_
}

tokens_underscore <- function(x) {
  nm <- get_basename_noext(x)
  toks <- strsplit(nm, "_", fixed = TRUE)[[1]]
  toks <- toks[nzchar(toks)]
  toks
}

filter_sample_like_names <- function(sample_names) {
  keep <- !vapply(sample_names, is_non_sample_col, logical(1))
  sample_names[keep]
}

auto_guess_group <- function(sample_name) {
  parent <- get_parent_folder(sample_name)
  nm     <- get_basename_noext(sample_name)
  dtok   <- extract_date_token(sample_name)
  
  if (!is.na(parent) && grepl("20min", parent, ignore.case = TRUE) && grepl("mix", parent, ignore.case = TRUE))
    return("Mix_20min")
  if (grepl("_20min\\b", nm, ignore.case = TRUE) && grepl("_Mix_", nm, ignore.case = TRUE))
    return("Mix_20min")
  
  if (!is.na(parent) && grepl("30min", parent, ignore.case = TRUE) && grepl("mix", parent, ignore.case = TRUE))
    return("Mix_30min")
  if (grepl("_Mix_", nm, ignore.case = TRUE) && !grepl("_20min\\b", nm, ignore.case = TRUE))
    return("Mix_30min")
  
  if (!is.na(parent) && grepl("krakatoa", parent, ignore.case = TRUE)) {
    if (!is.na(dtok)) return(paste0("Kra_", dtok))
    return("Kra")
  }
  if (grepl("_Kra_", nm, ignore.case = TRUE)) {
    if (!is.na(dtok)) return(paste0("Kra_", dtok))
    return("Kra")
  }
  
  if (!is.na(parent) && grepl("trypsin", parent, ignore.case = TRUE)) {
    if (!is.na(dtok)) return(paste0("Try_", dtok))
    return("Try")
  }
  if (grepl("_Try_", nm, ignore.case = TRUE) || grepl("trypsin|tryptic", nm, ignore.case = TRUE)) {
    if (!is.na(dtok)) return(paste0("Try_", dtok))
    return("Try")
  }
  if (grepl("_BOB_", nm, ignore.case = TRUE)) {
    if (!is.na(dtok)) return(paste0("Try_", dtok))
    return("Try")
  }
  
  "Unknown"
}

build_group_map_auto <- function(samples) {
  data.table(
    Sample = samples,
    Group  = vapply(samples, auto_guess_group, character(1)),
    Source = "auto"
  )
}

learn_rules_from_examples <- function(unknown_samples, k = 3) {
  unknown_samples <- unique(unknown_samples)
  
  show_n <- min(length(unknown_samples), 15)
  cat("\nUnknown sample preview (first ", show_n, "):\n", sep = "")
  cat(paste0("  [", seq_len(show_n), "] ", unknown_samples[seq_len(show_n)]), sep = "\n")
  cat("\nWe'll learn rules from a few examples and apply to the rest.\n")
  
  labs <- data.table(Sample = character(), Group = character())
  for (i in seq_len(min(k, show_n))) {
    idx <- readline(sprintf("Pick an index to label (1-%d) [default %d]: ", show_n, i))
    idx <- trimws(idx)
    if (!nzchar(idx)) idx <- as.character(i)
    idx <- suppressWarnings(as.integer(idx))
    if (!is.finite(idx) || idx < 1 || idx > show_n) {
      cat("  Invalid index, skipping.\n")
      next
    }
    s <- unknown_samples[idx]
    g <- trimws(readline(sprintf("  Group name for '%s': ", s)))
    if (!nzchar(g)) {
      cat("  Empty group name, skipping.\n")
      next
    }
    labs <- rbind(labs, data.table(Sample = s, Group = g))
  }
  
  if (nrow(labs) == 0) return(data.table(Group = character(), Mode = character(), Pattern = character()))
  
  rules <- data.table(Group = character(), Mode = character(), Pattern = character())
  for (g in unique(labs$Group)) {
    ex <- labs[Group == g, Sample][1]
    toks <- tokens_underscore(ex)
    toks2 <- toks[!grepl("^\\d{6,}$", toks)]
    if (length(toks2) == 0) toks2 <- toks
    tok <- toks2[which.max(nchar(toks2))]
    rules <- rbind(rules, data.table(Group = g, Mode = "token", Pattern = tok))
  }
  
  rules
}

apply_rules_token <- function(group_map, rules) {
  if (nrow(rules) == 0) return(group_map)
  
  for (i in seq_len(nrow(rules))) {
    g <- rules$Group[i]
    tok <- rules$Pattern[i]
    
    unk <- group_map[Group == "Unknown", Sample]
    if (length(unk) == 0) break
    
    hit <- vapply(unk, function(s) {
      any(tolower(tokens_underscore(s)) == tolower(tok))
    }, logical(1))
    
    to_assign <- unk[hit]
    if (length(to_assign) > 0) {
      group_map[Sample %in% to_assign & Group == "Unknown",
                `:=`(Group = g, Source = paste0("learned:token:", tok))]
    }
  }
  group_map
}

assist_group_unknowns_learn_scheme <- function(group_map, k = 3) {
  unk <- group_map[Group == "Unknown", Sample]
  if (length(unk) == 0) return(group_map)
  
  cat("\n=== Unknown grouping: scheme learning ===\n")
  cat("Unknowns remaining: ", length(unk), "\n", sep = "")
  
  rules <- learn_rules_from_examples(unk, k = k)
  if (nrow(rules) > 0) {
    cat("\nLearned rules:\n")
    print(rules)
    group_map <- apply_rules_token(group_map, rules)
  }
  
  unk2 <- group_map[Group == "Unknown", Sample]
  if (length(unk2) > 0) {
    cat("\nAfter learned rules, Unknowns remaining: ", length(unk2), "\n", sep = "")
    fb <- trimws(readline("Assign all remaining Unknowns to a fallback group? (Enter name, or press Enter to leave Unknown): "))
    if (nzchar(fb)) group_map[Group == "Unknown", `:=`(Group = fb, Source = "fallback")]
  }
  
  group_map
}


# ==============================
# OPTIONAL MANUAL GROUP ASSIGNMENT TABLE
# ==============================
write_manual_group_template <- function(group_map, out_dir) {
  template_path <- file.path(out_dir, "sample_group_mapping_MANUAL_EDIT.csv")
  gm <- as.data.table(group_map)
  gm[, Instructions := "Edit Group values as needed, save the CSV, then press Enter in R."]
  fwrite(gm[, .(Sample, Group, Source, Instructions)], template_path)
  template_path
}

read_manual_group_template <- function(template_path, original_samples) {
  if (!file.exists(template_path)) stop("Manual group table was not found: ", template_path)
  gm <- fread(template_path, data.table = TRUE)
  required <- c("Sample", "Group")
  missing <- setdiff(required, names(gm))
  if (length(missing) > 0) stop("Manual group table is missing column(s): ", paste(missing, collapse = ", "))
  gm <- gm[, .(Sample = as.character(Sample), Group = as.character(Group))]
  gm[, `:=`(Sample = trimws(Sample), Group = trimws(Group))]
  gm <- gm[Sample %in% original_samples]
  gm <- gm[!is.na(Sample) & nzchar(Sample) & !is.na(Group) & nzchar(Group)]
  gm <- unique(gm, by = "Sample")
  missing_samples <- setdiff(original_samples, gm$Sample)
  if (length(missing_samples) > 0) {
    stop("Manual group table does not assign every sample. Missing: ", paste(head(missing_samples, 20), collapse = "; "),
         ifelse(length(missing_samples) > 20, " ...", ""))
  }
  gm[, Source := "manual_table"]
  gm[]
}

maybe_use_manual_group_table <- function(group_map, out_dir) {
  gm <- as.data.table(group_map)
  cat("\nCurrent sample-to-group assignment preview:\n")
  print(gm[, .N, by = Group][order(-N)])
  cat("\nA manual assignment template can be written if the auto-assignment is not correct.\n")
  ans <- trimws(readline("Use a manual CSV table to assign groups? [y/N]: "))
  if (!tolower(ans) %in% c("y", "yes")) return(gm)
  
  template_path <- write_manual_group_template(gm, out_dir)
  cat("\nWrote manual group table:\n  ", template_path, "\n", sep = "")
  cat("Open this CSV, edit the Group column, save it, then return here.\n")
  invisible(readline("Press Enter after saving the manual group table..."))
  
  gm2 <- read_manual_group_template(template_path, gm$Sample)
  cat("\nManual group counts:\n")
  print(gm2[, .N, by = Group][order(-N)])
  gm2[]
}

# ==============================
# SAMPLE KEY (biological sample) for matching across groups
# ==============================
extract_sample_key <- function(sample_path) {
  nm <- get_basename_noext(sample_path)
  
  # Normalize separators (just for matching)
  x <- gsub("[^A-Za-z0-9]+", "_", nm)
  
  # 1) Prefer HMN tokens like HMN849926A / HMN849931B (case-insensitive)
  # Allow optional trailing letter(s) A/B/etc.
  m <- regexpr("(?i)HMN\\d+[A-Z]?", x, perl = TRUE)
  if (m[1] != -1) {
    key <- toupper(regmatches(x, m))
    return(key)
  }
  
  # 2) If you also have BRH#### style keys, keep them as fallback
  m <- regexpr("(?i)BRH\\d+", x, perl = TRUE)
  if (m[1] != -1) {
    key <- toupper(regmatches(x, m))
    return(key)
  }
  
  # 3) Final fallback: full basename without extension
  nm
}


make_sample_meta <- function(samples, group_map) {
  dt <- data.table(Sample = samples)
  dt[, SampleKey := vapply(Sample, extract_sample_key, character(1))]
  dt <- merge(dt, group_map[, .(Sample, Group)], by = "Sample", all.x = TRUE)
  dt
}

plot_group_all_ids_combined <- function(id_df_pep, id_df_prot, group_map,
                                        out_png = NULL, return_plot = FALSE,
                                        alpha_groups = TRUE) {
  id_df_pep  <- as.data.table(id_df_pep)
  id_df_prot <- as.data.table(id_df_prot)
  group_map  <- as.data.table(group_map)
  
  id_df <- rbind(id_df_pep, id_df_prot, fill = TRUE)
  if (nrow(id_df) == 0) return(invisible(NULL))
  
  id_df <- merge(id_df, group_map[, .(Sample, Group)], by = "Sample", all.x = TRUE, allow.cartesian = TRUE)
  
  # ---- Group order (alphabetical to match your request)
  ord <- sort(unique(as.character(id_df$Group)))
  id_df[, Group := factor(as.character(Group), levels = ord)]
  
  # ---- Search order: strategies then Union at end
  search_levels <- c(setdiff(unique(id_df$Search), ALL_SEARCHES_LABEL), ALL_SEARCHES_LABEL)
  id_df[, Search := factor(Search, levels = search_levels)]
  
  # ---- Level order (Peptide on top, Protein on bottom)
  id_df[, Level := factor(Level, levels = c("Peptide", "Protein"))]
  
  # ---- vertical dividers between groups
  vlines <- data.table(xint = seq_along(ord) + 0.5)
  
  p <- ggplot(id_df, aes(x = Group, y = IDs, fill = Search)) +
    geom_violin(position = position_dodge(width = 0.85), trim = TRUE, scale = "width", alpha = 0.75) +
    # points (make easier to see)
    geom_point(
      aes(fill = Search, group = interaction(Group, Search)),
      position = position_jitterdodge(jitter.width = 0.10, dodge.width = 0.85),
      shape = 21, color = "black", stroke = 0.25,
      size = 1.7, alpha = 0.75
    ) +
    guides(color = "none") + 
    # median PER Group x Search
    stat_summary(
      aes(group = Search),
      fun = median,
      geom = "point",
      position = position_dodge(width = 0.85),
      shape = 23,
      fill = "white",
      color = "black",
      size = 2.7
    ) +
    # mean PER Group x Search
    stat_summary(
      aes(group = Search),
      fun = mean,
      geom = "point",
      position = position_dodge(width = 0.85),
      shape = 4,
      color = "black",
      stroke = 1.0,
      size = 2.8
    ) +
    
    # ---- Dividers
    geom_vline(
      data = vlines, aes(xintercept = xint),
      linewidth = 0.35, color = "black")+
    facet_grid(Level ~ ., scales = "free_y", switch = "y") +
    labs(title = NULL, x = "Group", y = NULL, fill = NULL, color = NULL) +  # y title removed
    pub_theme(base_size = 16, strip_size = 16) +
    theme(
      strip.placement = "outside",
      plot.margin = margin(6, 6, 6, 6)
    )
  
  if (isTRUE(return_plot)) return(p)
  
  # ---- Make it taller
  if (!is.null(out_png)) save_pub(p, out_png, width = 15.5, height = 8.2)
  
  invisible(p)
}


# ==============================
# Overlap heatmaps (combined triangles) - red->green
#  - VISUAL FIX: reverse y-axis so "upper triangle" is actually visually above diagonal
# ==============================
fmt_int <- function(x) formatC(x, format = "d", big.mark = ",")

plot_search_overlap_combined <- function(peptide_pres, protein_pres, out_png = NULL, return_plot = FALSE) {
  searches <- colnames(peptide_pres)
  if (!identical(searches, colnames(protein_pres))) {
    stop("Peptide/protein presence matrices have different search columns.")
  }
  
  pep <- pairwise_overlap_pct(peptide_pres)
  pro <- pairwise_overlap_pct(protein_pres)
  
  n <- length(searches)
  df <- CJ(Row = searches, Col = searches)
  df[, `:=`(Value = NA_real_, Label = "")]
  
  for (i in seq_len(n)) for (j in seq_len(n)) {
    r <- searches[i]; c <- searches[j]
    
    if (i == j) {
      df[Row == r & Col == c, `:=`(
        Value = NA_real_,
        Label = paste0(fmt_int(pro$diag_n[r]), " proteins\n", fmt_int(pep$diag_n[r]), " peptides")
      )]
    } else if (i < j) {
      pct <- pro$pct[r, c]
      nn  <- pro$n[r, c]
      df[Row == r & Col == c, `:=`(
        Value = pct,
        Label = paste0(fmt_int(nn), " proteins")
      )]
    } else {
      pct <- pep$pct[r, c]
      nn  <- pep$n[r, c]
      df[Row == r & Col == c, `:=`(
        Value = pct,
        Label = paste0(fmt_int(nn), " peptides")
      )]
    }
  }
  
  df[, Row := factor(Row, levels = rev(searches))]
  df[, Col := factor(Col, levels = searches)]
  
  p <- make_overlap_heatmap_plot(df, fill_limits = c(0,100), text_size = 3.4)
  
  if (isTRUE(return_plot)) return(p)
  if (!is.null(out_png)) save_pub(p, out_png, width = 10.5, height = 6.0)
  invisible(NULL)
  
}

plot_group_overlap_combined_from_union <- function(m_union_pep, m_union_prot, group_map,
                                                   out_png = NULL, return_plot = FALSE) {
  gm <- as.data.table(group_map)
  groups <- sort(unique(gm$Group))
  if (length(groups) < 2) return(invisible(NULL))
  
  sets_from_presence <- function(pres_mat) {
    out <- setNames(vector("list", length(groups)), groups)
    diag_counts <- integer(length(groups)); names(diag_counts) <- groups
    for (g in groups) {
      smps <- gm[Group == g, Sample]
      smps <- intersect(smps, colnames(pres_mat))
      if (length(smps) == 0) {
        out[[g]] <- character(0); diag_counts[g] <- 0L; next
      }
      present <- rowSums(pres_mat[, smps, drop = FALSE]) > 0
      out[[g]] <- rownames(pres_mat)[present]
      diag_counts[g] <- sum(present)
    }
    list(sets = out, diag = diag_counts)
  }
  
  jaccard_from_sets <- function(sets) {
    g <- names(sets)
    n <- length(g)
    pct <- matrix(NA_real_, n, n, dimnames = list(g, g))
    nn  <- matrix(NA_integer_, n, n, dimnames = list(g, g))
    for (i in seq_len(n)) for (j in seq_len(n)) {
      a <- sets[[i]]; b <- sets[[j]]
      inter <- length(intersect(a, b))
      uni <- length(union(a, b))
      nn[i, j] <- inter
      pct[i, j] <- if (uni == 0) NA_real_ else 100 * inter / uni
    }
    list(pct = pct, n = nn)
  }
  
  pep_pres <- (!is.na(m_union_pep))
  pro_pres <- (!is.na(m_union_prot))
  rownames(pep_pres) <- rownames(m_union_pep); colnames(pep_pres) <- colnames(m_union_pep)
  rownames(pro_pres) <- rownames(m_union_prot); colnames(pro_pres) <- colnames(m_union_prot)
  
  pep_sets <- sets_from_presence(pep_pres)
  pro_sets <- sets_from_presence(pro_pres)
  
  pep <- jaccard_from_sets(pep_sets$sets)
  pro <- jaccard_from_sets(pro_sets$sets)
  
  df <- CJ(Row = groups, Col = groups)
  df[, `:=`(Value = NA_real_, Label = "")]
  
  nG <- length(groups)
  for (i in seq_len(nG)) for (j in seq_len(nG)) {
    r <- groups[i]; c <- groups[j]
    if (i == j) {
      df[Row == r & Col == c, `:=`(
        Value = NA_real_,
        Label = paste0(fmt_int(pro_sets$diag[r]), " proteins\n", fmt_int(pep_sets$diag[r]), " peptides")
      )]
    } else if (i < j) {
      pct <- pro$pct[r, c]; nn <- pro$n[r, c]
      df[Row == r & Col == c, `:=`(Value = pct, Label = paste0(fmt_int(nn), " proteins"))]
    } else {
      pct <- pep$pct[r, c]; nn <- pep$n[r, c]
      df[Row == r & Col == c, `:=`(Value = pct, Label = paste0(fmt_int(nn), " peptides"))]
    }
  }
  
  df[, Row := factor(Row, levels = rev(groups))]
  df[, Col := factor(Col, levels = groups)]
  
  p <- make_overlap_heatmap_plot(df, fill_limits = c(0,100), text_size = 3.4)
  
  if (isTRUE(return_plot)) return(p)
  if (!is.null(out_png)) save_pub(p, out_png, width = 10.5, height = 6.0)
  invisible(NULL)
  
}

# ==============================
# Scatter / regression
# ==============================
orthogonal_fit <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3) return(list(m = NA_real_, b = NA_real_))
  
  xbar <- mean(x); ybar <- mean(y)
  xc <- x - xbar; yc <- y - ybar
  
  Sxx <- mean(xc^2)
  Syy <- mean(yc^2)
  Sxy <- mean(xc*yc)
  
  if (!is.finite(Sxy) || abs(Sxy) < .Machine$double.eps) {
    return(list(m = NA_real_, b = NA_real_))
  }
  
  disc <- (Syy - Sxx)^2 + 4 * Sxy^2
  m <- ((Syy - Sxx) + sqrt(disc)) / (2 * Sxy)
  b <- ybar - m * xbar
  list(m = m, b = b)
}

fit_y_eq_mx_plus_b <- function(x, y, method = c("lm", "orthogonal")) {
  method <- match.arg(method)
  if (method == "lm") {
    fit <- lm(y ~ x)
    m <- unname(coef(fit)[["x"]])
    b <- unname(coef(fit)[["(Intercept)"]])
    r2 <- summary(fit)$r.squared
    return(list(m = m, b = b, r2 = r2))
  } else {
    fo <- orthogonal_fit(x, y)
    r <- suppressWarnings(cor(x, y, method = "pearson", use = "complete.obs"))
    r2 <- r^2
    return(list(m = fo$m, b = fo$b, r2 = r2))
  }
}

get_matched_sample_pairs <- function(sample_meta, g1, g2) {
  a <- sample_meta[Group == g1 & !is.na(SampleKey) & nzchar(SampleKey), .(Sample1 = Sample, SampleKey)]
  b <- sample_meta[Group == g2 & !is.na(SampleKey) & nzchar(SampleKey), .(Sample2 = Sample, SampleKey)]
  merge(a, b, by = "SampleKey", allow.cartesian = TRUE)
}

maybe_sample_ids <- function(m, max_ids) {
  if (!is.finite(max_ids) || max_ids <= 0) return(m)
  if (nrow(m) <= max_ids) return(m)
  set.seed(1)
  keep <- sample.int(nrow(m), max_ids)
  m[keep, , drop = FALSE]
}

collect_scatter_points_all_ids <- function(m_any, sample_meta, g1, g2,
                                           use_log2 = TRUE,
                                           max_ids = 500000) {
  pairs <- get_matched_sample_pairs(sample_meta, g1, g2)
  if (nrow(pairs) == 0) return(NULL)
  
  # speed: sample IDs if enormous
  m_use <- maybe_sample_ids(m_any, max_ids = max_ids)
  
  keys <- unique(pairs$SampleKey)
  g1_samps <- pairs[, .(S1 = unique(Sample1)), by = SampleKey]
  g2_samps <- pairs[, .(S2 = unique(Sample2)), by = SampleKey]
  
  out_list <- vector("list", length(keys))
  
  for (i in seq_along(keys)) {
    k <- keys[i]
    s1 <- intersect(g1_samps[SampleKey == k, S1][[1]], colnames(m_use))
    s2 <- intersect(g2_samps[SampleKey == k, S2][[1]], colnames(m_use))
    if (length(s1) == 0 || length(s2) == 0) next
    
    x <- if (length(s1) == 1) m_use[, s1] else rowMeans(m_use[, s1, drop = FALSE], na.rm = TRUE)
    y <- if (length(s2) == 1) m_use[, s2] else rowMeans(m_use[, s2, drop = FALSE], na.rm = TRUE)
    
    ok <- is.finite(x) & is.finite(y) & x > 0 & y > 0
    if (sum(ok) < 10) next
    x <- x[ok]; y <- y[ok]
    
    if (use_log2) {
      x <- log2(x); y <- log2(y)
      ok2 <- is.finite(x) & is.finite(y)
      if (sum(ok2) < 10) next
      x <- x[ok2]; y <- y[ok2]
    }
    
    out_list[[i]] <- data.table(x = x, y = y, SampleKey = k)
  }
  
  pts <- rbindlist(out_list, use.names = TRUE, fill = TRUE)
  if (nrow(pts) == 0) return(NULL)
  pts
}

scatter_plot_with_stats <- function(pts, out_png,
                                    xlab, ylab,
                                    method = c("lm", "orthogonal"),
                                    max_points_plot = 250000,
                                    bins = 350) {
  method <- match.arg(method)
  if (is.null(pts) || nrow(pts) < 50) return(invisible(NULL))
  
  # Fit on all collected points (already maybe ID-sampled upstream for speed)
  fit <- fit_y_eq_mx_plus_b(pts$x, pts$y, method = method)
  m <- fit$m; b <- fit$b; r2 <- fit$r2
  n_all <- nrow(pts)
  
  # Plot downsample (visual only)
  pts_plot <- pts
  if (is.finite(max_points_plot) && nrow(pts_plot) > max_points_plot) {
    set.seed(1)
    pts_plot <- pts_plot[sample.int(nrow(pts_plot), max_points_plot)]
  }
  
  lab <- sprintf("y = %.3f x %s %.3f\nR² = %.3f\n%s points",
                 m, ifelse(b >= 0, "+", "-"), abs(b),
                 r2, format(n_all, big.mark = ","))
  
  # label anchor in upper-left
  lx <- quantile(pts_plot$x, 0.02, na.rm = TRUE)
  ly <- quantile(pts_plot$y, 0.98, na.rm = TRUE)
  
  p <- ggplot(pts_plot, aes(x = x, y = y)) +
    stat_bin2d(bins = bins) +
    scale_fill_gradient(
      low = scatter_fill_low,
      high = scatter_fill_high,
      trans = scatter_density_trans
    ) +
    guides(fill = "none") +
    {
      if (method == "lm") geom_smooth(method = "lm", se = FALSE, linewidth = 1.05)
      else geom_abline(slope = m, intercept = b, linewidth = 1.05)
    } +
    labs(title = NULL, x = xlab, y = ylab) +
    geom_label(
      data = data.frame(x = lx, y = ly, lab = lab),
      aes(x = x, y = y, label = lab),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1,
      size = 5.0,
      label.size = 0.35,
      label.padding = unit(0.30, "lines"),
      color = "black",
      fill = "white",
      fontface = "bold"
    ) +
    pub_theme(base_size = 13)
  
  # large export size = less pixelated
  save_pub(p, out_png, width = 9.2, height = 7.3)
  invisible(list(m = m, b = b, r2 = r2, n = n_all, method = method))
}

compute_r2_per_samplepair <- function(m_any, sample_meta, g1, g2,
                                      min_points = 50,
                                      use_log2 = TRUE,
                                      max_ids = 500000) {
  pairs <- get_matched_sample_pairs(sample_meta, g1, g2)
  if (nrow(pairs) == 0) return(data.table())
  
  m_use <- maybe_sample_ids(m_any, max_ids = max_ids)
  
  keys <- unique(pairs$SampleKey)
  g1_samps <- pairs[, .(S1 = unique(Sample1)), by = SampleKey]
  g2_samps <- pairs[, .(S2 = unique(Sample2)), by = SampleKey]
  
  out <- vector("list", length(keys))
  for (i in seq_along(keys)) {
    k <- keys[i]
    s1 <- intersect(g1_samps[SampleKey == k, S1][[1]], colnames(m_use))
    s2 <- intersect(g2_samps[SampleKey == k, S2][[1]], colnames(m_use))
    if (length(s1) == 0 || length(s2) == 0) next
    
    x <- if (length(s1) == 1) m_use[, s1] else rowMeans(m_use[, s1, drop = FALSE], na.rm = TRUE)
    y <- if (length(s2) == 1) m_use[, s2] else rowMeans(m_use[, s2, drop = FALSE], na.rm = TRUE)
    
    ok <- is.finite(x) & is.finite(y) & x > 0 & y > 0
    if (sum(ok) < min_points) next
    x <- x[ok]; y <- y[ok]
    
    if (use_log2) {
      x <- log2(x); y <- log2(y)
      ok2 <- is.finite(x) & is.finite(y)
      if (sum(ok2) < min_points) next
      x <- x[ok2]; y <- y[ok2]
    }
    
    r <- suppressWarnings(cor(x, y, method = "pearson"))
    out[[i]] <- data.table(SampleKey = k, N = length(x), R2 = r^2)
  }
  rbindlist(out, fill = TRUE)
}

plot_r2_violin_group_pairs <- function(r2_dt, out_png) {
  if (nrow(r2_dt) == 0) return(invisible(NULL))
  
  ord <- r2_dt[, .(Med = median(R2, na.rm = TRUE)), by = Pair][order(-Med)]$Pair
  r2_dt[, Pair := factor(Pair, levels = ord)]
  
  p <- ggplot(r2_dt, aes(x = Pair, y = R2)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.75) +
    geom_point(position = position_jitter(width = 0.12), size = 1.0, alpha = 0.35) +
    stat_summary(fun = median, geom = "point", shape = 23, fill = "white", color = "black", size = 2.6) +
    stat_summary(fun = mean, geom = "point", shape = 4, color = "black", stroke = 1.0, size = 2.8) +
    labs(title = NULL, x = NULL, y = expression(R^2)) +
    pub_theme(base_size = 12, strip_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  save_pub(p, out_png, width = 13.5, height = 6.5)
}

# ==============================
# IDs by group violin (only)
# ==============================
plot_group_all_ids <- function(id_df, group_map, out_png) {
  id_df <- as.data.table(id_df)
  group_map <- as.data.table(group_map)
  id_df <- merge(id_df, group_map[, .(Sample, Group)], by = "Sample", all.x = TRUE, allow.cartesian = TRUE)
  
  ord <- id_df[Search == ALL_SEARCHES_LABEL,
               .(Med = median(IDs, na.rm = TRUE)), by = Group][order(-Med)]$Group
  if (length(ord) == 0) ord <- sort(unique(id_df$Group))
  id_df[, Group := factor(Group, levels = ord)]
  
  search_levels <- c(setdiff(unique(id_df$Search), ALL_SEARCHES_LABEL), ALL_SEARCHES_LABEL)
  id_df[, Search := factor(Search, levels = search_levels)]
  
  p <- ggplot(id_df, aes(x = Group, y = IDs, fill = Search)) +
    geom_violin(position = position_dodge(width = 0.85), trim = TRUE, scale = "width", alpha = 0.75) +
    geom_point(aes(color = Search),
               position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.85),
               size = 1.1, alpha = 0.55) +
    stat_summary(fun = median, geom = "point", position = position_dodge(width = 0.85), shape = 23, fill = "white", color = "black", size = 2.5) +
    stat_summary(fun = mean, geom = "point", position = position_dodge(width = 0.85), shape = 4, color = "black", stroke = 1.0, size = 2.7) +
    labs(title = NULL, x = "Group", y = "# IDs", fill = NULL, color = NULL) +
    pub_theme(base_size = 12, strip_size = 16)
  
  save_pub(p, out_png, width = 12.5, height = 6.5)
  
  id_df[, .(
    N = sum(is.finite(IDs)),
    Median = median(IDs, na.rm = TRUE),
    Mean = mean(IDs, na.rm = TRUE),
    Q25 = quantile(IDs, 0.25, na.rm = TRUE),
    Q75 = quantile(IDs, 0.75, na.rm = TRUE)
  ), by = .(Search, Group)]
}

# ==============================
# Completeness (combined Protein+Peptide)
# ==============================
bin_completeness <- function(pct) {
  if (!is.finite(pct)) return(NA_character_)
  if (pct >= 100) return("100%")
  if (pct >= 95)  return("95-100%")
  if (pct >= 80)  return("80-95%")
  if (pct >= 50)  return("50-80%")
  if (pct >= 20)  return("20-50%")
  if (pct > 0)    return("0-20%")
  NA_character_  # 0% excluded
}

compute_completeness_dt <- function(mats_named, group_map, level_name) {
  gm <- as.data.table(group_map)
  searches <- names(mats_named)
  
  out <- vector("list", length(searches))
  for (i in seq_along(searches)) {
    s <- searches[i]
    m <- mats_named[[s]]
    pres <- !is.na(m)
    
    tmp <- vector("list", length(unique(gm$Group)))
    gnames <- sort(unique(gm$Group))
    for (k in seq_along(gnames)) {
      g <- gnames[k]
      smps <- gm[Group == g, Sample]
      smps <- intersect(smps, colnames(pres))
      if (length(smps) == 0) next
      pct <- 100 * rowMeans(pres[, smps, drop = FALSE], na.rm = TRUE)
      bins <- vapply(pct, bin_completeness, character(1))
      keep <- !is.na(bins)
      if (!any(keep)) next
      dt <- data.table(Bin = bins[keep])
      tmp[[k]] <- dt[, .(Count = .N), by = Bin][, `:=`(Search = s, Level = level_name, Group = g)]
    }
    out[[i]] <- rbindlist(tmp, fill = TRUE)
  }
  dt_search <- rbindlist(out, fill = TRUE)
  
  union_pres <- build_union_presence(mats_named)
  tmpU <- vector("list", length(unique(gm$Group)))
  gnames <- sort(unique(gm$Group))
  for (k in seq_along(gnames)) {
    g <- gnames[k]
    smps <- gm[Group == g, Sample]
    smps <- intersect(smps, colnames(union_pres))
    if (length(smps) == 0) next
    pct <- 100 * rowMeans(union_pres[, smps, drop = FALSE], na.rm = TRUE)
    bins <- vapply(pct, bin_completeness, character(1))
    keep <- !is.na(bins)
    if (!any(keep)) next
    dt <- data.table(Bin = bins[keep])
    tmpU[[k]] <- dt[, .(Count = .N), by = Bin][, `:=`(Search = ALL_SEARCHES_LABEL, Level = level_name, Group = g)]
  }
  dt_union <- rbindlist(tmpU, fill = TRUE)
  
  res <- rbind(dt_search, dt_union, fill = TRUE)
  if (nrow(res) == 0) return(res)
  
  # Completeness desired bottom -> top
  comp_levels <- c("100%", "95-100%", "80-95%", "50-80%", "20-50%", "0-20%")
  
  res[, Bin_chr := as.character(Bin)]
  res[, StackRank := match(Bin_chr, comp_levels)]   # 1=bottom (100%), 6=top (0-20%)
  res[, Bin := factor(Bin_chr, levels = comp_levels, ordered = TRUE)]
}

plot_completeness_combined <- function(dt_all, out_png, return_plot = FALSE) {
  dt_all <- as.data.table(dt_all)
  if (nrow(dt_all) == 0) return(invisible(NULL))
  
  # Desired bottom -> top
  comp_levels <- c("100%", "95-100%", "80-95%", "50-80%", "20-50%", "0-20%")
  
  # Search order: strategies then Union at end
  search_levels <- c(setdiff(unique(dt_all$Search), ALL_SEARCHES_LABEL), ALL_SEARCHES_LABEL)
  dt_all[, Search := factor(Search, levels = search_levels)]
  
  # Group order alphabetically
  ord <- sort(unique(as.character(dt_all$Group)))
  dt_all[, Group := factor(Group, levels = ord)]
  
  # Stack rank
  dt_all[, Bin_chr := as.character(Bin)]
  dt_all[, StackRank := match(Bin_chr, comp_levels)]
  dt_all <- dt_all[!is.na(StackRank)]
  dt_all[, Bin := factor(Bin_chr, levels = comp_levels, ordered = TRUE)]
  
  # Count labels (hide unreadable) — thresholds adapt by Level
  dt_all[, Total := sum(Count, na.rm = TRUE), by = .(Level, Search, Group)]
  
  # Peptide panels are ~10x bigger, so use stricter absolute threshold there
  dt_all[, min_abs := ifelse(Level == "Peptide", 100L, 70L)]
  dt_all[, min_frac := ifelse(Level == "Peptide", 0.06, 0.06)]  # peptide can use smaller fraction, protein a bit bigger
  
  dt_all[, LabelCnt := ifelse(
    Count >= min_abs & Count >= (min_frac * Total),
    fmt_int(Count),
    ""
  )]
  
  
  # ---- MANUAL STACK
  setorder(dt_all, Level, Search, Group, StackRank)
  dt_all[, `:=`(
    ymax = cumsum(Count),
    ymin = c(0, head(cumsum(Count), -1))
  ), by = .(Level, Search, Group)]
  
  dt_all[, ylab := (ymin + ymax) / 2]
  
  dt_all[, x := as.numeric(Group)]
  bar_w <- 0.85
  dt_all[, `:=`(xmin = x - bar_w/2, xmax = x + bar_w/2)]
  
  pal <- comp_palette()
  
  p <- ggplot(dt_all) +
    geom_rect(
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Bin),
      color = "black", linewidth = 0.25
    ) +
    facet_grid(Level ~ Search, scales = "free_y", switch = "y") +
    scale_x_continuous(breaks = seq_along(ord), labels = ord) +
    labs(title = NULL, x = "Group", y = NULL, fill = NULL) +
    scale_fill_manual(values = pal, breaks = rev(comp_levels), drop = FALSE) +
    pub_theme(base_size = 12, strip_size = 16) +
    theme(
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.title.y = element_blank()
    )
  
  if (isTRUE(return_plot)) return(p)
  save_pub(p, out_png, width = 15.5, height = max(7.5, 2.6 * length(unique(dt_all$Level))))
}


# ==============================
# Peptides-per-protein stacked bar (per group) + UNION
# ==============================
infer_pr_protein_peptide_cols <- function(dt_names) {
  prot_candidates <- c("Protein.Group", "ProteinGroup", "Protein.Groups", "PG.ProteinGroups", "Protein")
  pep_candidates  <- c("Stripped.Sequence", "StrippedSequence", "ModifiedPeptideSequence", "Modified.Sequence", "Modified sequence", "Sequence")
  prot_col <- prot_candidates[prot_candidates %in% dt_names]
  pep_col  <- pep_candidates[pep_candidates %in% dt_names]
  list(
    prot_col = if (length(prot_col)) prot_col[1] else NA_character_,
    pep_col  = if (length(pep_col))  pep_col[1]  else NA_character_
  )
}

bin_peptides_per_protein <- function(npep) {
  if (!is.finite(npep)) return(NA_character_)
  if (npep <= 0) return(NA_character_)
  if (npep == 1) return("1")
  if (npep == 2) return("2")
  if (npep == 3) return("3")
  if (npep == 4) return("4")
  if (npep == 5) return("5")
  if (npep <= 10) return("6-10")
  ">10"
}

peptides_per_protein_stacked_by_group <- function(pr_dt, group_map, search_name) {
  cols <- infer_pr_protein_peptide_cols(names(pr_dt))
  if (is.na(cols$prot_col) || is.na(cols$pep_col)) return(list(bins = data.table(), median = data.table(), dist = data.table()))
  
  prot_col <- cols$prot_col
  pep_col  <- cols$pep_col
  
  exclude <- c(prot_col, pep_col)
  sample_cols <- numeric_cols_only(pr_dt, exclude = exclude,
                                   diag_tag = paste0("PepsPerProt_", search_name),
                                   diag_dir = NULL)
  sample_cols <- intersect(sample_cols, group_map$Sample)
  if (length(sample_cols) == 0) return(list(bins = data.table(), median = data.table(), dist = data.table()))
  
  pr_dt[[prot_col]] <- sanitize_ids(pr_dt[[prot_col]])
  pr_dt[[pep_col]]  <- sanitize_ids(pr_dt[[pep_col]])
  pr_dt <- pr_dt[!is.na(pr_dt[[prot_col]]) & nzchar(pr_dt[[prot_col]]) &
                   !is.na(pr_dt[[pep_col]])  & nzchar(pr_dt[[pep_col]])]
  
  for (cname in sample_cols) pr_dt[[cname]] <- suppressWarnings(as.numeric(pr_dt[[cname]]))
  
  gm <- as.data.table(group_map)
  groups <- sort(unique(gm$Group))
  
  out_bins <- vector("list", length(groups))
  out_med  <- vector("list", length(groups))
  out_dist <- vector("list", length(groups))
  
  for (i in seq_along(groups)) {
    g <- groups[i]
    smps <- gm[Group == g, Sample]
    smps <- intersect(smps, sample_cols)
    if (length(smps) == 0) next
    
    mm <- as.matrix(pr_dt[, ..smps])
    mode(mm) <- "numeric"
    present <- rowSums(is.finite(mm) & mm > 0) > 0
    if (!any(present)) next
    
    dt <- pr_dt[present, .(Protein = get(prot_col), Peptide = get(pep_col))]
    dt <- unique(dt)
    
    cnt <- dt[, .(Npep = uniqueN(Peptide)), by = Protein]
    cnt <- cnt[is.finite(Npep) & Npep > 0]
    
    if (nrow(cnt) == 0) next
    
    # bins
    cnt[, Bin := vapply(Npep, bin_peptides_per_protein, character(1))]
    cnt <- cnt[!is.na(Bin)]
    bins_dt <- cnt[, .(Count = .N), by = Bin][, `:=`(Search = search_name, Level = "Protein", Group = g)]
    bins_dt[, Bin := factor(as.character(Bin), levels = ppp_levels, ordered = TRUE)]
    out_bins[[i]] <- bins_dt
    
    # median overlay
    out_med[[i]] <- data.table(
      Search = search_name,
      Group = g,
      MedianNpep = median(cnt$Npep, na.rm = TRUE),
      NProteins = nrow(cnt)
    )
    
    # distribution for boxplot
    out_dist[[i]] <- cnt[, .(Search = search_name, Group = g, Protein, Npep)]
  }
  
  list(
    bins   = rbindlist(out_bins, fill = TRUE),
    median = rbindlist(out_med,  fill = TRUE),
    dist   = rbindlist(out_dist, fill = TRUE)
  )
}

build_ppp_union_from_searches <- function(pr_raw_list, group_map) {
  gm <- as.data.table(group_map)
  groups <- sort(unique(gm$Group))
  
  out_bins <- vector("list", length(groups))
  out_med  <- vector("list", length(groups))
  out_dist <- vector("list", length(groups))
  
  for (i in seq_along(groups)) {
    g <- groups[i]
    smps <- gm[Group == g, Sample]
    
    pairs_all <- list()
    
    for (s in names(pr_raw_list)) {
      pr_dt <- pr_raw_list[[s]]
      cols <- infer_pr_protein_peptide_cols(names(pr_dt))
      if (is.na(cols$prot_col) || is.na(cols$pep_col)) next
      
      prot_col <- cols$prot_col
      pep_col  <- cols$pep_col
      
      sample_cols <- intersect(smps, names(pr_dt))
      sample_cols <- sample_cols[!vapply(sample_cols, is_non_sample_col, logical(1))]
      if (length(sample_cols) == 0) next
      
      for (cname in sample_cols) pr_dt[[cname]] <- suppressWarnings(as.numeric(pr_dt[[cname]]))
      
      mm <- as.matrix(pr_dt[, ..sample_cols])
      mode(mm) <- "numeric"
      present <- rowSums(is.finite(mm) & mm > 0) > 0
      if (!any(present)) next
      
      dt <- pr_dt[present, .(Protein = sanitize_ids(get(prot_col)),
                             Peptide = sanitize_ids(get(pep_col)))]
      dt <- dt[!is.na(Protein) & nzchar(Protein) & !is.na(Peptide) & nzchar(Peptide)]
      if (nrow(dt) == 0) next
      
      pairs_all[[length(pairs_all) + 1]] <- unique(dt)
    }
    
    if (length(pairs_all) == 0) next
    pairs <- unique(rbindlist(pairs_all))
    
    cnt <- pairs[, .(Npep = uniqueN(Peptide)), by = Protein]
    cnt <- cnt[is.finite(Npep) & Npep > 0]
    if (nrow(cnt) == 0) next
    
    cnt[, Bin := vapply(Npep, bin_peptides_per_protein, character(1))]
    cnt <- cnt[!is.na(Bin)]
    
    bins_dt <- cnt[, .(Count = .N), by = Bin][, `:=`(Search = ALL_SEARCHES_LABEL, Level = "Protein", Group = g)]
    bins_dt[, Bin := factor(as.character(Bin), levels = ppp_levels, ordered = TRUE)]
    out_bins[[i]] <- bins_dt
    
    out_med[[i]] <- data.table(
      Search = ALL_SEARCHES_LABEL,
      Group = g,
      MedianNpep = median(cnt$Npep, na.rm = TRUE),
      NProteins = nrow(cnt)
    )
    
    out_dist[[i]] <- cnt[, .(Search = ALL_SEARCHES_LABEL, Group = g, Protein, Npep)]
  }
  
  list(
    bins   = rbindlist(out_bins, fill = TRUE),
    median = rbindlist(out_med,  fill = TRUE),
    dist   = rbindlist(out_dist, fill = TRUE)
  )
}

plot_peptides_per_protein_combined <- function(ppp_dt2, out_png,
                                               return_plot = FALSE,
                                               show_search_strip = TRUE) {
  ppp_dt2 <- as.data.table(ppp_dt2)
  if (nrow(ppp_dt2) == 0) return(invisible(NULL))
  
  # Desired bottom -> top
  ppp_levels <- c(">10", "6-10", "5", "4", "3", "2", "1")
  
  # Search order: strategies then Union at end
  search_levels <- c(setdiff(unique(ppp_dt2$Search), ALL_SEARCHES_LABEL), ALL_SEARCHES_LABEL)
  ppp_dt2[, Search := factor(Search, levels = search_levels)]
  
  # Group order by UNION total proteins
  ord <- sort(unique(as.character(ppp_dt2$Group)))
  ppp_dt2[, Group := factor(Group, levels = ord)]
  
  # Rank bins in the desired stack order
  ppp_dt2[, Bin_chr := as.character(Bin)]
  ppp_dt2[, StackRank := match(Bin_chr, ppp_levels)]
  ppp_dt2 <- ppp_dt2[!is.na(StackRank)]
  
  # Count labels (hide unreadable)
  ppp_dt2[, Total := sum(Count, na.rm = TRUE), by = .(Search, Group)]
  
  min_abs  <- 100      # at least 100 proteins in the segment
  min_frac <- 0.05     # OR at least 5% of the bar total
  ppp_dt2[, LabelCnt := ifelse(
    Count >= min_abs & Count >= (min_frac * Total),
    fmt_int(Count),
    ""
  )]
  
  
  # ---- MANUAL STACK: compute ymin/ymax in StackRank order
  setorder(ppp_dt2, Search, Group, StackRank)
  ppp_dt2[, `:=`(
    ymax = cumsum(Count),
    ymin = c(0, head(cumsum(Count), -1))
  ), by = .(Search, Group)]

  # label y position
  ppp_dt2[, ylab := (ymin + ymax) / 2]
  
  # x positions for geom_rect
  ppp_dt2[, x := as.numeric(Group)]
  bar_w <- 0.85
  ppp_dt2[, `:=`(xmin = x - bar_w/2, xmax = x + bar_w/2)]
  
  # palette mapped by bin name
  pal <- ppp_palette()
  # ensure fill levels exist even if missing in some facets
  ppp_dt2[, Bin_plot := factor(Bin_chr, levels = ppp_levels, ordered = TRUE)]
  
  p <- ggplot(ppp_dt2) +
    geom_rect(
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Bin_plot),
      color = "black", linewidth = 0.25
    ) +
    facet_wrap(~ Search, nrow = 1) +   # shared y-axis
    scale_x_continuous(breaks = seq_along(ord), labels = ord) +
    labs(title = NULL, x = "Group", y = "Proteins (count)", fill = NULL) +
    scale_fill_manual(values = pal, breaks = rev(ppp_levels), drop = FALSE) +
    pub_theme(base_size = 12, strip_size = 16) +
    theme(strip.background = element_rect(fill = "white", color = "white"))

  if (!isTRUE(show_search_strip)) {
    p <- p + theme(
      strip.text = element_blank(),
      strip.background = element_blank(),
      axis.title.y = element_blank()
    )
  }
  
  
  if (isTRUE(return_plot)) return(p)
  save_pub(p, out_png, width = 15.5, height = 7.0)
  
}

plot_completeness_and_ppp_combined <- function(comp_dt, ppp_dt, out_png) {
  
  # need these packages for perfect axis + legend alignment
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("Package 'cowplot' is required for axis/legend alignment. Install with install.packages('cowplot').")
  }
  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("Package 'grid' is required.")
  }
  
  # Build plots (no saving yet)
  p_comp <- plot_completeness_combined(comp_dt, out_png = NULL, return_plot = TRUE) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),     # remove "IDs" label as requested
      legend.position = "right"
    )
  
  p_ppp <- plot_peptides_per_protein_combined(
    ppp_dt, out_png = NULL, return_plot = TRUE,
    show_search_strip = FALSE              # removes the strip on top of PPP if you added that option
  ) +
    theme(
      legend.position = "right"
    )
  
  # ---- Extract legends
  leg_comp <- cowplot::get_legend(p_comp)
  leg_ppp  <- cowplot::get_legend(p_ppp)
  
  # ---- Remove legends from plots (so plots align perfectly)
  p_comp_noleg <- p_comp + theme(legend.position = "none")
  p_ppp_noleg  <- p_ppp  + theme(legend.position = "none")
  
  # ---- FORCE identical plot margins (helps keep left columns consistent)
  common_margins <- theme(plot.margin = margin(6, 6, 6, 10))
  p_comp_noleg <- p_comp_noleg + common_margins
  p_ppp_noleg  <- p_ppp_noleg  + common_margins
  
  # ---- Align the left axes/panels (THIS fixes the y-axis offset)
  aligned <- cowplot::align_plots(
    p_comp_noleg, p_ppp_noleg,
    align = "v", axis = "l"
  )
  p_comp_al <- aligned[[1]]
  p_ppp_al  <- aligned[[2]]
  
  # ---- Stack plots: completeness height 2x PPP so PPP ~= "Protein row" height
  left_stack <- cowplot::plot_grid(
    p_comp_al, p_ppp_al,
    ncol = 1,
    rel_heights = c(2, 1),
    align = "v",
    axis = "l"
  )
  
  # ---- Stack legends into one right column (perfectly aligned)
  right_legs <- cowplot::plot_grid(
    leg_comp, leg_ppp,
    ncol = 1,
    rel_heights = c(2, 1)
  )
  
  # ---- Combine left plots + right legends
  final <- cowplot::plot_grid(
    left_stack, right_legs,
    ncol = 2,
    rel_widths = c(1, 0.16)
  )
  
  # Save
  ggsave(out_png, final, width = 15.5, height = 12.0, dpi = PUB_DPI, units = "in")
  ggsave(sub("\\.png$", ".pdf", out_png, ignore.case = TRUE),
         final, width = 15.5, height = 12.0, device = cairo_pdf, units = "in")
  
  invisible(final)
}

# ==============================
# Pair selection (prompt ONCE)
# ==============================
prompt_group_pair_selection <- function(groups) {
  groups <- sort(unique(groups))
  if (length(groups) < 2) return(data.table(g1 = character(), g2 = character()))
  
  pairs <- CJ(i = seq_along(groups), j = seq_along(groups))[i < j]
  pairs[, g1 := groups[i]]
  pairs[, g2 := groups[j]]
  pairs[, idx := .I]
  
  cat("\nSelect group comparisons (used for BOTH scatters and R² violin)\n")
  cat("Type:\n")
  cat("  all                -> include all pairs\n")
  cat("  1,2,5-7            -> choose specific pairs\n")
  cat("  (blank)            -> auto-pick top pairs by matched keys\n\n")
  
  for (k in seq_len(nrow(pairs))) {
    cat(sprintf("  [%d] %s vs %s\n", pairs$idx[k], pairs$g1[k], pairs$g2[k]))
  }
  
  ans <- trimws(readline("\nEnter selection: "))
  if (!nzchar(ans)) {
    return(data.table(g1 = character(), g2 = character()))  # signals "auto"
  }
  if (tolower(ans) == "all") {
    return(pairs[, .(g1, g2)])
  }
  
  toks <- unlist(strsplit(gsub("\\s+", "", ans), ",", fixed = TRUE))
  idxs <- integer()
  for (t in toks) {
    if (!nzchar(t)) next
    if (grepl("^[0-9]+-[0-9]+$", t)) {
      ab <- as.integer(strsplit(t, "-", fixed = TRUE)[[1]])
      if (length(ab) == 2 && all(is.finite(ab))) idxs <- c(idxs, seq(min(ab), max(ab)))
    } else if (grepl("^[0-9]+$", t)) {
      idxs <- c(idxs, as.integer(t))
    }
  }
  idxs <- sort(unique(idxs))
  idxs <- idxs[idxs >= 1 & idxs <= nrow(pairs)]
  
  if (length(idxs) == 0) {
    cat("No valid indices detected; using auto selection.\n")
    return(data.table(g1 = character(), g2 = character()))
  }
  
  pairs[idx %in% idxs, .(g1, g2)]
}

# ==============================
# Hardcoded group-comparison plan (2-protease paper) -- PI-specified 2026-06
# ------------------------------------------------------------------------------
# Replaces the interactive pair picker so every run is deterministic. Detected
# groups are produced by auto_guess_group() as:
#     Mix_20min, Mix_30min   (Combined / mixed-protease digest; 20- vs 30-min)
#     Kra_<date>             (Krakatoa single-protease; per acquisition batch)
#     Try_<date>             (Trypsin  single-protease; per acquisition batch)
#
# Comparison rules (from the PI):
#   * Trypsin is always compared to Krakatoa AND to Combined (Mix).
#   * Krakatoa is always compared to Trypsin AND to Combined (Mix).
#   * "The timed runs are the same" -> condition/batch-matched:
#       - Try vs Kra are matched on the shared acquisition-batch token (date)
#         when one exists (e.g. Try_250807 vs Kra_250807); else all-vs-all.
#       - Combined (Mix) carries gradient tokens (20min/30min) that do NOT match
#         the single-protease batch tokens, so every protease batch is compared
#         to every Combined gradient.
#       - The two Combined gradients are compared to each other (20min vs 30min).
#   * "...etc" -> each protease is also compared across its own batches
#     (within-protease reproducibility).
# Robust to whichever exact dates/gradients are present; returns data.table(g1,g2).
# ==============================
classify_group_2p <- function(g) {
  prot <- if (grepl("^Try", g, ignore.case = TRUE)) "Try" else
          if (grepl("^Kra", g, ignore.case = TRUE)) "Kra" else
          if (grepl("^Mix", g, ignore.case = TRUE)) "Mix" else "Other"
  # Condition / set token = the trailing token after the last space or underscore.
  # Robust to BOTH the raw auto-labels (Try_250807 -> "250807", Mix_20min -> "20min")
  # AND the manuscript relabels (Trypsin 2 -> "2", Mix 30min -> "30min"). Cross-
  # protease batch-matching (Try vs Kra) keys on this token, so it must agree
  # between proteases in either scheme (251117/251117 or 2/2).
  cond <- sub(".*[ _]", "", g)
  list(prot = prot, cond = cond)
}

build_hardcoded_group_pairs <- function(groups) {
  groups <- sort(unique(groups))
  if (length(groups) < 2) return(data.table(g1 = character(), g2 = character()))
  info <- lapply(groups, classify_group_2p)
  G <- data.table(
    g    = groups,
    prot = vapply(info, `[[`, character(1), "prot"),
    cond = vapply(info, `[[`, character(1), "cond")
  )

  acc <- list()
  push <- function(a, b) {
    if (length(a) == 1L && length(b) == 1L && !is.na(a) && !is.na(b) &&
        nzchar(a) && nzchar(b) && a != b) {
      ab <- sort(c(a, b))
      acc[[length(acc) + 1L]] <<- data.table(g1 = ab[1], g2 = ab[2])
    }
  }

  protease_set <- c("Try", "Kra")

  # 1) Trypsin vs Krakatoa -- batch-matched on shared token, else all-vs-all.
  try_g <- G[prot == "Try"]; kra_g <- G[prot == "Kra"]
  shared <- intersect(try_g$cond, kra_g$cond)
  if (length(shared) > 0) {
    for (cc in shared) push(try_g[cond == cc, g][1], kra_g[cond == cc, g][1])
  } else {
    for (a in try_g$g) for (b in kra_g$g) push(a, b)
  }

  # 2) Each protease vs Combined (Mix) -- no shared token => all-vs-all.
  mix_g <- G[prot == "Mix"]
  for (a in G[prot %in% protease_set, g]) for (b in mix_g$g) push(a, b)

  # 3) Combined gradient comparison (the "20 min vs the other 20 min").
  if (nrow(mix_g) >= 2) {
    mc <- mix_g$g
    for (i in seq_along(mc)) for (j in seq_along(mc)) if (i < j) push(mc[i], mc[j])
  }

  # 4) Same protease across its own batches (reproducibility) -- the "etc".
  for (p in protease_set) {
    gg <- G[prot == p, g]
    if (length(gg) >= 2)
      for (i in seq_along(gg)) for (j in seq_along(gg)) if (i < j) push(gg[i], gg[j])
  }

  if (length(acc) == 0) return(data.table(g1 = character(), g2 = character()))
  unique(rbindlist(acc))
}

# ==============================
# Search selection for group scatter/R² (NEW)
# ==============================
prompt_search_selection <- function(search_names, include_union = TRUE, purpose = "scatter/R2") {
  opts <- search_names
  if (isTRUE(include_union)) opts <- c(opts, ALL_SEARCHES_LABEL)
  
  cat("\nSelect which search(es) to generate ", purpose, " for:\n", sep = "")
  cat("Type:\n")
  cat("  all                -> all searches (and Union if listed)\n")
  cat("  1,3,5-6            -> choose specific entries\n")
  cat("  (blank)            -> default: ", ALL_SEARCHES_LABEL, "\n\n", sep = "")
  
  for (i in seq_along(opts)) cat(sprintf("  [%d] %s\n", i, opts[i]))
  
  ans <- trimws(readline("\nEnter selection: "))
  if (!nzchar(ans)) return(ALL_SEARCHES_LABEL)
  if (tolower(ans) == "all") return(opts)
  
  toks <- unlist(strsplit(gsub("\\s+", "", ans), ",", fixed = TRUE))
  idxs <- integer()
  for (t in toks) {
    if (!nzchar(t)) next
    if (grepl("^[0-9]+-[0-9]+$", t)) {
      ab <- as.integer(strsplit(t, "-", fixed = TRUE)[[1]])
      if (length(ab) == 2 && all(is.finite(ab))) idxs <- c(idxs, seq(min(ab), max(ab)))
    } else if (grepl("^[0-9]+$", t)) {
      idxs <- c(idxs, as.integer(t))
    }
  }
  idxs <- sort(unique(idxs))
  idxs <- idxs[idxs >= 1 & idxs <= length(opts)]
  if (length(idxs) == 0) return(ALL_SEARCHES_LABEL)
  opts[idxs]
}

# ==============================
# HTML index
# ==============================
write_index_html <- function(path, sections) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  
  cat("<html><head><meta charset='utf-8'><title>DIA-NN Search Comparison</title></head><body>\n", file = con)
  cat("<h1>DIA-NN Search Comparison Report</h1>\n", file = con)
  cat("<p>Generated: ", format(Sys.time()), "</p>\n", file = con)
  cat("<p><b>Note:</b> All figures are exported as high-res PNG (600 dpi) and vector PDF.</p>\n", file = con)
  
  for (sec in names(sections)) {
    cat("<h2>", sec, "</h2>\n<ul>\n", file = con, sep = "")
    
    items <- sections[[sec]]
    if (is.null(items)) items <- list()
    
    # Normalize items into a list of 2-element lists
    norm_items <- list()
    
    if (is.list(items)) {
      # Case 1: items already looks like list(list("a","b"), ...)
      if (length(items) > 0 && is.list(items[[1]])) {
        norm_items <- items
      } else {
        # Case 2: items is list("a","b") or list(list("a","b")) weirdness
        if (length(items) >= 2 && all(vapply(items[1:2], function(x) is.character(x) && length(x) == 1, logical(1)))) {
          norm_items <- list(list(items[[1]], items[[2]]))
        } else {
          norm_items <- list()
        }
      }
    }
    
    for (it in norm_items) {
      if (!is.list(it) || length(it) < 2) next
      title <- it[[1]]
      href  <- it[[2]]
      if (!is.character(title) || !is.character(href) || length(title) != 1 || length(href) != 1) next
      
      cat("<li><a href='", href, "'>", title, "</a></li>\n", file = con, sep = "")
    }
    
    cat("</ul>\n", file = con)
  }
  
  cat("</body></html>\n", file = con)
}

# ==============================
# MAIN
# ==============================
tryCatch({
  
  message("Discovering search folders under: ", base_dir)
  search_paths <- discover_search_folders(base_dir, pg_file, pr_file)
  search_paths <- make_unique_names(search_paths)
  
  if (length(search_paths) < 2) stop("Need at least 2 valid search folders.")
  if (is.finite(max_searches) && length(search_paths) > max_searches) {
    message("Found ", length(search_paths), " searches; using first ", max_searches, " (alphabetical).")
    search_paths <- search_paths[seq_len(max_searches)]
  }
  
  search_paths <- prompt_search_short_names(search_paths)
  fwrite(data.table(FolderName = basename(search_paths), ShortName = names(search_paths), Path = as.character(search_paths)),
         file.path(out_dir, "search_name_mapping.csv"))
  
  message("\nLoading searches...")
  peptide_mats    <- list()  # collapsed by Modified.Sequence  (TRUE peptides)
  precursor_mats  <- list()  # uncollapsed pr_matrix rows      (precursors)
  protein_mats    <- list()
  pr_raw_list     <- list()
  
  sample_diag <- data.table()
  diag_dir_s <- file.path(out_dir, "Diagnostics")
  dir.create(diag_dir_s, showWarnings = FALSE)
  
  for (s in names(search_paths)) {
    sp <- search_paths[[s]]
    pr_path <- file.path(sp, pr_file)
    pg_path <- file.path(sp, pg_file)
    
    stop_if_missing(pr_path)
    stop_if_missing(pg_path)
    
    pr <- fread(pr_path, sep = "\t", header = TRUE, data.table = TRUE, showProgress = FALSE)
    pr_raw_list[[s]] <- pr
    
    # NEW: precursor-level (uncollapsed) from pr_matrix
    prec_m <- build_precursor_matrix_from_pr(pr, zero_as_na = zero_as_na,
                                             diag_tag = paste0("Precursor_", s), diag_dir = diag_dir_s)
    colnames(prec_m) <- make_unique_samples(colnames(prec_m), "Precursor", s, out_dir)
    precursor_mats[[s]] <- prec_m
    
    # NEW: peptide-level (collapsed) from pr_matrix by Modified.Sequence
    pep_m <- collapse_to_peptides(
      pr,
      collapse_method = peptide_collapse_method,
      zero_as_na = zero_as_na,
      peptide_id_candidates = c("Modified.Sequence"),
      diag_tag = paste0("Peptide_", s),
      diag_dir = diag_dir_s
    )
    # --- DIAGNOSTICS: peptide collapse check ---
    ms_col <- c("Modified.Sequence","ModifiedPeptideSequence","Modified.Sequence")[
      c("Modified.Sequence","ModifiedPeptideSequence","Modified.Sequence") %in% names(pr)
    ][1]
    
    if (!is.na(ms_col)) {
      n_rows_raw <- nrow(pr)
      n_unique_ms <- uniqueN(pr[[ms_col]])
      cat("\n[", s, "] pr rows=", n_rows_raw,
          " unique(", ms_col, ")=", n_unique_ms,
          " collapsed pep_m rows=", nrow(pep_m), "\n", sep = "")
    } else {
      cat("\n[", s, "] Could not find Modified.Sequence-like column in pr_matrix\n", sep = "")
    }
    colnames(pep_m) <- make_unique_samples(colnames(pep_m), "Peptide", s, out_dir)
    peptide_mats[[s]] <- pep_m
    
    
    pg <- read_diann_matrix(pg_path, id_col_guess = "Protein.Group",
                            diag_tag = paste0("Protein_", s), diag_dir = diag_dir_s)
    prot_m <- matrix_from_dt(pg$dt, pg$id_col, pg$num_cols, zero_as_na = zero_as_na)
    colnames(prot_m) <- make_unique_samples(colnames(prot_m), "Protein", s, out_dir)
    rn2 <- sanitize_ids(rownames(prot_m))
    prot_m <- prot_m[!is.na(rn2) & nzchar(rn2), , drop = FALSE]
    protein_mats[[s]] <- prot_m
    
    message("  ", s, " | peptides: ", nrow(pep_m), " | proteins: ", nrow(prot_m))
    sample_diag <- rbind(sample_diag, data.table(Search = s, Level = "Peptide", N_Samples = ncol(pep_m)))
    sample_diag <- rbind(sample_diag, data.table(Search = s, Level = "Protein", N_Samples = ncol(prot_m)))
  }
  
  fwrite(sample_diag, file.path(out_dir, "sample_counts_per_search.csv"))
  
  peptide_aligned <- align_samples_union(peptide_mats); peptide_mats <- peptide_aligned$mats
  protein_aligned <- align_samples_union(protein_mats); protein_mats <- protein_aligned$mats
  
  common_samples <- colnames(peptide_mats[[1]])
  common_samples <- filter_sample_like_names(common_samples)
  
  peptide_mats <- lapply(peptide_mats, function(m) m[, common_samples, drop = FALSE])
  protein_mats <- lapply(protein_mats, function(m) m[, common_samples, drop = FALSE])
  
  fwrite(data.table(Sample = common_samples), file.path(out_dir, "sample_list_common.csv"))
  message("\nFinal union sample columns: ", length(common_samples))
  
  # Groups
  group_map <- build_group_map_auto(common_samples)
  group_map <- assist_group_unknowns_learn_scheme(group_map, k = 3)
  group_map <- maybe_use_manual_group_table(group_map, out_dir)
  # Hardcoded manuscript labels (Trypsin 1/2, Krakatoa 1/2, Mix 20min/30min)
  # in place of the interactive prompt_group_renames().
  group_map <- apply_manuscript_group_labels(group_map)
  fwrite(group_map, file.path(out_dir, "sample_group_mapping.csv"))
  
  cat("\nGroup counts:\n")
  print(as.data.table(group_map)[, .N, by = Group][order(-N)])
  
  # ---- Search overlap combined triangle (proteins one side, peptides other)
  peptide_pres <- id_presence(peptide_mats)
  protein_pres <- id_presence(protein_mats)
  plot_search_overlap_combined(peptide_pres, protein_pres, file.path(out_dir, "SearchOverlap_Proteins_vs_Peptides_Triangles.png"))
  
  # ---- per level: IDs-by-group only
  analyze_level <- function(level_name, mats_named) {
    level_dir <- file.path(out_dir, level_name)
    dir.create(level_dir, showWarnings = FALSE)
    
    pres <- id_presence(mats_named)
    ov <- pairwise_overlap_pct(pres)
    fwrite(as.data.table(ov$pct, keep.rownames = "Search"), file.path(level_dir, "overlap_percent.csv"))
    fwrite(as.data.table(ov$n, keep.rownames = "Search"), file.path(level_dir, "overlap_n.csv"))
    
    id_df <- build_idcount_df_with_all_searches(mats_named, level_name)
    fwrite(id_df, file.path(level_dir, "id_counts_per_sample_with_all_searches.csv"))
    
    sum_ids <- plot_group_all_ids(
      id_df,
      group_map,
      file.path(level_dir, "ids_by_group_violin.png")
    )
    fwrite(sum_ids, file.path(level_dir, "ids_by_group_summary.csv"))
    
    list(level_dir = level_dir)
  }
  
  message("\nAnalyzing peptide level...")
  pep_res <- analyze_level("Peptide", peptide_mats)
  
  message("Analyzing protein level...")
  prot_res <- analyze_level("Protein", protein_mats)
  
  plot_group_all_ids_combined(
    id_df_pep  = fread(file.path(out_dir, "Peptide", "id_counts_per_sample_with_all_searches.csv")),
    id_df_prot = fread(file.path(out_dir, "Protein", "id_counts_per_sample_with_all_searches.csv")),
    group_map  = group_map,
    out_png    = file.path(out_dir, "IDs_by_group_Violin_Combined_Peptide_Protein.png")
  )
  
  # ---- Completeness (combined)
  comp_pep <- compute_completeness_dt(peptide_mats, group_map, "Peptide")
  comp_pro <- compute_completeness_dt(protein_mats, group_map, "Protein")
  comp_all <- rbind(comp_pep, comp_pro, fill = TRUE)
  if (nrow(comp_all) > 0) {
    fwrite(comp_all, file.path(out_dir, "Completeness_by_group_bins.csv"))
    plot_completeness_combined(comp_all, file.path(out_dir, "Completeness_Combined_PeptideProtein_byGroup.png"))
  }
  
  # ---- Peptides-per-protein (combined, + UNION)  [FIXED y-axis shared]
  ppp_res_list <- lapply(names(pr_raw_list), function(s) peptides_per_protein_stacked_by_group(pr_raw_list[[s]], group_map, s))
  names(ppp_res_list) <- names(pr_raw_list)
  
  ppp_bins   <- rbindlist(lapply(ppp_res_list, `[[`, "bins"), fill = TRUE)
  ppp_median <- rbindlist(lapply(ppp_res_list, `[[`, "median"), fill = TRUE)
  ppp_dist   <- rbindlist(lapply(ppp_res_list, `[[`, "dist"), fill = TRUE)
  
  ppp_union <- build_ppp_union_from_searches(pr_raw_list, group_map)
  
  ppp_bins2   <- rbind(ppp_bins,   ppp_union$bins,   fill = TRUE)
  ppp_median2 <- rbind(ppp_median, ppp_union$median, fill = TRUE)
  ppp_dist2   <- rbind(ppp_dist,   ppp_union$dist,   fill = TRUE)
  
  if (nrow(ppp_bins2) > 0) {
    fwrite(ppp_bins2,   file.path(out_dir, "PeptidesPerProtein_bins_byGroup.csv"))
    fwrite(ppp_median2, file.path(out_dir, "PeptidesPerProtein_median_byGroup.csv"))
    fwrite(ppp_dist2,   file.path(out_dir, "PeptidesPerProtein_distribution_perProtein.csv"))
    
    plot_peptides_per_protein_combined(ppp_bins2, file.path(out_dir, "PeptidesPerProtein_Combined_byGroup.png"))
    
  }
  
  plot_completeness_and_ppp_combined(
    comp_dt = comp_all,
    ppp_dt  = ppp_bins2,
    out_png = file.path(out_dir, "StackedBars_CompletenessTop_PeptidesPerProteinBottom.png")
  )
  
  
  
  # ---- Group scatters + R² violin (SELECT SEARCHES via interface)  [NEW]
  if (isTRUE(enable_group_plots)) {
    gm <- as.data.table(group_map)
    groups <- sort(unique(gm$Group))
    
    # Hardcoded, deterministic comparison plan (PI-specified). Falls back to the
    # interactive picker only if the hardcoded builder yields nothing.
    selected_pairs <- build_hardcoded_group_pairs(groups)
    if (nrow(selected_pairs) == 0) selected_pairs <- prompt_group_pair_selection(groups)
    cat(sprintf("\n[group-pairs] Hardcoded comparison plan (%d pairs):\n", nrow(selected_pairs)))
    if (nrow(selected_pairs) > 0)
      for (k in seq_len(nrow(selected_pairs)))
        cat(sprintf("   %-12s vs  %s\n", selected_pairs$g1[k], selected_pairs$g2[k]))

    rank_pairs_by_matched_keys <- function(sample_meta) {
      groups <- sort(unique(sample_meta$Group))
      pairs <- CJ(i = seq_along(groups), j = seq_along(groups))[i < j]
      pairs[, g1 := groups[i]]
      pairs[, g2 := groups[j]]
      pairs[, MatchedKeys := {
        pp <- get_matched_sample_pairs(sample_meta, g1, g2)
        length(unique(pp$SampleKey))
      }, by = .(g1, g2)]
      pairs[MatchedKeys > 0][order(-MatchedKeys)]
    }
    
    # choose searches once (applies to BOTH peptide and protein levels)
    # Hardcoded: every individual library search AND the union, so each comparison
    # is shown within each search and across the combined ID space (deterministic).
    search_choices <- c(names(search_paths), ALL_SEARCHES_LABEL)
    cat("\n[search-choices] ", paste(search_choices, collapse = ", "), "\n", sep = "")
    
    run_group_plots_for_level <- function(level_name, mats_named, level_dir, search_choices) {
      message("\n[", level_name, "] Generating selected group scatters + R²...")
      
      for (search_label in search_choices) {
        
        # choose matrix
        if (identical(search_label, ALL_SEARCHES_LABEL)) {
          m_any <- build_union_matrix(mats_named)
        } else {
          if (!search_label %in% names(mats_named)) {
            message("  Skipping missing search: ", search_label)
            next
          }
          m_any <- mats_named[[search_label]]
        }
        
        sample_meta <- make_sample_meta(colnames(m_any), group_map)
        
        # pick pairs
        pairs_ranked <- rank_pairs_by_matched_keys(sample_meta)
        if (nrow(selected_pairs) == 0) {
          pairs_use <- pairs_ranked[, .(g1, g2)][seq_len(min(top_group_pairs_scatter, .N))]
        } else {
          pairs_use <- copy(selected_pairs)
        }
        
        # output dirs
        comp_dir <- file.path(level_dir, paste0("GroupComparisons_", sanitize_filename(search_label)))
        dir.create(comp_dir, showWarnings = FALSE)
        dir.create(file.path(comp_dir, "Scatter_SelectedPairs"), showWarnings = FALSE)
        
        fwrite(sample_meta, file.path(comp_dir, "sample_meta_with_keys.csv"))
        fwrite(pairs_use, file.path(comp_dir, "selected_group_pairs.csv"))
        
        # scatters
        scatter_dir <- file.path(comp_dir, "Scatter_SelectedPairs")
        for (k in seq_len(nrow(pairs_use))) {
          g1 <- pairs_use$g1[k]; g2 <- pairs_use$g2[k]
          
          pts <- collect_scatter_points_all_ids(
            m_any, sample_meta, g1, g2,
            use_log2 = use_log2_for_cor,
            max_ids  = max_ids_for_scatter_calc
          )
          if (is.null(pts)) next
          
          out_png <- file.path(scatter_dir, paste0("scatter_", sanitize_filename(search_label), "_", g1, "_vs_", g2, ".png"))
          scatter_plot_with_stats(
            pts = pts,
            out_png = out_png,
            xlab = g1,
            ylab = g2,
            method = regression_method,
            max_points_plot = max_scatter_points_plot,
            bins = scatter_bins
          )
        }
        
        # R² violin
        r2_all <- rbindlist(lapply(seq_len(nrow(pairs_use)), function(k) {
          g1 <- pairs_use$g1[k]; g2 <- pairs_use$g2[k]
          r2dt <- compute_r2_per_samplepair(
            m_any, sample_meta, g1, g2,
            min_points = min_points_for_cor,
            use_log2   = use_log2_for_cor,
            max_ids    = max_ids_for_r2_calc
          )
          if (nrow(r2dt) == 0) return(NULL)
          r2dt[, `:=`(Search = search_label, Group1 = g1, Group2 = g2, Pair = paste0(g1, " vs ", g2))]
          r2dt
        }), fill = TRUE)
        
        if (nrow(r2_all) > 0) {
          fwrite(r2_all, file.path(comp_dir, "r2_per_matched_samplepair_selected.csv"))
          plot_r2_violin_group_pairs(
            r2_all,
            file.path(comp_dir, paste0("r2_violin_selected_", sanitize_filename(search_label), ".png"))
          )
        }
      }
    }
    
    run_group_plots_for_level("Peptide", peptide_mats, pep_res$level_dir, search_choices)
    run_group_plots_for_level("Protein", protein_mats, prot_res$level_dir, search_choices)
    df_search_hm <- make_df_search_heatmap(peptide_pres, protein_pres)
    # Combined UNION group overlap heatmap in same triangle style
    m_union_pep <- build_union_matrix(peptide_mats)
    m_union_pro <- build_union_matrix(protein_mats)
    df_group_hm <- make_df_group_heatmap_from_union(m_union_pep, m_union_pro, group_map)
    # smaller = fewer rows
    if (nlevels(df_search_hm$Col) * nlevels(df_search_hm$Row) <= nlevels(df_group_hm$Col) * nlevels(df_group_hm$Row)) {
      df_left  <- df_search_hm
      df_right <- df_group_hm
    } else {
      df_left  <- df_group_hm
      df_right <- df_search_hm
    }
    plot_group_overlap_combined_from_union(
      m_union_pep, m_union_pro, group_map,
      file.path(out_dir, "UnionGroupOverlap_Proteins_vs_Peptides_Triangles.png")
    )
    p_search_hm <- plot_search_overlap_combined(peptide_pres, protein_pres, out_png = NULL, return_plot = TRUE)
    p_group_hm  <- plot_group_overlap_combined_from_union(m_union_pep, m_union_pro, group_map, out_png = NULL, return_plot = TRUE)
    p_search_hm <- p_search_hm + shared_fill
    p_group_hm  <- p_group_hm  + shared_fill
    p_search_hm <- p_search_hm +
      coord_fixed(ratio = 0.75) +   # <— shorter rectangles
      theme(
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.height = unit(0.3, "cm"),
        legend.key.width  = unit(1.4, "cm")
      )
    
    p_group_hm <- p_group_hm +
      coord_fixed(ratio = 0.75) +
      theme(legend.position = "none")   # legend only once
    
    # df_left = smaller heatmap (fewer rows)
    # df_right = larger heatmap
    
    
    plot_two_heatmaps_side_by_side_legend_over_left(
      df_left  = df_left,
      df_right = df_right,
      out_png  = file.path(out_dir, "Heatmaps_SideBySide_LegendOverLeft.png"),
      width = 16.5,
      height = 5.0,
      text_size = 4.2,
      cell_ratio = 0.38,
      left_height_rel = 0.85
    )
    
  }
  
  # ==============================
  # INDEX
  # ==============================
  sections <- list(
    "Run-level tables & maps" = list(
      list("search_name_mapping.csv", "search_name_mapping.csv"),
      list("sample_counts_per_search.csv", "sample_counts_per_search.csv"),
      list("sample_list_common.csv", "sample_list_common.csv"),
      list("sample_group_mapping.csv", "sample_group_mapping.csv")
    ),
    
    "Search overlap" = list(
      list("Search overlap (proteins vs peptides triangles)",
           "SearchOverlap_Proteins_vs_Peptides_Triangles.png")
    ),
    
    "IDs by group" = list(
      list("IDs by group violin (protein + peptide combined)",
           "IDs_by_group_Violin_Combined_Peptide_Protein.png")
    ),
    
    "Completeness & peptides-per-protein" = list(
      list("Stacked bars combined (completeness top / peptides-per-protein bottom)",
           "StackedBars_CompletenessTop_PeptidesPerProteinBottom.png"),
      list("Completeness (peptide+protein combined, by group, + Union)",
           "Completeness_Combined_PeptideProtein_byGroup.png"),
      list("Peptides per protein (by group, + Union)",
           "PeptidesPerProtein_Combined_byGroup.png"),
      list("Union group overlap (proteins vs peptides triangles)",
           "UnionGroupOverlap_Proteins_vs_Peptides_Triangles.png")
    ),
    
    "Peptide level" = list(
      list("IDs by group violin", file.path("Peptide", "ids_by_group_violin.png")),
      list("Selected scatters (Union)",
           file.path("Peptide",
                     paste0("GroupComparisons_", sanitize_filename(ALL_SEARCHES_LABEL)),
                     "Scatter_SelectedPairs")),
      list("Selected R² violin (Union)",
           file.path("Peptide",
                     paste0("GroupComparisons_", sanitize_filename(ALL_SEARCHES_LABEL)),
                     paste0("r2_violin_selected_", sanitize_filename(ALL_SEARCHES_LABEL), ".png")))
    ),
    
    "Protein level" = list(
      list("IDs by group violin", file.path("Protein", "ids_by_group_violin.png")),
      list("Selected scatters (Union)",
           file.path("Protein",
                     paste0("GroupComparisons_", sanitize_filename(ALL_SEARCHES_LABEL)),
                     "Scatter_SelectedPairs")),
      list("Selected R² violin (Union)",
           file.path("Protein",
                     paste0("GroupComparisons_", sanitize_filename(ALL_SEARCHES_LABEL)),
                     paste0("r2_violin_selected_", sanitize_filename(ALL_SEARCHES_LABEL), ".png")))
    )
  )
  
  
  index_path <- file.path(out_dir, "index.html")
  write_index_html(index_path, sections)
  
  message("\nDONE.")
  message("Report folder: ", out_dir)
  message("Open this file in your browser: ", index_path)
  cat("\nDONE.\nReport folder: ", out_dir, "\nIndex: ", index_path, "\n", sep = "")
  
}, error = function(e) {
  msg <- paste0("ERROR: ", conditionMessage(e))
  cat("\n", msg, "\n", sep = "")
  writeLines(msg, file.path(out_dir, "ERROR.txt"))
  stop(e)
})


