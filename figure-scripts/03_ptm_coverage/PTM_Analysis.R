# --------------------------------------------------------------------------------------------------
# Two-Protease-Paper | PTM report from the discovery searches
# Manuscript: Fig 3A (PTM sites per sample, PTM-site overlap triangle) and the PTM counts in the text.
# Authors: Simion Kreimer (original); Mario M. Alba (manuscript group labels and comparison plan).
# Inputs: TWO_PROTEASE_DATA/<search>/report.pr_matrix.tsv; TWO_PROTEASE_FASTA for site numbering;
#         unimod_name_library.tsv (next to this script). Downloads UniMod names from unimod.org if needed.
# Outputs: TWO_PROTEASE_OUT/03_ptm_coverage/DIANN_PTM_Report_<timestamp>/
# Run: Rscript figure-scripts/03_ptm_coverage/PTM_Analysis.R (prompts take defaults).
# Packages: data.table, ggplot2, viridisLite, ggnewscale, Biostrings (Bioconductor).
# Released as used for the manuscript: only file paths, line endings and this header were changed;
# the analysis code is untouched. Paths come from environment variables, see ../../README.md.
# --------------------------------------------------------------------------------------------------

# ============================================================
# DIA-NN multi-search PTM report from report.pr_matrix.tsv
# - PTM stacked bars (Top N, UniMod name + ID)
# - PTM-site mapping to proteins (FASTA) for MODIFIED peptides only
# - Group overlap triangle heatmap:
#     upper triangle = intersection size (tile color + numbers)
#     lower triangle = PTM makeup barcode (Top N mods ONLY, same colors as stacked bars)
#     NO "Other"
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Biostrings)
})

# ------------------------------
# USER SETTINGS
# ------------------------------
pg_file <- "report.pg_matrix.tsv"      # used only for folder discovery
pr_file <- "report.pr_matrix.tsv"

peptide_collapse_method <- "sum"       # "sum" or "mean"
zero_as_na <- TRUE

exclude_unimod_ids <- c("4","35","2","7","385","23","28","30","6","193","421","5")  # excluded from PTM counting + site mapping
TOP_N_PTMS <- 20
max_searches <- 100

PUB_DPI <- 600

# ------------------------------
# THEME + SAVE
# ------------------------------
pub_theme <- function(base_size = 14, strip_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", color = "black", size = base_size + 2),
      axis.title = element_text(face = "bold", color = "black", size = base_size + 1),
      axis.text  = element_text(color = "black", face = "plain", size = base_size),
      legend.text  = element_text(color = "black", size = base_size - 1),
      strip.text   = element_text(color = "black", face = "bold", size = strip_size),
      axis.text.x  = element_text(angle = 0, hjust = 0.5),
      panel.grid   = element_blank(),
      strip.background = element_blank()
    )
}

save_pub <- function(plot_obj, out_png, width, height, dpi = PUB_DPI) {
  ggsave(out_png, plot_obj, width = width, height = height, dpi = dpi, units = "in")
  out_pdf <- sub("\\.png$", ".pdf", out_png, ignore.case = TRUE)
  ggsave(out_pdf, plot_obj, width = width, height = height, device = cairo_pdf, units = "in")
}

# ------------------------------
# WINDOWS PICKERS
# ------------------------------
pick_parent_folder <- function() {
  if (interactive() && .Platform$OS.type == "windows") {
    p <- utils::choose.dir(caption = "Select parent folder containing DIA-NN search result folders")
    if (is.na(p) || !nzchar(p)) stop("No folder selected.")
    return(p)
  }
  stop("Folder picker is Windows interactive only. Set base_dir manually.")
}

pick_fasta_file <- function() {
  if (interactive() && .Platform$OS.type == "windows") {
    f <- utils::choose.files(caption = "Select FASTA used for DIA-NN search", multi = FALSE)
    if (length(f) == 0 || is.na(f) || !nzchar(f)) stop("No FASTA selected.")
    return(f)
  }
  stop("FASTA picker is Windows interactive only. Set fasta_path manually.")
}

# ------------------------------
# BASIC HELPERS
# ------------------------------
stop_if_missing <- function(path) if (!file.exists(path)) stop("Missing file: ", path)
# Given site_dt, return a named numeric vector of counts by ModLabel for a set of SiteKeys
region_mod_counts <- function(site_dt, sitekeys, mod_keep) {
  dt <- as.data.table(site_dt)
  dt <- dt[SiteKey %chin% sitekeys]
  dt <- dt[ModLabel %chin% mod_keep]
  if (nrow(dt) == 0) return(setNames(numeric(0), character(0)))
  v <- dt[, .N, by = ModLabel][order(ModLabel)]
  setNames(v$N, v$ModLabel)
}

# Turn counts into stacked rectangles for a circular region (simple horizontal bars clipped to a circle)
make_circle_barcodes <- function(counts, cols, cx, cy, r, n_stripes = 40) {
  
  # counts: named numeric vector
  if (is.null(counts) || length(counts) == 0) return(data.table())
  
  counts <- as.numeric(counts)
  names(counts) <- names(counts)
  counts <- counts[!is.na(counts)]
  counts[is.na(counts)] <- 0
  
  counts[is.na(counts)] <- 0
  counts <- counts[is.finite(counts) & counts > 0]
  
  if (length(counts) == 0) return(data.table())
  
  tot <- sum(counts)
  if (!is.finite(tot) || tot <= 0) return(data.table())
  
  props <- counts / tot
  props <- props[is.finite(props) & props > 0]
  if (length(props) == 0) return(data.table())
  
  # allocate stripes
  times <- pmax(1L, as.integer(round(props * n_stripes)))
  if (any(!is.finite(times)) || any(times < 1)) return(data.table())
  
  lab <- rep(names(props), times)
  if (length(lab) == 0) return(data.table())
  
  lab <- lab[seq_len(min(length(lab), n_stripes))]
  
  # build stripe rectangles
  x0 <- cx - r; x1 <- cx + r
  y0 <- cy - r; y1 <- cy + r
  h <- (y1 - y0) / length(lab)
  
  data.table(
    xmin = x0, xmax = x1,
    ymin = y0 + (seq_along(lab) - 1) * h,
    ymax = y0 + seq_along(lab) * h,
    ModLabel = lab,
    fill = cols[lab]
  )
}


# ---- geometry: circle intersection area ----
.circle_intersection_area <- function(r1, r2, d) {
  # no overlap
  if (d >= r1 + r2) return(0)
  # one inside the other
  if (d <= abs(r1 - r2)) return(pi * min(r1, r2)^2)
  
  # general case
  a1 <- r1^2 * acos((d^2 + r1^2 - r2^2) / (2*d*r1))
  a2 <- r2^2 * acos((d^2 + r2^2 - r1^2) / (2*d*r2))
  a3 <- 0.5 * sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
  a1 + a2 - a3
}

.solve_distance_for_overlap <- function(r1, r2, target_area) {
  maxA <- pi * min(r1, r2)^2
  target_area <- max(0, min(maxA, target_area))
  
  if (target_area <= 0) return(r1 + r2)              # no overlap
  if (abs(target_area - maxA) < 1e-9) return(abs(r1 - r2))  # full containment overlap
  
  f <- function(d) .circle_intersection_area(r1, r2, d) - target_area
  lo <- abs(r1 - r2) + 1e-9
  hi <- (r1 + r2) - 1e-9
  
  # monotonic decreasing in d, so uniroot works
  uniroot(f, lower = lo, upper = hi)$root
}

# ---- PTM composition (by ModLabel, matching your palette names) ----
.region_mod_counts <- function(site_dt, sitekeys, mod_keep) {
  dt <- as.data.table(site_dt)[SiteKey %chin% sitekeys]
  dt <- dt[ModLabel %chin% mod_keep]
  if (nrow(dt) == 0) return(setNames(numeric(0), character(0)))
  v <- dt[, .N, by = ModLabel]
  setNames(v$N, v$ModLabel)
}

# ---- make vertical "barcode" stripes in a rectangle ----
.make_barcode_rects <- function(counts_named, cols_named, x0, x1, y0, y1, n_stripes = 40) {
  if (is.null(counts_named) || length(counts_named) == 0) return(data.table())
  
  counts <- counts_named
  counts[is.na(counts)] <- 0
  counts <- counts[is.finite(counts) & counts > 0]
  if (length(counts) == 0) return(data.table())
  
  tot <- sum(counts)
  if (!is.finite(tot) || tot <= 0) return(data.table())
  
  props <- counts / tot
  props[is.na(props)] <- 0
  props <- props[is.finite(props) & props > 0]
  if (length(props) == 0) return(data.table())
  
  times <- pmax(1L, as.integer(round(props * n_stripes)))
  times[!is.finite(times) | is.na(times) | times < 1] <- 1L
  
  lab <- rep(names(props), times)
  if (length(lab) == 0) return(data.table())
  lab <- lab[seq_len(min(length(lab), n_stripes))]
  
  w <- (x1 - x0) / length(lab)
  
  data.table(
    xmin = x0 + (seq_along(lab) - 1) * w,
    xmax = x0 + seq_along(lab) * w,
    ymin = y0,
    ymax = y1,
    ModLabel = lab,
    fill = cols_named[lab]
  )
}

# ---- proportional 2-set Venn with PTM barcodes for A, B, and overlap ----
plot_venn2_ptmfilled_proportional <- function(site_dt, site_sets, mod_keep, mod_cols, out_png,
                                              title = "PTM-site overlap by search strategy (Top20 PTMs)") {
  stopifnot(length(site_sets) == 2)
  nm <- names(site_sets)
  A <- unique(site_sets[[1]])
  B <- unique(site_sets[[2]])
  
  I <- intersect(A, B)
  nA <- length(A); nB <- length(B); nI <- length(I)
  nAonly <- length(setdiff(A, I))
  nBonly <- length(setdiff(B, I))
  
  # --- radii proportional to sqrt(n)
  rA <- sqrt(nA)
  rB <- sqrt(nB)
  # scale radii to a nice drawing size
  scale <- 1.2 / max(rA, rB)
  r1 <- rA * scale
  r2 <- rB * scale
  
  # --- target overlap area fraction ~ nI / min(nA,nB)
  frac <- if (min(nA, nB) > 0) (nI / min(nA, nB)) else 0
  frac <- max(0, min(1, frac))
  target_area <- frac * (pi * min(r1, r2)^2)
  
  d <- .solve_distance_for_overlap(r1, r2, target_area)
  
  # centers
  c1x <- -d/2
  c2x <-  d/2
  c1y <- 0
  c2y <- 0
  
  circle_pts <- function(cx, cy, r, n = 500) {
    th <- seq(0, 2*pi, length.out = n)
    data.table(x = cx + r*cos(th), y = cy + r*sin(th))
  }
  c1 <- circle_pts(c1x, c1y, r1)
  c2 <- circle_pts(c2x, c2y, r2)
  
  # --- use exact palette names
  cols <- mod_cols
  mod_keep2 <- intersect(mod_keep, names(cols))
  
  # compositions for full A, full B, and overlap I (NOT A-only/B-only)
  cA <- .region_mod_counts(site_dt, A, mod_keep2)
  cB <- .region_mod_counts(site_dt, B, mod_keep2)
  cI <- .region_mod_counts(site_dt, I, mod_keep2)
  
  # barcode rectangles:
  # inside each circle: a tall thin rectangle centered in each circle
  w1 <- 0.55 * r1; h1 <- 1.55 * r1
  w2 <- 0.55 * r2; h2 <- 1.55 * r2
  
  fillA <- .make_barcode_rects(cA[mod_keep2], cols, c1x - w1/2, c1x + w1/2, -h1/2, h1/2, n_stripes = 40)
  fillB <- .make_barcode_rects(cB[mod_keep2], cols, c2x - w2/2, c2x + w2/2, -h2/2, h2/2, n_stripes = 40)
  
  # overlap barcode: rectangle centered at 0
  wI <- 0.65 * min(r1, r2)
  hI <- 1.60 * min(r1, r2)
  fillI <- .make_barcode_rects(cI[mod_keep2], cols, -wI/2, wI/2, -hI/2, hI/2, n_stripes = 40)
  
  # label positions
  x_left_num  <- c1x - 0.55*r1
  x_right_num <- c2x + 0.55*r2
  x_mid_num   <- 0
  
  # plot limits
  xlim <- c(min(c1x - r1, c2x - r2) - 0.8, max(c1x + r1, c2x + r2) + 0.8)
  ylim <- c(-max(r1, r2) - 0.8, max(r1, r2) + 0.9)
  
  png(out_png, width = 2000, height = 1400, res = 200)
  on.exit(dev.off(), add = TRUE)
  
  p <- ggplot() +
    # barcodes first (so outlines/labels sit on top)
    { if (nrow(fillA) > 0) geom_rect(data = fillA, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
                                     fill = fillA$fill, alpha = 0.35, color = NA) } +
    { if (nrow(fillB) > 0) geom_rect(data = fillB, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
                                     fill = fillB$fill, alpha = 0.35, color = NA) } +
    { if (nrow(fillI) > 0) geom_rect(data = fillI, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
                                     fill = fillI$fill, alpha = 0.60, color = NA) } +
    # circle outlines
    geom_path(data = c1, aes(x, y), linewidth = 1.0, color = "black") +
    geom_path(data = c2, aes(x, y), linewidth = 1.0, color = "black") +
    # numbers (bold)
    annotate("text", x = x_left_num,  y = 0, label = as.character(nAonly), fontface = "bold", size = 12) +
    annotate("text", x = x_mid_num,   y = 0, label = as.character(nI),     fontface = "bold", size = 12) +
    annotate("text", x = x_right_num, y = 0, label = as.character(nBonly), fontface = "bold", size = 12) +
    # set labels
    annotate("text", x = c1x, y = max(r1, r2) + 0.35, label = paste0(nm[1], "\n(n=", nA, ")"),
             fontface = "bold", size = 7) +
    annotate("text", x = c2x, y = max(r1, r2) + 0.35, label = paste0(nm[2], "\n(n=", nB, ")"),
             fontface = "bold", size = 7) +
    labs(
      title = title,
      subtitle = paste0("Intersection=", nI,
                        " | ", nm[1], " only=", nAonly,
                        " | ", nm[2], " only=", nBonly,
                        " | circles ~ set size | barcodes = Top20 PTM makeup")
    ) +
    coord_equal(xlim = xlim, ylim = ylim, expand = FALSE) +
    theme_void(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  print(p)
  invisible(TRUE)
}


write_stackedbar_order_csv <- function(group_map, out_samples_csv, out_keys_csv) {
  gm <- as.data.table(group_map)
  
  # Ensure Group exists for all rows (don’t drop anything)
  gm[, Group := as.character(Group)]
  gm[is.na(Group) | !nzchar(Group), Group := "Unknown"]
  
  gm[, SampleKey := vapply(Sample, extract_sample_key, character(1))]
  gm[, SampleKey := as.character(SampleKey)]
  
  sample_meta <- unique(gm[, .(Sample, Group, SampleKey)])
  
  # EXACT same ordering logic used by the plot
  key_levels <- make_matched_samplekey_levels(sample_meta, group_col = "Group", key_col = "SampleKey")
  
  # Key-level order table
  key_dt <- data.table(
    KeyOrder = match(key_levels, key_levels),
    SampleKey = key_levels
  )
  
  # Sample-level table (this is the one you care about for “324 samples”)
  sample_dt <- merge(sample_meta, key_dt, by = "SampleKey", all.x = TRUE)
  
  # Also include position within each facet (same as KeyOrder, repeated)
  sample_dt[, FacetPosition := KeyOrder]
  
  # Sort as they appear left->right within each Group facet
  data.table::setorder(sample_dt, Group, FacetPosition, Sample)
  
  fwrite(sample_dt, out_samples_csv)
  fwrite(key_dt, out_keys_csv)
}


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

sanitize_ids <- function(ids) {
  ids <- as.character(ids)
  ids <- trimws(ids)
  ids[ids %in% c("", "NA", "NaN")] <- NA_character_
  ids
}

make_unique_names <- function(x) {
  nm <- names(x)
  if (anyDuplicated(nm)) names(x) <- make.unique(nm, sep = "_")
  x
}
# Build strategy sets using UnimodID (never mismatches label text)
build_strategy_site_sets_by_unimodid <- function(site_dt, unimod_keep) {
  dt <- as.data.table(site_dt)
  dt[, UnimodID := as.character(UnimodID)]
  unimod_keep <- as.character(unique(unimod_keep))
  
  dt <- dt[UnimodID %chin% unimod_keep]
  dt <- unique(dt[, .(Search, SiteKey)])
  
  keep <- dt[, .N, by = Search][N > 0, Search]
  sets <- lapply(keep, function(s) unique(dt[Search == s, SiteKey]))
  names(sets) <- keep
  sets
}

plot_venn2_gg <- function(A, B, nameA, nameB, out_png, title = "PTM-site overlap (Top20)") {
  # A, B are vectors of SiteKey
  A <- unique(A); B <- unique(B)
  nA <- length(A)
  nB <- length(B)
  nI <- length(intersect(A, B))
  nAonly <- nA - nI
  nBonly <- nB - nI
  
  # circle points
  circle_pts <- function(cx, cy, r, n = 400) {
    th <- seq(0, 2*pi, length.out = n)
    data.table(x = cx + r*cos(th), y = cy + r*sin(th))
  }
  
  r <- 1
  c1 <- circle_pts(-0.6, 0, r)
  c2 <- circle_pts( 0.6, 0, r)
  
  # label positions tuned for readability
  dt_txt <- data.table(
    x = c(-1.2, 0, 1.2),
    y = c(0, 0, 0),
    lab = c(as.character(nAonly), as.character(nI), as.character(nBonly)),
    font = 2
  )
  
  # set labels
  dt_sets <- data.table(
    x = c(-0.6, 0.6),
    y = c(1.25, 1.25),
    lab = c(paste0(nameA, "\n(n=", nA, ")"), paste0(nameB, "\n(n=", nB, ")"))
  )
  
  png(out_png, width = 1600, height = 1200, res = 200)
  on.exit(dev.off(), add = TRUE)
  
  p <- ggplot() +
    geom_path(data = c1, aes(x, y), linewidth = 2) +
    geom_path(data = c2, aes(x, y), linewidth = 2) +
    annotate("text", x = dt_txt$x, y = dt_txt$y, label = dt_txt$lab, fontface = "bold", size = 10) +
    annotate("text", x = dt_sets$x, y = dt_sets$y, label = dt_sets$lab, fontface = "bold", size = 6) +
    labs(title = title, subtitle = paste0("Intersection=", nI, " | ", nameA, " only=", nAonly, " | ", nameB, " only=", nBonly)) +
    coord_equal(xlim = c(-2.2, 2.2), ylim = c(-1.6, 1.8), expand = FALSE) +
    theme_void(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  print(p)
  invisible(TRUE)
}

# ------------------------------
# SEARCH DISCOVERY + SHORT NAMES
# ------------------------------
discover_search_folders <- function(parent_dir, pg_file = "report.pg_matrix.tsv", pr_file = "report.pr_matrix.tsv") {
  subs <- list.dirs(parent_dir, full.names = TRUE, recursive = FALSE)
  ok <- vapply(subs, function(d) {
    file.exists(file.path(d, pg_file)) && file.exists(file.path(d, pr_file))
  }, logical(1))
  found <- subs[ok]
  names(found) <- basename(found)
  found <- found[order(names(found))]
  found
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

# ------------------------------
# SAMPLE COLUMN DETECTION
# ------------------------------
normalize_colname <- function(x) tolower(gsub("[^a-z0-9]+", "", x))

is_non_sample_col <- function(colname) {
  nrm <- normalize_colname(colname)
  if (nrm %in% c(
    "proteotypic",
    "nsequences","nsequence","nseq","nseqs",
    "precursorcharge","precursorchargeid",
    "qvalue","pgqvalue","pep","pgpep","libqvalue",
    "decoy","decoyevidence",
    "rt","rtstart","rtstop",
    "ms1area","ms1profilecorr"
  )) return(TRUE)
  FALSE
}

looks_like_sample_colname <- function(colname) {
  cn <- as.character(colname)
  grepl("\\.(raw|d|mzml)$", cn, ignore.case = TRUE) ||
    grepl("^[A-Za-z]:\\\\", cn) ||
    grepl("^[A-Za-z]:/", cn) ||
    grepl("[/\\\\].*\\.(raw|d|mzml)$", cn, ignore.case = TRUE)
}

numeric_cols_only <- function(dt, exclude = character()) {
  dt <- as.data.table(dt)
  cols_all <- setdiff(names(dt), exclude)
  cols_all <- cols_all[!vapply(cols_all, is_non_sample_col, logical(1))]
  
  sample_like <- cols_all[vapply(cols_all, looks_like_sample_colname, logical(1))]
  if (length(sample_like) > 0) return(sample_like)
  
  keep <- character()
  for (cname in cols_all) {
    v <- suppressWarnings(as.numeric(dt[[cname]]))
    if (sum(!is.na(v)) > 0) keep <- c(keep, cname)
  }
  unique(keep)
}

# ------------------------------
# GROUPING (your heuristics)
# ------------------------------
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
# Mirrors apply_manuscript_group_labels() in Identification_Analysis.R so PTM
# figures (Fig 3) carry the same nomenclature: Try_<date> -> Trypsin 1/2,
# Kra_<date> -> Krakatoa 1/2 (earliest date = set 1), Mix_20min/30min -> Mix
# 20min/30min. Set rank derived from sorted distinct dates (no hardcoded dates).
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

# de-identifying x-axis key (keep your BRH logic)
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


make_matched_samplekey_levels <- function(sample_meta, group_col = "Group", key_col = "SampleKey") {
  sm <- as.data.table(sample_meta)
  sm <- sm[!is.na(get(key_col)) & nzchar(get(key_col)) &
             !is.na(get(group_col)) & nzchar(get(group_col))]
  if (nrow(sm) == 0) return(character())
  keys_by_group <- split(sm[[key_col]], sm[[group_col]])
  keys_by_group <- lapply(keys_by_group, unique)
  matched <- if (length(keys_by_group) >= 2) Reduce(intersect, keys_by_group) else unique(sm[[key_col]])
  matched <- sort(unique(matched))
  remaining <- setdiff(sort(unique(sm[[key_col]])), matched)
  c(matched, remaining)
}

# ------------------------------
# PTM parsing
# ------------------------------
extract_unimod_ids <- function(mod_seq, exclude_ids = c("4")) {
  if (is.na(mod_seq) || !nzchar(mod_seq)) return(character())
  hits <- gregexpr("UniMod:\\s*([0-9]+)", mod_seq, perl = TRUE)
  m <- regmatches(mod_seq, hits)[[1]]
  if (length(m) == 0) return(character())
  ids <- gsub(".*UniMod:\\s*([0-9]+).*", "\\1", m, perl = TRUE)
  ids <- ids[nzchar(ids)]
  ids <- ids[!(ids %in% as.character(exclude_ids))]
  ids
}

# Collapse pr_matrix rows -> peptide matrix by Modified.Sequence
collapse_pr_to_peptides_by_modified_sequence <- function(pr_dt,
                                                         collapse_method = "sum",
                                                         zero_as_na = TRUE,
                                                         id_col = "Modified.Sequence") {
  pr_dt <- as.data.table(pr_dt)
  if (!id_col %in% names(pr_dt)) {
    cand <- c("Modified.Sequence", "ModifiedPeptideSequence", "Modified.Peptide.Sequence", "ModifiedPeptide.Sequence")
    cand <- cand[cand %in% names(pr_dt)]
    if (length(cand) == 0) stop("pr_matrix missing Modified.Sequence-like column.")
    id_col <- cand[1]
  }
  
  pr_dt[[id_col]] <- sanitize_ids(pr_dt[[id_col]])
  pr_dt <- pr_dt[!is.na(pr_dt[[id_col]]) & nzchar(pr_dt[[id_col]])]
  
  num_cols <- numeric_cols_only(pr_dt, exclude = id_col)
  if (length(num_cols) == 0) stop("No numeric/sample-like columns detected in pr_matrix.")
  
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

# align sample columns across matrices (union)
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

# union matrix across searches = average across available
align_ids_union <- function(mats_named) {
  all_ids <- sort(unique(unlist(lapply(mats_named, rownames))))
  lapply(mats_named, function(m) {
    miss <- setdiff(all_ids, rownames(m))
    if (length(miss) > 0) {
      add <- matrix(NA_real_, nrow = length(miss), ncol = ncol(m))
      rownames(add) <- miss
      colnames(add) <- colnames(m)
      m <- rbind(m, add)
    }
    m[all_ids, , drop = FALSE]
  })
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

# ------------------------------
# UniMod name library (fast cache)
# ------------------------------
html_unescape_basic <- function(x) {
  x <- gsub("&nbsp;", " ", x, fixed = TRUE)
  x <- gsub("&amp;", "&", x, fixed = TRUE)
  x <- gsub("&lt;", "<", x, fixed = TRUE)
  x <- gsub("&gt;", ">", x, fixed = TRUE)
  x <- gsub("&#39;", "'", x, fixed = TRUE)
  x <- gsub("&quot;", "\"", x, fixed = TRUE)
  x
}

fetch_unimod_name_one <- function(unimod_id) {
  url <- paste0("https://www.unimod.org/modifications_view.php?editid1=", unimod_id)
  tmp <- tempfile(fileext = ".html")
  ok <- tryCatch({
    utils::download.file(url, destfile = tmp, quiet = TRUE, mode = "wb")
    TRUE
  }, error = function(e) FALSE, warning = function(w) TRUE)
  if (!ok || !file.exists(tmp)) return(NA_character_)
  
  x <- tryCatch(readLines(tmp, warn = FALSE, encoding = "UTF-8"), error = function(e) character())
  if (length(x) == 0) return(NA_character_)
  txt <- paste(x, collapse = " ")
  txt <- gsub("[\r\n\t ]+", " ", txt)
  txt <- html_unescape_basic(txt)
  
  grab_td_value <- function(label_regex) {
    pat <- paste0("(?i)", label_regex, "\\s*</td>\\s*<td[^>]*>\\s*([^<]+?)\\s*</td>")
    m <- regexpr(pat, txt, perl = TRUE)
    if (m[1] == -1) return(NA_character_)
    val <- sub(pat, "\\1", regmatches(txt, m), perl = TRUE)
    val <- trimws(val)
    if (!nzchar(val)) NA_character_ else val
  }
  
  v <- grab_td_value("PSI\\-MS\\s*Name"); if (!is.na(v)) return(v)
  v <- grab_td_value("Interim\\s*Name");  if (!is.na(v)) return(v)
  v <- grab_td_value("^Name$");           if (!is.na(v)) return(v)
  
  mt <- regexpr("(?i)<title>\\s*([^<]+?)\\s*</title>", txt, perl = TRUE)
  if (mt[1] != -1) {
    tval <- sub("(?i)<title>\\s*([^<]+?)\\s*</title>", "\\1", regmatches(txt, mt), perl = TRUE)
    tval <- trimws(tval)
    tval <- sub("(?i)^unimod\\s*[:-]\\s*", "", tval, perl = TRUE)
    if (nzchar(tval)) return(tval)
  }
  NA_character_
}

# persistent library file: only fetch IDs not already in library
load_unimod_library <- function(lib_tsv) {
  if (file.exists(lib_tsv)) {
    dt <- tryCatch(fread(lib_tsv, sep = "\t"), error = function(e) data.table())
    dt <- as.data.table(dt)
    if (nrow(dt) > 0 && all(c("UnimodID","UnimodName") %in% names(dt))) {
      dt[, UnimodID := as.character(UnimodID)]
      dt[, UnimodName := as.character(UnimodName)]
      dt <- unique(dt, by = "UnimodID")
      return(dt)
    }
  }
  data.table(UnimodID = character(), UnimodName = character())
}

write_unimod_library <- function(dt, lib_tsv) {
  dt <- as.data.table(dt)
  dt[, UnimodID := as.character(UnimodID)]
  dt[, UnimodName := as.character(UnimodName)]
  dt <- unique(dt, by = "UnimodID")
  dt[, UnimodID_int := suppressWarnings(as.integer(UnimodID))]
  setorder(dt, UnimodID_int, na.last = TRUE)
  dt[, UnimodID_int := NULL]
  fwrite(dt, lib_tsv, sep = "\t")
}

get_unimod_name_map_fast <- function(unimod_ids, lib_tsv) {
  unimod_ids <- sort(unique(as.character(unimod_ids)))
  unimod_ids <- unimod_ids[nzchar(unimod_ids)]
  if (length(unimod_ids) == 0) return(data.table(UnimodID = character(), UnimodName = character()))
  
  lib <- load_unimod_library(lib_tsv)
  missing <- setdiff(unimod_ids, lib$UnimodID)
  
  if (length(missing) > 0) {
    message("Translating UniMod IDs via unimod.org (new IDs: ", length(missing), ") ...")
    new <- data.table(UnimodID = missing, UnimodName = NA_character_)
    for (i in seq_along(missing)) {
      id <- missing[i]
      nm <- fetch_unimod_name_one(id)
      new[UnimodID == id, UnimodName := nm]
      if (i %% 25 == 0) message("  translated ", i, "/", length(missing))
    }
    lib <- unique(rbind(lib, new, fill = TRUE), by = "UnimodID")
    write_unimod_library(lib, lib_tsv)
  }
  
  lib[UnimodID %in% unimod_ids]
}

# ------------------------------
# PTM counts (from UNION peptide presence)
# ------------------------------
build_ptm_counts_from_union_peptides <- function(m_union_pep, exclude_ids) {
  pep_ids <- rownames(m_union_pep)
  pep_mods <- lapply(pep_ids, extract_unimod_ids, exclude_ids = exclude_ids)
  
  samples <- colnames(m_union_pep)
  out <- vector("list", length(samples))
  
  for (i in seq_along(samples)) {
    s <- samples[i]
    present <- which(!is.na(m_union_pep[, s]))
    if (length(present) == 0) {
      out[[i]] <- data.table(Sample = s, UnimodID = character(), Count = integer())
      next
    }
    mods <- unlist(pep_mods[present], use.names = FALSE)
    mods <- mods[nzchar(mods)]
    if (length(mods) == 0) {
      out[[i]] <- data.table(Sample = s, UnimodID = character(), Count = integer())
      next
    }
    tab <- as.data.table(table(mods))
    setnames(tab, c("UnimodID", "Count"))
    tab[, `:=`(Sample = s, Count = as.integer(Count))]
    out[[i]] <- tab[, .(Sample, UnimodID, Count)]
  }
  
  rbindlist(out, fill = TRUE)
}



# Count only a specific set of UniMod IDs.
# This is used to make separate plots for modifications that are excluded
# from the main PTM analysis, without changing the original main plots.
build_ptm_counts_for_unimod_ids <- function(m_union_pep, keep_ids) {
  keep_ids <- as.character(unique(keep_ids))
  keep_ids <- keep_ids[nzchar(keep_ids)]

  pep_ids <- rownames(m_union_pep)
  # Extract all UniMod IDs first, then keep only the requested excluded IDs.
  pep_mods <- lapply(pep_ids, extract_unimod_ids, exclude_ids = character())
  pep_mods <- lapply(pep_mods, function(x) x[x %chin% keep_ids])

  samples <- colnames(m_union_pep)
  out <- vector("list", length(samples))

  for (i in seq_along(samples)) {
    s <- samples[i]
    present <- which(!is.na(m_union_pep[, s]))
    if (length(present) == 0) {
      out[[i]] <- data.table(Sample = s, UnimodID = character(), Count = integer())
      next
    }

    mods <- unlist(pep_mods[present], use.names = FALSE)
    mods <- mods[nzchar(mods)]
    mods <- mods[mods %chin% keep_ids]

    if (length(mods) == 0) {
      out[[i]] <- data.table(Sample = s, UnimodID = character(), Count = integer())
      next
    }

    tab <- as.data.table(table(mods))
    setnames(tab, c("UnimodID", "Count"))
    tab[, `:=`(Sample = s, Count = as.integer(Count))]
    out[[i]] <- tab[, .(Sample, UnimodID, Count)]
  }

  rbindlist(out, fill = TRUE)
}


# Convert each search-specific peptide matrix into long peptide/sample/search presence table.
build_presence_dt_with_search <- function(peptide_mats) {
  out <- list()

  for (s in names(peptide_mats)) {
    m <- peptide_mats[[s]]
    if (is.null(m) || nrow(m) == 0 || ncol(m) == 0) next

    idx <- which(!is.na(m), arr.ind = TRUE)
    if (nrow(idx) == 0) next

    out[[length(out) + 1]] <- data.table(
      Search = s,
      Modified.Sequence = rownames(m)[idx[, "row"]],
      Sample = colnames(m)[idx[, "col"]]
    )
  }

  if (length(out) == 0) {
    return(data.table(
      Search = character(),
      Modified.Sequence = character(),
      Sample = character()
    ))
  }

  unique(rbindlist(out, fill = TRUE))
}

# compute TopN labels + exact colors (this is the single source of truth)
get_top_mod_palette <- function(ptm_dt, unimod_map, top_n = 20) {
  if (!requireNamespace("viridisLite", quietly = TRUE)) {
    stop("Please install viridisLite: install.packages('viridisLite')")
  }
  
  dt <- as.data.table(copy(ptm_dt))
  um <- as.data.table(copy(unimod_map))
  
  dt[, UnimodID := as.character(UnimodID)]
  um[, UnimodID := as.character(UnimodID)]
  
  # Join names
  dt <- merge(dt, um, by = "UnimodID", all.x = TRUE)
  
  # Build label with Unimod included
  dt[, ModName := fifelse(
    !is.na(UnimodName) & nzchar(UnimodName),
    paste0(UnimodName, " (UniMod:", UnimodID, ")"),
    paste0("Unimod:", UnimodID)
  )]
  
  dt <- dt[is.finite(Count) & Count > 0 & !is.na(UnimodID) & nzchar(UnimodID)]
  
  # Rank by TOTAL using UnimodID (stable)
  top_dt <- dt[, .(Total = sum(Count, na.rm = TRUE), ModName = ModName[1]), by = UnimodID][
    order(-Total)
  ][1:min(top_n, .N)]
  
  top_ids  <- top_dt$UnimodID
  top_mods <- top_dt$ModName
  
  cols <- viridisLite::viridis(length(top_mods), option = "turbo", begin = 0.05, end = 0.95)
  names(cols) <- top_mods
  
  return(list(
    top_mods = top_mods,
    cols = cols,
    top_ids = top_ids
  ))
}




# Separate palette for excluded modifications so they do not reuse the included PTM colors.
get_excluded_mod_palette <- function(ptm_dt, unimod_map, top_n = 50) {
  if (!requireNamespace("viridisLite", quietly = TRUE)) {
    stop("Please install viridisLite: install.packages('viridisLite')")
  }
  pal <- get_top_mod_palette(ptm_dt = ptm_dt, unimod_map = unimod_map, top_n = top_n)
  if (length(pal$top_mods) > 0) {
    cols <- viridisLite::viridis(length(pal$top_mods), option = "magma", begin = 0.15, end = 0.95)
    names(cols) <- pal$top_mods
    pal$cols <- cols
  }
  pal
}
plot_ptm_counts_topn <- function(ptm_dt, unimod_map, group_map, out_png, top_mods, mod_cols) {
  dt <- as.data.table(copy(ptm_dt))
  gm <- as.data.table(copy(group_map))
  gm <- unique(gm[, .(Sample, Group)])
  
  # add Group ONCE
  dt <- merge(dt, gm, by = "Sample", all.x = TRUE)
  dt <- dt[!is.na(Group) & nzchar(Group)]
  
  dt[, UnimodID := as.character(UnimodID)]
  um <- as.data.table(copy(unimod_map))
  um[, UnimodID := as.character(UnimodID)]
  dt <- merge(dt, um, by = "UnimodID", all.x = TRUE)
  
  # exact label format used everywhere
  dt[, ModName := fifelse(
    !is.na(UnimodName) & nzchar(UnimodName),
    paste0(UnimodName, " (UniMod:", UnimodID, ")"),
    paste0("Unimod:", UnimodID)
  )]
  
  dt <- dt[is.finite(Count) & Count > 0]
  dt <- dt[ModName %chin% top_mods]
  
  # SampleKey
  samples_used <- sort(unique(dt$Sample))
  sample_meta <- data.table(Sample = samples_used)
  sample_meta[, SampleKey := vapply(Sample, extract_sample_key, character(1))]
  sample_meta <- merge(sample_meta, gm, by = "Sample", all.x = TRUE)
  
  dt <- merge(dt, sample_meta[, .(Sample, SampleKey)], by = "Sample", all.x = TRUE)
  dt <- dt[!is.na(SampleKey) & nzchar(SampleKey)]
  
  # aggregate
  dt <- dt[, .(Count = sum(Count, na.rm = TRUE)), by = .(Group, SampleKey, ModName)]
  
  key_levels <- make_matched_samplekey_levels(sample_meta)
  dt[, SampleKey := factor(SampleKey, levels = key_levels)]
  grp_levels <- sort(unique(as.character(dt$Group)))
  dt[, Group := factor(Group, levels = grp_levels)]
  dt[, ModName := factor(ModName, levels = top_mods)]
  
  cols <- mod_cols[top_mods]
  
  p <- ggplot(dt, aes(x = SampleKey, y = Count, fill = ModName)) +
    geom_col(width = 0.9, color = "black", linewidth = 0.2) +
    facet_wrap(~ Group, nrow = 1, scales = "free_x") +
    scale_fill_manual(values = cols, drop = FALSE) +
    labs(x = NULL, y = NULL, fill = NULL) +
    pub_theme(base_size = 12, strip_size = 14) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing.x = unit(0.4, "lines"),
      legend.position = "right"
    )
  
  save_pub(p, out_png, width = 15.5, height = 6.0)
  invisible(p)
}

# ------------------------------
# FASTA loading + PTM site mapping
# ------------------------------
# FASTA headers are like "sp|P12345|..." — we want to map accession tokens like P12345
load_fasta_map <- function(fasta_path) {
  aa <- Biostrings::readAAStringSet(fasta_path)
  hdr <- names(aa)
  # primary key = middle token of sp|ACC|NAME or tr|ACC|NAME
  acc <- ifelse(grepl("\\|", hdr), sub("^.*\\|([^|]+)\\|.*$", "\\1", hdr), sub("\\s.*$", "", hdr))
  seqs <- as.character(aa)
  names(seqs) <- acc
  seqs
}

pick_modseq_col <- function(dt) {
  cn <- names(dt)
  cands <- c("Modified.Sequence", "ModifiedPeptideSequence", "Modified.Peptide.Sequence",
             "ModifiedPeptide.Sequence", "ModifiedSequence", "Modified")
  hit <- cands[cands %in% cn][1]
  if (is.na(hit)) stop("Could not find Modified.Sequence-like column in pr_matrix. Columns are: ", paste(cn, collapse = ", "))
  hit
}

pick_protein_col <- function(dt) {
  cn <- names(dt)
  cands <- c("Protein.Ids","Protein.IDs","Protein.Group","ProteinGroup",
             "Protein.Names","Proteins","Protein")
  hit <- cands[cands %in% cn][1]
  if (is.na(hit)) stop("Could not find a protein identifier column in pr_matrix. Columns are: ", paste(cn, collapse = ", "))
  hit
}

# Parse Modified.Sequence where UniMod tag is AFTER the residue: A(UniMod:1)AD...
parse_modified_sequence_to_sites <- function(modseq) {
  if (is.na(modseq) || !nzchar(modseq)) return(list(base = NA_character_, mods = data.table(idx=integer(), UnimodID=character())))
  s <- gsub("\\s+", "", modseq)
  chars <- strsplit(s, "")[[1]]
  
  base_chars <- character()
  mods <- list()
  pep_pos <- 0L
  
  i <- 1L
  while (i <= length(chars)) {
    ch <- chars[i]
    if (grepl("^[A-Z]$", ch)) {
      pep_pos <- pep_pos + 1L
      base_chars <- c(base_chars, ch)
      
      # lookahead for "(UniMod:###)" right after residue
      if (i < length(chars) && chars[i+1] == "(") {
        j <- i+1L
        while (j <= length(chars) && chars[j] != ")") j <- j + 1L
        if (j <= length(chars)) {
          tag <- paste0(chars[(i+2L):(j-1L)], collapse = "")
          m <- regmatches(tag, regexpr("UniMod:\\s*[0-9]+", tag, perl = TRUE))
          if (length(m) == 1 && nzchar(m)) {
            uid <- gsub(".*UniMod:\\s*([0-9]+).*", "\\1", m, perl = TRUE)
            mods[[length(mods)+1]] <- list(idx = pep_pos, UnimodID = uid)
          }
          i <- j
        }
      }
    }
    i <- i + 1L
  }
  
  base <- paste0(base_chars, collapse = "")
  mod_dt <- rbindlist(lapply(mods, as.data.table), fill = TRUE)
  if (nrow(mod_dt) == 0) mod_dt <- data.table(idx=integer(), UnimodID=character())
  list(base = base, mods = mod_dt)
}

# Map peptide (base) to protein positions for each UniMod site
map_peptide_to_protein_positions <- function(protein_id_raw, pep_base, pep_mods, fasta_map) {
  if (is.na(protein_id_raw) || !nzchar(protein_id_raw) || is.na(pep_base) || !nzchar(pep_base)) return(data.table())
  ids <- unlist(strsplit(protein_id_raw, "[;,\t ]+"))
  ids <- ids[nzchar(ids)]
  out <- list()
  for (pid in ids) {
    if (!pid %in% names(fasta_map)) next
    prot <- fasta_map[[pid]]
    hits <- gregexpr(pep_base, prot, fixed = TRUE)[[1]]
    if (length(hits) == 1 && hits[1] == -1) next
    for (start in hits) {
      if (nrow(pep_mods) == 0) next
      site_pos <- start + pep_mods$idx - 1L
      out[[length(out)+1]] <- data.table(
        ProteinID  = pid,
        ProteinPos = site_pos,
        Residue    = substring(prot, site_pos, site_pos),
        UnimodID   = as.character(pep_mods$UnimodID)
      )
    }
  }
  if (length(out) == 0) return(data.table())
  unique(rbindlist(out))
}

# Apply UniMod map with EXACT label format used in stacked bars
apply_unimod_name_map <- function(site_dt, unimod_map) {
  um <- as.data.table(copy(unimod_map))
  um[, UnimodID := as.character(UnimodID)]
  s <- as.data.table(copy(site_dt))
  s[, UnimodID := as.character(UnimodID)]
  s <- merge(s, um, by = "UnimodID", all.x = TRUE)
  s[, ModLabel := fifelse(
    !is.na(UnimodName) & nzchar(UnimodName),
    paste0(UnimodName, " (UniMod:", UnimodID, ")"),
    paste0("Unimod:", UnimodID)
  )]
  s
}

# ------------------------------
# TRIANGLE HEATMAP (upper = intersection tiles, lower = barcodes)
# ------------------------------
plot_triangle_overlap_with_barcodes <- function(
    site_dt, top_mods, top_cols, out_png,
    title = "PTM-site overlap (Intersection) + PTM makeup (Top 20)",
    cell_h = 1.50, cell_w = 1.10,
    diag_tile_alpha = 0.35,
    text_size = 4.6
) {
  if (!requireNamespace("viridisLite", quietly = TRUE)) stop("Install viridisLite")
  if (!requireNamespace("ggnewscale", quietly = TRUE)) stop("Install ggnewscale: install.packages('ggnewscale')")
  
  dt <- data.table::as.data.table(site_dt)
  need <- c("Group","SiteKey","ModLabel")
  miss <- setdiff(need, names(dt))
  if (length(miss) > 0) stop("site_dt missing columns: ", paste(miss, collapse = ", "))
  
  dt[, ModLabel := as.character(ModLabel)]
  dt <- dt[ModLabel %chin% top_mods]
  dt <- unique(dt[, .(Group, SiteKey, ModLabel)])
  
  groups <- sort(unique(dt$Group))
  ng <- length(groups)
  if (ng < 2) stop("Need >=2 groups with sites.")
  
  # sites per group
  sites_by_group <- split(dt$SiteKey, dt$Group)
  
  # --- helper: convert (idx) -> plotting coords with Y flipped
  # idx: 1..ng ; y_plot puts idx=1 at TOP and idx=ng at BOTTOM (like matrix/heatmap)
  idx_to_x <- function(idx) idx * cell_w
  idx_to_y <- function(idx) (ng - idx + 1) * cell_h
  
  # -------------------------
  # Cell grid + intersection N
  # -------------------------
  cell <- data.table::CJ(Group1 = groups, Group2 = groups, unique = TRUE)
  cell[, `:=`(
    i = match(Group1, groups),
    j = match(Group2, groups)
  )]
  
  cell[, IntersectN := {
    si <- unique(sites_by_group[[Group1]])
    sj <- unique(sites_by_group[[Group2]])
    length(intersect(si, sj))
  }, by = .(Group1, Group2)]
  
  cell[, `:=`(
    x = idx_to_x(j),
    y = idx_to_y(i),
    diag = (i == j),
    tri  = fifelse(i <= j, "upper", "lower")   # by index (matrix convention)
  )]
  
  cell[, `:=`(
    xmin = x - 0.5 * cell_w, xmax = x + 0.5 * cell_w,
    ymin = y - 0.5 * cell_h, ymax = y + 0.5 * cell_h
  )]
  
  cell_u <- cell[tri == "upper"]
  cell_u_diag <- cell_u[diag == TRUE]
  cell_u_off  <- cell_u[diag == FALSE]
  
  maxN <- max(cell_u$IntersectN, na.rm = TRUE)
  if (!is.finite(maxN) || maxN <= 0) maxN <- 1
  overlap_cols <- viridisLite::viridis(256)
  
  # -------------------------
  # Barcode composition for LOWER triangle intersections + diagonal (within-group)
  # -------------------------
  comp_list <- list()
  
  # lower triangle (i > j)
  lower_cells <- cell[tri == "lower" & IntersectN > 0]
  if (nrow(lower_cells) > 0) {
    for (k in seq_len(nrow(lower_cells))) {
      g1 <- lower_cells$Group1[k]
      g2 <- lower_cells$Group2[k]
      
      inter_sites <- intersect(unique(sites_by_group[[g1]]), unique(sites_by_group[[g2]]))
      if (length(inter_sites) == 0) next
      
      sub <- unique(dt[Group %chin% c(g1, g2) & SiteKey %chin% inter_sites, .(SiteKey, ModLabel)])
      if (nrow(sub) == 0) next
      
      comp <- sub[, .N, by = ModLabel]
      comp[, Prop := N / sum(N)]
      comp[, `:=`(Group1 = g1, Group2 = g2)]
      comp_list[[length(comp_list)+1]] <- comp
    }
  }
  
  # diagonal (within group)
  for (g in groups) {
    g_sites <- unique(sites_by_group[[g]])
    if (length(g_sites) == 0) next
    
    sub <- unique(dt[Group == g & SiteKey %chin% g_sites, .(SiteKey, ModLabel)])
    if (nrow(sub) == 0) next
    
    comp <- sub[, .N, by = ModLabel]
    comp[, Prop := N / sum(N)]
    comp[, `:=`(Group1 = g, Group2 = g)]
    comp_list[[length(comp_list)+1]] <- comp
  }
  
  comp_all <- data.table::rbindlist(comp_list, fill = TRUE)
  if (nrow(comp_all) == 0) stop("No barcode composition available (unexpected).")
  
  comp_all[, `:=`(
    i = match(Group1, groups),
    j = match(Group2, groups),
    ModLabel = factor(ModLabel, levels = top_mods)
  )]
  
  comp_all[, `:=`(
    x = idx_to_x(j),
    y = idx_to_y(i)
  )]
  
  data.table::setorder(comp_all, Group1, Group2, ModLabel)
  
  # vertical stacking inside each cell (bottom->top within the cell)
  comp_all[, cum := cumsum(Prop), by = .(Group1, Group2)]
  comp_all[, prev := data.table::shift(cum, fill = 0), by = .(Group1, Group2)]
  
  comp_all[, `:=`(
    xmin = x - 0.5 * cell_w,
    xmax = x + 0.5 * cell_w,
    ymin = y - 0.5 * cell_h + prev * cell_h,
    ymax = y - 0.5 * cell_h + cum  * cell_h
  )]
  
  # -------------------------
  # Plot
  # -------------------------
  p <- ggplot2::ggplot() +
    # 1) barcodes (lower + diagonal background)
    ggplot2::geom_rect(
      data = comp_all,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = ModLabel),
      color = NA, linewidth = 0
    ) +
    # Draw only one subtle border per barcode cell instead of a black border around every PTM segment.
    ggplot2::geom_rect(
      data = cell[tri == "lower" | diag == TRUE],
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = NA, color = "grey35", linewidth = 0.25
    ) +
    ggplot2::scale_fill_manual(values = top_cols[top_mods], name = "PTM makeup", drop = FALSE) +
    ggnewscale::new_scale_fill() +
    
    # 2) upper triangle tiles (off-diagonal opaque)
    ggplot2::geom_rect(
      data = cell_u_off,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = IntersectN),
      color = "grey35", linewidth = 0.25
    ) +
    
    # 3) diagonal tile (transparent overlay so barcode shows)
    ggplot2::geom_rect(
      data = cell_u_diag,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = IntersectN),
      color = "grey35", linewidth = 0.25, alpha = diag_tile_alpha
    ) +
    
    ggplot2::scale_fill_gradientn(colours = overlap_cols, limits = c(0, maxN), name = "Intersect\nsites") +
    
    # 4) bold numbers on upper (including diagonal)
    ggplot2::geom_label(
      data = cell_u,
      ggplot2::aes(x = x, y = y, label = as.integer(IntersectN)),
      fontface = "bold",
      size = text_size,
      fill = "white",
      color = "black",
      label.size = 0.15,
      label.padding = grid::unit(0.12, "lines")
    ) +
    
    ggplot2::coord_fixed() +
    ggplot2::scale_x_continuous(
      breaks = idx_to_x(seq_len(ng)),
      labels = groups,
      expand = ggplot2::expansion(add = 0.03)
    ) +
    ggplot2::scale_y_continuous(
      breaks = idx_to_y(seq_len(ng)),
      labels = groups,
      expand = ggplot2::expansion(add = 0.03)
    ) +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    pub_theme(base_size = 14, strip_size = 14) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "right",
      legend.box = "vertical"
    )
  
  w <- max(14, 1.15 * ng * cell_w)
  h <- max(7, 1.05 * ng * cell_h)
  save_pub(p, out_png, width = w, height = h)
  invisible(p)
}

# ------------------------------
# MAIN
# ------------------------------
# fasta_path <- pick_fasta_file()                       # interactive picker (disabled)
# Paths come from environment variables (see the repository README): TWO_PROTEASE_DATA = the DIA-NN
# search results (MassIVE MSV000103457), TWO_PROTEASE_OUT = where reports are written.
# TWO_PROTEASE_FASTA = UniProtKB human reviewed proteome used for site numbering (see README).
fasta_path <- Sys.getenv("TWO_PROTEASE_FASTA", "uniprotkb_Human_AND_reviewed_true_AND_m_2025_10_28.fasta")
# base_dir <- pick_parent_folder()                      # interactive picker (disabled)
base_dir <- Sys.getenv("TWO_PROTEASE_DATA", "data")   # parent folder of the DIA-NN search folders

# ---- Paste/Rscript-safe prompts (2-protease): take all default answers --------
# Every readline() returns "" so each prompt uses its DEFAULT branch and can never
# consume a following pasted line. Set ASSUME_DEFAULTS <- FALSE to answer prompts
# yourself when running via source() in an interactive console.
ASSUME_DEFAULTS <- TRUE
if (ASSUME_DEFAULTS) readline <- function(prompt = "") { if (nzchar(prompt)) cat(prompt); "" }

timestamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
out_dir   <- file.path(Sys.getenv("TWO_PROTEASE_OUT", "output"), "03_ptm_coverage",
                       paste0("DIANN_PTM_Report_", timestamp))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "run_log.txt")
sink(log_file, split = TRUE)
cat("DIA-NN PTM report started: ", format(Sys.time()), "\n")
cat("Base dir: ", base_dir, "\n")
cat("Output dir: ", out_dir, "\n\n")

on.exit({
  cat("\nRun finished: ", format(Sys.time()), "\n")
  sink()
}, add = TRUE)

# Where the UniMod library is stored (persistent per base_dir)
unimod_lib_tsv <- file.path(base_dir, "unimod_library.tsv")

tryCatch({
  
  message("Discovering search folders under: ", base_dir)
  search_paths <- discover_search_folders(base_dir, pg_file = pg_file, pr_file = pr_file)
  search_paths <- make_unique_names(search_paths)
  
  if (length(search_paths) < 2) stop("Need at least 2 valid search folders containing ", pg_file, " and ", pr_file, ".")
  if (is.finite(max_searches) && length(search_paths) > max_searches) {
    message("Found ", length(search_paths), " searches; using first ", max_searches, " (alphabetical).")
    search_paths <- search_paths[seq_len(max_searches)]
  }
  
  search_paths <- prompt_search_short_names(search_paths)
  
  fwrite(
    data.table(FolderName = basename(search_paths), ShortName = names(search_paths), Path = as.character(search_paths)),
    file.path(out_dir, "search_name_mapping.csv")
  )
  
  message("\nLoading pr_matrix.tsv + collapsing to peptides by Modified.Sequence ...")
  peptide_mats <- list()
  
  # also build mapping ModifiedSequence -> ProteinIDRaw from ALL searches (only for modified peptides)
  map_list <- list()
  
  for (s in names(search_paths)) {
    sp <- search_paths[[s]]
    pr_path <- file.path(sp, pr_file)
    stop_if_missing(pr_path)
    
    pr <- fread(pr_path, sep = "\t", header = TRUE, data.table = TRUE, showProgress = FALSE)
    
    # collapse to peptides (rows = Modified.Sequence)
    pep_m <- collapse_pr_to_peptides_by_modified_sequence(
      pr_dt = pr,
      collapse_method = peptide_collapse_method,
      zero_as_na = zero_as_na,
      id_col = "Modified.Sequence"
    )
    peptide_mats[[s]] <- pep_m
    message("  ", s, " | peptides: ", nrow(pep_m), " | samples: ", ncol(pep_m))
    
    # build modseq->protein mapping rows from this search (modified peptides only)
    mod_col <- pick_modseq_col(pr)
    prot_col <- pick_protein_col(pr)
    
    md <- pr[, .(
      ModifiedSequence = sanitize_ids(get(mod_col)),
      ProteinIDRaw     = as.character(get(prot_col))
    )]
    md <- md[!is.na(ModifiedSequence) & nzchar(ModifiedSequence)]
    md <- md[grepl("UniMod:", ModifiedSequence, fixed = TRUE)]
    md <- md[!is.na(ProteinIDRaw) & nzchar(ProteinIDRaw)]
    md <- unique(md, by = "ModifiedSequence")
    map_list[[s]] <- md
  }
  
  # Align sample columns across searches
  pep_aligned <- align_samples_union(peptide_mats)
  peptide_mats <- pep_aligned$mats
  common_samples <- pep_aligned$samples
  fwrite(data.table(Sample = common_samples), file.path(out_dir, "sample_list_common.csv"))
  message("\nUnion sample columns: ", length(common_samples))
  
  # Group mapping
  group_map <- build_group_map_auto(common_samples)
  # Hardcoded manuscript labels (Trypsin 1/2, Krakatoa 1/2, Mix 20min/30min).
  group_map <- apply_manuscript_group_labels(group_map)
  fwrite(group_map, file.path(out_dir, "sample_group_mapping.csv"))
  cat("\nGroup counts:\n")
  print(as.data.table(group_map)[, .N, by = Group][order(-N)])
  
  # Build UNION peptide matrix
  message("\nBuilding UNION peptide matrix...")
  m_union_pep <- build_union_matrix(peptide_mats)
  
  # Modified peptides only for site mapping (we still use all for union presence matrix, but sites only need modified)
  all_peps <- rownames(m_union_pep)
  mod_peps <- all_peps[grepl("UniMod:", all_peps, fixed = TRUE)]
  message("Modified peptides only: ", length(mod_peps), " / ", length(all_peps))
  
  # PTM counts (presence-based)
  message("\nCounting PTMs from UNION peptide presence...")
  ptm_dt <- build_ptm_counts_from_union_peptides(m_union_pep, exclude_ids = exclude_unimod_ids)
  
  fwrite(ptm_dt, file.path(out_dir, "PTM_counts_per_sample_raw.csv"))
  
  # UniMod names (fast cached library)
  unimod_ids_seen <- unique(ptm_dt$UnimodID)
  unimod_ids_seen <- unimod_ids_seen[nzchar(unimod_ids_seen)]
  unimod_map <- get_unimod_name_map_fast(unimod_ids_seen, lib_tsv = unimod_lib_tsv)
  fwrite(unimod_map, file.path(out_dir, "PTM_unimod_id_to_name.csv"))
  
  # Top20 palette = single source of truth
  pal <- get_top_mod_palette(ptm_dt, unimod_map, top_n = TOP_N_PTMS)
  top20_mods <- pal$top_mods
  top20_cols <- pal$cols
  top20_ids  <- pal$top_ids
  
  
  # PTM stacked bars
  plot_path <- file.path(out_dir, "PTM_Counts_Top20_byGroup.png")
  plot_ptm_counts_topn(
    ptm_dt = ptm_dt,
    unimod_map = unimod_map,
    group_map = group_map,
    out_png = plot_path,
    top_mods = top20_mods,
    mod_cols = top20_cols
  )
  write_stackedbar_order_csv(
    group_map = group_map,
    out_samples_csv = file.path(out_dir, "StackedBar_Sample_Order.csv"),
    out_keys_csv    = file.path(out_dir, "StackedBar_SampleKey_Order.csv")
  )
  message("Wrote: StackedBar_Sample_Order.csv + StackedBar_SampleKey_Order.csv")

  # ------------------------------------------------------------
  # EXCLUDED MODIFICATIONS: counts + stacked bar plot
  # ------------------------------------------------------------
  message("\nCounting EXCLUDED PTMs from UNION peptide presence...")
  ptm_excluded_dt <- build_ptm_counts_for_unimod_ids(
    m_union_pep = m_union_pep,
    keep_ids = exclude_unimod_ids
  )
  fwrite(ptm_excluded_dt, file.path(out_dir, "Excluded_PTM_counts_per_sample_raw.csv"))

  excluded_ids_seen <- sort(unique(as.character(ptm_excluded_dt$UnimodID)))
  excluded_ids_seen <- excluded_ids_seen[nzchar(excluded_ids_seen)]
  message("Excluded UniMod IDs found in peptide matrix: ",
          ifelse(length(excluded_ids_seen) == 0, "none", paste(excluded_ids_seen, collapse = ", ")))
  message("Excluded PTM count rows: ", nrow(ptm_excluded_dt))

  excluded_unimod_map <- data.table(UnimodID = character(), UnimodName = character())
  excluded_mods <- character()
  excluded_cols <- character()
  excluded_ids  <- character()

  if (nrow(ptm_excluded_dt) > 0 && length(excluded_ids_seen) > 0) {
    excluded_unimod_map <- get_unimod_name_map_fast(excluded_ids_seen, lib_tsv = unimod_lib_tsv)
    fwrite(excluded_unimod_map, file.path(out_dir, "Excluded_PTM_unimod_id_to_name.csv"))

    excluded_pal <- get_excluded_mod_palette(
      ptm_dt = ptm_excluded_dt,
      unimod_map = excluded_unimod_map,
      top_n = length(excluded_ids_seen)
    )
    excluded_mods <- excluded_pal$top_mods
    excluded_cols <- excluded_pal$cols
    excluded_ids  <- excluded_pal$top_ids

    if (length(excluded_mods) > 0) {
      excluded_plot_path <- file.path(out_dir, "Excluded_PTM_Counts_byGroup.png")
      plot_ptm_counts_topn(
        ptm_dt = ptm_excluded_dt,
        unimod_map = excluded_unimod_map,
        group_map = group_map,
        out_png = excluded_plot_path,
        top_mods = excluded_mods,
        mod_cols = excluded_cols
      )
      message("Wrote: ", basename(excluded_plot_path))
    }
  } else {
    fwrite(data.table(UnimodID = character(), UnimodName = character()),
           file.path(out_dir, "Excluded_PTM_unimod_id_to_name.csv"))
    message("No excluded PTMs were detected in the union peptide matrix, so excluded PTM count plot was skipped.")
  }
  
  
  # FASTA + site mapping
  message("\nLoading FASTA...")
  fasta_map <- load_fasta_map(fasta_path)
  message("FASTA proteins loaded: ", length(fasta_map))
  
  message("\nBuilding PTM-site table (Protein position mapping)...")
  
  # union mapping table across searches
  map_dt <- unique(rbindlist(map_list, fill = TRUE), by = "ModifiedSequence")
  setkey(map_dt, ModifiedSequence)
  
  # keep only modified peptides we will map
  mod_peps <- mod_peps[mod_peps %chin% map_dt$ModifiedSequence]
  
  message("Modified peptides: ", length(mod_peps), " | mapped to ProteinIDRaw: ", sum(mod_peps %chin% map_dt$ModifiedSequence))
  
  # parse modified peptides (once)
  parsed <- lapply(mod_peps, parse_modified_sequence_to_sites)
  pep_base <- vapply(parsed, `[[`, character(1), "base")
  pep_mods <- lapply(parsed, `[[`, "mods")
  
  # build peptide->protein site list (once per peptide)
  site_pep_list <- list()
  excluded_site_pep_list <- list()
  for (i in seq_along(mod_peps)) {
    ms <- mod_peps[i]
    pb <- pep_base[i]
    pm <- pep_mods[[i]]
    if (!nzchar(pb) || nrow(pm) == 0) next
    
    pid_raw <- map_dt[J(ms), ProteinIDRaw]
    if (length(pid_raw) == 0 || is.na(pid_raw) || !nzchar(pid_raw)) next
    
    # Split included vs excluded UniMod IDs at the site level.
    # Main analysis keeps the original behavior: excluded IDs are removed.
    pm_included <- pm[!(pm$UnimodID %chin% as.character(exclude_unimod_ids))]
    pm_excluded <- pm[ pm$UnimodID %chin% as.character(exclude_unimod_ids)]

    if (nrow(pm_included) > 0) {
      sites <- map_peptide_to_protein_positions(pid_raw, pb, pm_included, fasta_map)
      if (nrow(sites) > 0) {
        sites[, ModifiedSequence := ms]
        site_pep_list[[length(site_pep_list)+1]] <- sites
      }
    }

    if (nrow(pm_excluded) > 0) {
      excluded_sites <- map_peptide_to_protein_positions(pid_raw, pb, pm_excluded, fasta_map)
      if (nrow(excluded_sites) > 0) {
        excluded_sites[, ModifiedSequence := ms]
        excluded_site_pep_list[[length(excluded_site_pep_list)+1]] <- excluded_sites
      }
    }
  }
  site_pep_dt <- rbindlist(site_pep_list, fill = TRUE)
  excluded_site_pep_dt <- rbindlist(excluded_site_pep_list, fill = TRUE)
  message("excluded_site_pep_dt rows: ", nrow(excluded_site_pep_dt))
  message("site_pep_dt rows: ", nrow(site_pep_dt))
  if (nrow(site_pep_dt) == 0) stop("PTM-site mapping produced 0 rows. Check FASTA IDs vs Protein IDs and peptide matching.")
  
  # expand per sample by presence (modified peptides only)
  # expand per sample by presence (modified peptides only)
  message("\nExpanding PTM sites per sample by peptide presence...")
  
  # Presence table WITH Search (stacked from per-search peptide matrices)
  pres_dt <- build_presence_dt_with_search(peptide_mats)
  fwrite(pres_dt, file.path(out_dir, "PTM_peptide_presence_with_search.csv"))
  
  # keep only modified peptides we mapped
  pres_dt <- pres_dt[Modified.Sequence %chin% mod_peps]
  
  message("pres_dt rows: ", nrow(pres_dt))
  
  # --- normalize join key names
  # site_pep_dt currently has ModifiedSequence; make it Modified.Sequence to match pres_dt
  if ("ModifiedSequence" %in% names(site_pep_dt) && !"Modified.Sequence" %in% names(site_pep_dt)) {
    setnames(site_pep_dt, "ModifiedSequence", "Modified.Sequence")
  }
  
  # Join: site-per-peptide sites x presence (Sample+Search)
  setkey(site_pep_dt, "Modified.Sequence")
  setkey(pres_dt, "Modified.Sequence")
  
  site_dt <- site_pep_dt[pres_dt, nomatch = 0, allow.cartesian = TRUE]
  
  # attach Group ONCE
  gm <- as.data.table(group_map)[, .(Sample, Group)]
  site_dt <- merge(site_dt, gm, by = "Sample", all.x = TRUE)
  site_dt <- site_dt[!is.na(Group) & nzchar(Group)]
  
  message("joined site_dt rows: ", nrow(site_dt))
  
  
  # UniMod names for site_dt (re-use cached library; only fetch unknown IDs)
  site_unimods <- sort(unique(as.character(site_dt$UnimodID)))
  site_unimods <- site_unimods[nzchar(site_unimods)]
  unimod_map2 <- get_unimod_name_map_fast(site_unimods, lib_tsv = unimod_lib_tsv)
  site_dt <- apply_unimod_name_map(site_dt, unimod_map2)
  
  # SiteKey (used for overlap)
  site_dt[, SiteKey := paste(ProteinID, ProteinPos, ModLabel, sep = "|")]
  
  fwrite(site_dt, file.path(out_dir, "PTM_sites_per_sample.csv"))
  message("Site table written: PTM_sites_per_sample.csv (rows: ", nrow(site_dt), ")")

  # ------------------------------------------------------------
  # EXCLUDED MODIFICATIONS: site table + overlap plots
  # ------------------------------------------------------------
  message("\nExpanding EXCLUDED PTM sites per sample by peptide presence...")
  excluded_site_dt <- data.table()

  if (nrow(excluded_site_pep_dt) > 0) {
    if ("ModifiedSequence" %in% names(excluded_site_pep_dt) && !"Modified.Sequence" %in% names(excluded_site_pep_dt)) {
      setnames(excluded_site_pep_dt, "ModifiedSequence", "Modified.Sequence")
    }

    setkey(excluded_site_pep_dt, "Modified.Sequence")
    setkey(pres_dt, "Modified.Sequence")
    excluded_site_dt <- excluded_site_pep_dt[pres_dt, nomatch = 0, allow.cartesian = TRUE]
    excluded_site_dt <- merge(excluded_site_dt, gm, by = "Sample", all.x = TRUE)
    excluded_site_dt <- excluded_site_dt[!is.na(Group) & nzchar(Group)]

    if (nrow(excluded_site_dt) > 0) {
      excluded_site_unimods <- sort(unique(as.character(excluded_site_dt$UnimodID)))
      excluded_site_unimods <- excluded_site_unimods[nzchar(excluded_site_unimods)]

      # Ensure the excluded label map/palette includes every excluded site ID.
      missing_excluded_ids <- setdiff(excluded_site_unimods, excluded_unimod_map$UnimodID)
      if (length(missing_excluded_ids) > 0) {
        excluded_unimod_map <- unique(rbind(
          excluded_unimod_map,
          get_unimod_name_map_fast(missing_excluded_ids, lib_tsv = unimod_lib_tsv),
          fill = TRUE
        ), by = "UnimodID")
        fwrite(excluded_unimod_map, file.path(out_dir, "Excluded_PTM_unimod_id_to_name.csv"))
      }

      excluded_site_dt <- apply_unimod_name_map(excluded_site_dt, excluded_unimod_map)
      excluded_site_dt[, SiteKey := paste(ProteinID, ProteinPos, ModLabel, sep = "|")]
      fwrite(excluded_site_dt, file.path(out_dir, "Excluded_PTM_sites_per_sample.csv"))
      message("Excluded site table written: Excluded_PTM_sites_per_sample.csv (rows: ", nrow(excluded_site_dt), ")")

      # If no count-based palette was made, make a site-count palette as fallback.
      if (length(excluded_mods) == 0) {
        excluded_site_counts <- excluded_site_dt[, .(Sample = Sample[1], UnimodID = as.character(UnimodID), Count = .N), by = UnimodID]
        excluded_pal <- get_excluded_mod_palette(
          ptm_dt = excluded_site_counts,
          unimod_map = excluded_unimod_map,
          top_n = length(unique(excluded_site_counts$UnimodID))
        )
        excluded_mods <- excluded_pal$top_mods
        excluded_cols <- excluded_pal$cols
        excluded_ids  <- excluded_pal$top_ids
      }

      excluded_groups_with_sites <- sort(unique(excluded_site_dt$Group))
      message("Excluded PTM groups with mapped sites: ", paste(excluded_groups_with_sites, collapse = ", "))

      if (length(excluded_groups_with_sites) >= 2 && length(excluded_mods) > 0) {
        excluded_tri_png <- file.path(out_dir, "Excluded_PTM_Site_Overlap_Triangle.png")
        plot_triangle_overlap_with_barcodes(
          site_dt = excluded_site_dt,
          top_mods = excluded_mods,
          top_cols = excluded_cols,
          out_png = excluded_tri_png,
          title = "Excluded PTM-site overlap + excluded PTM makeup"
        )
        message("Wrote: ", basename(excluded_tri_png))
      } else {
        message("Skipped excluded overlap triangle: need at least 2 groups with mapped excluded sites and at least 1 excluded modification label.")
      }

      excluded_site_sets <- build_strategy_site_sets_by_unimodid(excluded_site_dt, unimod_keep = excluded_ids)
      message("Excluded strategy site-set sizes: ",
              ifelse(length(excluded_site_sets) == 0, "none", paste(names(excluded_site_sets), sapply(excluded_site_sets, length), sep="=", collapse=" | ")))

      if (length(excluded_site_sets) == 2 && length(excluded_mods) > 0) {
        excluded_venn_png <- file.path(out_dir, "Excluded_PTM_Site_Overlap_ByStrategy_Venn.png")
        plot_venn2_ptmfilled_proportional(
          site_dt   = excluded_site_dt,
          site_sets = excluded_site_sets,
          mod_keep  = excluded_mods,
          mod_cols  = excluded_cols,
          out_png   = excluded_venn_png,
          title     = "Excluded PTM-site overlap by search strategy"
        )
        message("Wrote: ", basename(excluded_venn_png))
      } else {
        message("Skipped excluded Venn: proportional Venn currently requires exactly 2 search strategies with excluded mapped sites.")
      }
    } else {
      fwrite(data.table(), file.path(out_dir, "Excluded_PTM_sites_per_sample.csv"))
      message("Skipped excluded site plots: excluded site mapping joined to 0 peptide/sample/search presence rows.")
    }
  } else {
    fwrite(data.table(), file.path(out_dir, "Excluded_PTM_sites_per_sample.csv"))
    message("Skipped excluded site plots: no excluded modification sites mapped to FASTA positions.")
  }

  message("\nBuilding Euler/Venn overlap by search strategy (Top20 mods only)...")
  site_sets <- build_strategy_site_sets_by_unimodid(site_dt, unimod_keep = top20_ids)
  
  message("Strategy site-set sizes: ",
          paste(names(site_sets), sapply(site_sets, length), sep="=", collapse=" | "))

  venn_png <- file.path(out_dir, "PTM_Site_Overlap_ByStrategy_Venn_Top20.png")
  
  plot_venn2_ptmfilled_proportional(
    site_dt   = site_dt,
    site_sets = site_sets,
    mod_keep  = top20_mods,   # SAME 20 labels as stacked bar
    mod_cols  = top20_cols,   # SAME colors as stacked bar
    out_png   = venn_png,
    title     = "PTM-site overlap by search strategy (Top20 PTMs)"
  )
  message("Wrote: ", basename(venn_png))

  
  # Triangle heatmap: Top20 mods ONLY, exact same colors as stacked bars, no Other
  tri_png <- file.path(out_dir, "PTM_Site_Overlap_Triangle_Top20.png")
  plot_triangle_overlap_with_barcodes(
    site_dt = site_dt,
    top_mods = top20_mods,
    top_cols = top20_cols,
    out_png = tri_png,
    title = "PTM-site overlap (Intersection) + PTM makeup (Top 20)"
  )
  message("Wrote: ", basename(tri_png))
  
  # simple index
  index_path <- file.path(out_dir, "index.html")
  writeLines(c(
    "<html><head><meta charset='utf-8'><title>DIA-NN PTM Report</title></head><body>",
    paste0("<h1>DIA-NN PTM Report</h1><p>Generated: ", format(Sys.time()), "</p>"),
    "<h2>Outputs</h2><ul>",
    "<li><a href='search_name_mapping.csv'>search_name_mapping.csv</a></li>",
    "<li><a href='sample_list_common.csv'>sample_list_common.csv</a></li>",
    "<li><a href='sample_group_mapping.csv'>sample_group_mapping.csv</a></li>",
    "<li><a href='PTM_counts_per_sample_raw.csv'>PTM_counts_per_sample_raw.csv</a></li>",
    "<li><a href='Excluded_PTM_counts_per_sample_raw.csv'>Excluded_PTM_counts_per_sample_raw.csv</a></li>",
    "<li><a href='PTM_unimod_id_to_name.csv'>PTM_unimod_id_to_name.csv</a></li>",
    "<li><a href='Excluded_PTM_unimod_id_to_name.csv'>Excluded_PTM_unimod_id_to_name.csv</a></li>",
    "<li><a href='PTM_Counts_Top20_byGroup.png'>PTM_Counts_Top20_byGroup.png</a></li>",
    "<li><a href='Excluded_PTM_Counts_byGroup.png'>Excluded_PTM_Counts_byGroup.png</a></li>",
    "<li><a href='PTM_sites_per_sample.csv'>PTM_sites_per_sample.csv</a></li>",
    "<li><a href='Excluded_PTM_sites_per_sample.csv'>Excluded_PTM_sites_per_sample.csv</a></li>",
    "<li><a href='PTM_Site_Overlap_Triangle_Top20.png'>PTM_Site_Overlap_Triangle_Top20.png</a></li>",
    "<li><a href='Excluded_PTM_Site_Overlap_Triangle.png'>Excluded_PTM_Site_Overlap_Triangle.png</a></li>",
    "<li><a href='StackedBar_SampleKey_Order.csv'>StackedBar_SampleKey_Order.csv</a></li>",
    "<li><a href='PTM_Site_Overlap_ByStrategy_Venn_Top20.png'>PTM_Site_Overlap_ByStrategy_Venn_Top20.png</a></li>",
    "<li><a href='Excluded_PTM_Site_Overlap_ByStrategy_Venn.png'>Excluded_PTM_Site_Overlap_ByStrategy_Venn.png</a></li>",
    "</ul></body></html>"
  ), con = index_path)
  
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

