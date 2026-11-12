# Shared config + helpers for the 07_verification analyses (Mix30 vs Trypsin-2 completeness).
# Sourced by completeness_mix30_vs_try2.R, discovery_lens.R and cross_protease_recovery.R after
# ../config.R, which defines SEARCH (= TWO_PROTEASE_DATA).

PG_TRY <- file.path(SEARCH, "Trypsin_Library_Search",  "trypsin_library_search_report.pg_matrix.tsv")
PG_KRA <- file.path(SEARCH, "Krakatoa_Library_Search", "report.pg_matrix.tsv")
RP_TRY <- file.path(SEARCH, "Reproducibility_Search", "Trypsin Library",  "report.pg_matrix.tsv")
RP_HTA <- file.path(SEARCH, "Reproducibility_Search", "Krakatoa Library", "report.pg_matrix.tsv")

ANN <- c("Genes","Protein.Group","Protein.Ids","Protein.Names",
         "First.Protein.Description","N.Sequences","N.Proteotypic.Sequences")
MIN_PAIRED <- 5L          # co-observed runs needed for a protein's offset (strict gate)
N_BOOT     <- 2000L
TAUS_REP   <- c(0.50, 0.75, 0.90, 1.00)
DPI        <- 600

COL <- c("Try2" = "#2C6FB3", "Try2-union" = "#8FB3D9", "Mix T-ch" = "#7A7A7A",
         "Mix deliv" = "#E08214", "Mix +imp" = "#1B7837", "Mix ceiling" = "#1B7837",
         "Try1" = "#B0B0B0")
COMP_LEVELS <- c("0-20%","20-50%","50-80%","80-95%","95-100%","100%")
COMP_COLS <- setNames(colorRampPalette(c("red2","gold","green3"))(6), COMP_LEVELS)   # palette from 02_identification/Identification_Analysis.R
bin_c <- function(pct) fcase(pct >= 100, "100%", pct >= 95, "95-100%", pct >= 80, "80-95%",
                             pct >= 50, "50-80%", pct >= 20, "20-50%", pct > 0, "0-20%", default = NA_character_)

# ------------------------------------------------------------------ loaders --
melt_pg <- function(dt) {
  rc <- setdiff(names(dt), ANN)
  dt[, (rc) := lapply(.SD, function(x) suppressWarnings(as.numeric(x))), .SDcols = rc]
  L <- melt(dt, id.vars = "Protein.Group", measure.vars = rc,
            variable.name = "Run", value.name = "Intensity", variable.factor = FALSE)
  L[!is.finite(Intensity) | Intensity <= 0, Intensity := NA_real_]
  L <- L[, .(Intensity = if (all(is.na(Intensity))) NA_real_ else sum(Intensity, na.rm = TRUE)),
         by = .(FeatureID = Protein.Group, Run)]
  L[!is.na(FeatureID) & FeatureID != ""]
}
load_disc <- function(path) {
  L <- melt_pg(fread(path, sep = "\t", header = TRUE, showProgress = FALSE))
  b <- sub("\\.raw$", "", basename(gsub("\\\\", "/", L$Run)), ignore.case = TRUE)
  L[, Cond := fcase(grepl("Plasma_Mix", b) &  grepl("20min", b), "Mix20",
                    grepl("Plasma_Mix", b) & !grepl("20min", b), "Mix30",
                    grepl("^251117_DIA_Plasma_Try", b), "Try2",
                    grepl("KreimerSimion_EC_BOB", b), "Try1",
                    grepl("^251117_DIA_Plasma_Kra", b), "Kra2",
                    grepl("^250807_DIA_Plasma_Kra", b), "Kra1", default = "OTHER")]
  # sample ID: BRH/HMN code for the paired sets; BOB run number for Try1 (unpaired)
  L[, SampleID := fifelse(Cond == "Try1", sub(".*_BOB_Plasma_", "BOB", b),
                          sub("_20min$", "", sub(".*Plasma_(Mix|Try|Kra)_", "", b)))][]
}

# ------------------------------------------------ Panel A guardrail (repro) --
assert_panel_a <- function() {
  cnt <- function(path) {
    L <- melt_pg(fread(path, sep = "\t", header = TRUE, showProgress = FALSE))
    b <- sub("\\.raw$", "", basename(gsub("\\\\", "/", L$Run)), ignore.case = TRUE)
    p <- tstrsplit(b, "_", fixed = TRUE); n <- length(p)
    L[, Day := as.integer(p[[n-2]])][, Group := p[[n-1]]]
    dm <- L[is.finite(Intensity), .(DayMean = mean(Intensity)), by = .(FeatureID, Group, Day)]
    ed <- L[, .(ED = uniqueN(Day)), by = Group]
    cv <- dm[, .(CV = if (.N < 2) NA_real_ else 100*sd(DayMean)/mean(DayMean), MD = .N), by = .(FeatureID, Group)]
    cv <- merge(cv, ed, by = "Group")
    so <- as.character(sort(as.integer(unique(L$Group))))
    w <- cv[MD == ED & MD >= 2 & is.finite(CV) & CV <= 20, .N, by = Group]
    v <- w$N[match(so, w$Group)]; v[is.na(v)] <- 0L; as.integer(v)
  }
  ok <- identical(cnt(RP_TRY), c(120L,107L,127L,123L,123L,123L,146L,131L)) &&
        identical(cnt(RP_HTA), c(92L,80L,88L,89L,108L,101L,91L,82L))
  cat("Panel A guardrail:", if (ok) "PASSED" else "FAILED", "\n"); if (!ok) stop("Panel A failed")
}

# ------------------------------------------------ presence-matrix builders --
# wide logical matrix: proteins (rows) x samples (cols, in `samples` order)
pres_mat <- function(L, samples) {
  w <- dcast(L[is.finite(Intensity), .(FeatureID, SampleID, p = TRUE)], FeatureID ~ SampleID,
             value.var = "p", fill = FALSE)
  for (s in setdiff(samples, names(w))) w[, (s) := FALSE]
  M <- as.matrix(w[, ..samples]); rownames(M) <- w$FeatureID; M
}
int_mat <- function(L, samples) {
  w <- dcast(L[is.finite(Intensity), .(FeatureID, SampleID, Intensity)], FeatureID ~ SampleID,
             value.var = "Intensity")
  for (s in setdiff(samples, names(w))) w[, (s) := NA_real_]
  M <- as.matrix(w[, ..samples]); rownames(M) <- w$FeatureID; M
}
align <- function(M, rows) {        # expand to a row set, FALSE/NA for absent proteins
  out <- matrix(if (is.logical(M)) FALSE else NA_real_, length(rows), ncol(M),
                dimnames = list(rows, colnames(M)))
  k <- intersect(rows, rownames(M)); out[k, ] <- M[k, ]; out
}

# strict cross-protease fill of `tgt` from `don` (log2 intensity matrices, same dims):
# per-protein median offset from >= MIN_PAIRED co-observed runs; target must have >= 1 real obs
strict_fill <- function(tgt, don) {
  lt <- log2(tgt); ld <- log2(don)
  co <- is.finite(lt) & is.finite(ld)
  npair <- rowSums(co)
  d <- lt - ld; d[!co] <- NA
  off <- apply(d, 1, median, na.rm = TRUE)
  ok  <- npair >= MIN_PAIRED & rowSums(is.finite(lt)) > 0
  fill <- !is.finite(lt) & is.finite(ld) & ok
  list(fill = fill, off = off, npair = npair)
}

# Build all Mix arms for one mixture condition
build_mix_arms <- function(TR, KR, cond, samples) {
  Ti <- int_mat(TR[Cond == cond], samples); Ki <- int_mat(KR[Cond == cond], samples)
  prots <- union(rownames(Ti), rownames(Ki))
  Ti <- align(Ti, prots); Ki <- align(Ki, prots)
  Tp <- is.finite(Ti); Kp <- is.finite(Ki)
  fT <- strict_fill(Ti, Ki); fK <- strict_fill(Ki, Ti)
  chosen <- ifelse(rowSums(Tp) >= rowSums(Kp), "Trypsin", "Krakatoa")   # ties -> trypsin
  C <- Tp; C[chosen == "Krakatoa", ] <- Kp[chosen == "Krakatoa", ]
  Fill <- fT$fill; Fill[chosen == "Krakatoa", ] <- fK$fill[chosen == "Krakatoa", ]
  D <- C | Fill
  list(B = Tp[rowSums(Tp) > 0, , drop = FALSE], C = C, D = D, E = Tp | Kp,
       chosen = chosen, n_filled = sum(Fill), Tp = Tp, Kp = Kp,
       sep9_union_imp = (Tp | fT$fill), sep9_union_impK = (Kp | fK$fill))
}

# counts of proteins with completeness >= k/ncol for k = 1..ncol
ge_counts <- function(M, idx = seq_len(ncol(M))) {
  n <- rowSums(M[, idx, drop = FALSE]); n <- n[n > 0]
  rev(cumsum(rev(tabulate(n, nbins = length(idx)))))
}

