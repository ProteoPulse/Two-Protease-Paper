# =============================================================================
# cv_common.R — shared data + counting code for the Fig 4A CV figures.
#
# Sourced by make_manuscript_figures.R and make_panel_figures.R so both draw
# from the same numbers. Loads the reproducibility-cohort pg_matrix files,
# computes per protein x sample inter-day CVs under two definitions, applies
# cross-protease gap filling, asserts the published Panel A counts, and writes
# tables/manuscript_figure_counts.csv.
#
#   published/  CV = 100*sd(3 day-means)/mean      filter = all 3 days detected
#               (the definition used for the manuscript's Fig 4A)
#   corrected/  CV = injection noise removed:      filter = all 24 injections
#               sigma2_day = var(day means) - MS_within/8, floored at 0
#
# Imputation: cross-protease gap fill into each protease's OWN matrix (never
# averaged); donor rescaled by the protein's median log2 offset over >= 5
# co-observed runs; measured values never altered.
# =============================================================================


suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

# Expects config.R to have been sourced (REPRO, step_dirs()).
OUT  <- step_dirs("05_cross_protease")$root
TBL  <- file.path(OUT, "tables")
for (d in c(file.path(OUT,"figures","corrected"), file.path(OUT,"figures","published"), TBL))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

CV_THRESHOLD <- 20; N_REP <- 8L; MIN_PAIRED <- 5L; PUB_DPI <- 600
COL   <- c(Krakatoa = "#E1812C", Trypsin = "#3C7DA6", Combined = "#7D74B8")
LIGHT <- c(Krakatoa = "#F4CDA4", Trypsin = "#AECBDE", Combined = "#CFC9E6")
PROTEASES <- c("Krakatoa", "Combined", "Trypsin")   # Panel A row order
TIERS <- c(strict = "every day keeps >=1 real measurement",
           onefill = "at most 1 day wholly imputed",
           partial = "at most 2 days wholly imputed",
           ungated = "no restriction (upper bound)")

# ------------------------------------------------------------------ data ----
load_pg <- function(path) {
  dt <- fread(path, sep = "\t", header = TRUE, showProgress = FALSE)
  ann <- c("Genes","Protein.Group","Protein.Ids","Protein.Names",
           "First.Protein.Description","N.Sequences","N.Proteotypic.Sequences")
  L <- melt(dt, id.vars = "Protein.Group", measure.vars = setdiff(names(dt), ann),
            variable.name = "Sample", value.name = "Intensity", variable.factor = FALSE)
  L[, Intensity := as.numeric(Intensity)]
  L[!is.finite(Intensity) | Intensity <= 0, Intensity := NA_real_]
  L <- L[, .(Intensity = if (all(is.na(Intensity))) NA_real_ else sum(Intensity, na.rm = TRUE)),
         by = .(FeatureID = Protein.Group, Sample)]
  L <- L[!is.na(FeatureID) & FeatureID != ""]
  b <- sub("\\.raw$", "", basename(gsub("\\\\", "/", L$Sample)), ignore.case = TRUE)
  p <- tstrsplit(b, "_", fixed = TRUE); n <- length(p)
  L[, `:=`(Day = as.integer(p[[n-2]]), SampleID = p[[n-1]], Replicate = as.integer(p[[n]]),
           imputed = FALSE)]
  L[is.finite(Intensity)]
}
make_offsets <- function(target, donor) {
  m <- merge(target[, .(FeatureID, Sample, tI = Intensity)],
             donor[,  .(FeatureID, Sample, dI = Intensity)], by = c("FeatureID","Sample"))
  m[, d := log2(tI) - log2(dI)]
  m[, .(offset = median(d), n_paired = .N), by = FeatureID][n_paired >= MIN_PAIRED]
}
impute <- function(target, donor, off, runs) {
  grid <- CJ(FeatureID = unique(target$FeatureID), Sample = runs$Sample, unique = TRUE)
  g <- merge(grid, target[, .(FeatureID, Sample, Intensity)], by = c("FeatureID","Sample"), all.x = TRUE)
  g <- merge(g, donor[, .(FeatureID, Sample, dI = Intensity)], by = c("FeatureID","Sample"), all.x = TRUE)
  g <- merge(g, off[, .(FeatureID, offset)], by = "FeatureID", all.x = TRUE)
  g[, imputed := is.na(Intensity) & is.finite(dI) & !is.na(offset)]
  g[imputed == TRUE, Intensity := 2^(log2(dI) + offset)]
  merge(g[is.finite(Intensity), .(FeatureID, Sample, Intensity, imputed)], runs, by = "Sample")
}
metrics <- function(L) {
  d <- L[, .(dm = mean(Intensity), vw = if (.N < 2) NA_real_ else var(Intensity),
             n = .N, n_imp = sum(imputed)), by = .(FeatureID, SampleID, Day)]
  d[, .(cv_published = if (.N < 2) NA_real_ else 100*sd(dm)/mean(dm),
        cv_corrected = if (.N < 3 || !all(n == N_REP)) NA_real_ else
                       100*sqrt(max(0, var(dm) - mean(vw)/N_REP))/mean(dm),
        days = .N, runs = sum(n), balanced = all(n == N_REP),
        dfi = sum(n_imp == n)), by = .(FeatureID, SampleID)]
}
# eligible population for a definition (before applying the 20% cut)
elig <- function(m, def) {
  if (def == "published") m[days == 3 & is.finite(cv_published)]
  else m[days == 3 & runs == 24 & balanced == TRUE & is.finite(cv_corrected)]
}
cvcol <- function(def) if (def == "published") "cv_published" else "cv_corrected"
# One protein x sample's CV under a definition. "Combined" is the UNION of the
# two proteases -- a protein counts if either protease quantifies it well -- and
# the CV shown is the one from the protease you would actually use (the lower of
# the two). The proteases are never averaged.
pdata <- function(pr, def, which = "b", tier = "strict") {
  cc <- cvcol(def)
  one <- function(p) {
    e <- elig(M[[p]][[which]], def)
    if (which == "a") e <- e[tier_ok(e, tier)]
    e[, .(FeatureID, SampleID, CV = get(cc))]
  }
  if (pr != "Combined") return(one(pr)[, .(SampleID, CV)])
  m <- merge(one("Krakatoa"), one("Trypsin"), by = c("FeatureID","SampleID"),
             all = TRUE, suffixes = c("_k","_t"))
  m[, CV := pmin(CV_k, CV_t, na.rm = TRUE)]
  m[is.finite(CV), .(SampleID, CV)]
}
tier_ok <- function(m, tier) switch(tier, strict = m$dfi == 0, onefill = m$dfi <= 1,
                                    partial = m$dfi <= 2, ungated = rep(TRUE, nrow(m)))

cat("Loading...\n")
hta <- load_pg(file.path(REPRO, "Krakatoa Library", "report.pg_matrix.tsv"))
try <- load_pg(file.path(REPRO, "Trypsin Library",  "report.pg_matrix.tsv"))
runs <- unique(rbind(hta, try)[, .(Sample, Day, SampleID, Replicate)])
M <- list(Krakatoa = list(b = metrics(hta), a = metrics(impute(hta, try, make_offsets(hta, try), runs))),
          Trypsin  = list(b = metrics(try), a = metrics(impute(try, hta, make_offsets(try, hta), runs))))
samples <- as.character(sort(as.integer(unique(runs$SampleID))))
vv <- function(d) { r <- d[, .N, by = SampleID]; x <- r$N[match(samples, r$SampleID)]
                    x[is.na(x)] <- 0L; as.integer(x) }

# gate: Panel A
PA <- list(Trypsin = c(120L,107L,127L,123L,123L,123L,146L,131L),
           Krakatoa = c(92L,80L,88L,89L,108L,101L,91L,82L))
cat("\n--- Panel A reproduction gate ---\n")
for (pr in names(PA)) {
  e <- elig(M[[pr]]$b, "published"); v <- vv(e[cv_published <= CV_THRESHOLD])
  cat(sprintf("  %-9s %s  %s\n", pr, paste(v, collapse=" "),
              if (identical(v, PA[[pr]])) "MATCH" else "MISMATCH"))
}

# ------------------------------------------------------------- counts -------
allc <- rbindlist(lapply(c("published","corrected"), function(def)
  rbindlist(lapply(PROTEASES, function(pr) rbindlist(c(
    list(data.table(Definition=def, Protease=pr, Tier="Measured only", SampleID=samples,
                    n = vv(pdata(pr, def, "b")[CV <= CV_THRESHOLD]))),
    lapply(names(TIERS), function(tr)
      data.table(Definition=def, Protease=pr, Tier=tr, SampleID=samples,
                 n = vv(pdata(pr, def, "a", tr)[CV <= CV_THRESHOLD])))))))))
fwrite(allc, file.path(TBL, "manuscript_figure_counts.csv"))

cat("\n=== Counts (total across 8 samples) ===\n")
print(dcast(allc[, .(total = sum(n)), by = .(Definition, Protease, Tier)],
            Definition + Tier ~ Protease, value.var = "total"))

