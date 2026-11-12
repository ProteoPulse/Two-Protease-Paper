# modification_inventory.R: every UniMod modification present in the DIA-NN libraries and
# precursor-filtered reports (Q.Value and Lib.Q.Value <= 0.01) of the four searches, with names,
# for the MassIVE "post-translational modifications" field.
# Usage: Rscript figure-scripts/07_verification/modification_inventory.R
suppressPackageStartupMessages({ library(duckdb); library(data.table) })
.sd <- function() { a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f)) normalizePath(dirname(f)) else getwd() }
source(file.path(.sd(), "..", "config.R"))   # DATA_DIR, OUT_ROOT, step_dirs()
B <- DATA_DIR
sets <- c(disc_krakatoa = "Krakatoa_Library_Search", disc_trypsin = "Trypsin_Library_Search",
          repro_krakatoa = "Reproducibility_Search/Krakatoa Library",
          repro_trypsin = "Reproducibility_Search/Trypsin Library")
con <- dbConnect(duckdb())
res <- rbindlist(lapply(names(sets), function(nm) rbindlist(lapply(c("lib", "report"), function(kind) {
  f <- file.path(B, sets[[nm]], paste0(kind, ".parquet"))
  filt <- if (kind == "report") "WHERE \"Q.Value\" <= 0.01 AND \"Lib.Q.Value\" <= 0.01" else ""
  s <- dbGetQuery(con, sprintf("SELECT DISTINCT \"Modified.Sequence\" s FROM read_parquet('%s') %s", f, filt))$s
  m <- regmatches(s, gregexpr("[A-Z]?\\(UniMod:[0-9]+\\)", s))
  dt <- data.table(seq = rep(s, lengths(m)), tok = unlist(m))
  dt[, `:=`(id = as.integer(sub(".*UniMod:([0-9]+)\\)", "\\1", tok)),
            site = fifelse(substr(tok, 1, 1) == "(", "N-term", substr(tok, 1, 1)))]
  dt[, .(n_peptides = uniqueN(seq)), by = .(id, site)][, `:=`(set = nm, kind = kind)]
}))))
nm_tbl <- fread(file.path(FIGSCRIPTS_DIR, "03_ptm_coverage", "unimod_name_library.tsv"))
# names missing from unimod_name_library.tsv (the four most common types)
nm_tbl <- rbind(nm_tbl, data.table(UnimodID = c(4L, 35L, 7L, 2L),
                UnimodName = c("Carbamidomethyl", "Oxidation", "Deamidated", "Amidated")))[!duplicated(UnimodID)]
rep <- res[kind == "report"]
# main sites = residues carrying >= 10% of a type's occurrences (site-level peptide counts)
sites <- rep[, .(n = sum(n_peptides)), by = .(id, site)][, frac := n / sum(n), by = id][
  order(-frac, site)][, .(main_sites = paste(site[frac >= 0.10], collapse = ","), n_sites = .N), by = id]   # ties broken by site name (deterministic)
w <- dcast(rep[, .(n_peptides = sum(n_peptides)), by = .(id, set)], id ~ set, value.var = "n_peptides", fill = 0)
w <- merge(merge(w, sites, by = "id"), nm_tbl, by.x = "id", by.y = "UnimodID", all.x = TRUE)
w[, total_peptides := disc_krakatoa + disc_trypsin + repro_krakatoa + repro_trypsin]
w[, max_any_search := pmax(disc_krakatoa, disc_trypsin, repro_krakatoa, repro_trypsin)]
setorder(w, -total_peptides)
w[, unimod_accession := paste0("UNIMOD:", id)]
setcolorder(w, c("unimod_accession", "UnimodName", "main_sites", "n_sites", "total_peptides", "max_any_search"))
fwrite(w, file.path(step_dirs("07_verification")$tables, "modifications_in_reports.csv"))
cat("Distinct UniMod types | libraries:", uniqueN(res[kind == "lib"]$id),
    "| filtered reports:", nrow(w), "| >=10 peptides in a search:", w[max_any_search >= 10, .N],
    "| unnamed:", w[is.na(UnimodName), .N], "\n")
print(w[1:40, .(unimod_accession, UnimodName, main_sites, n_sites, total_peptides, max_any_search)])
