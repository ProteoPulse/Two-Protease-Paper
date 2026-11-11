# classify.R — apply a (recommended or human-locked) threshold to per-protein metrics.
# Binary only: concordant iff the PRIMARY statistic (stat_value; kendall tau-b by default,
# spearman rho if selected) >= threshold, else discordant. No intermediate class.
# Adds classification_uncertain (bootstrap CI straddles the threshold) and inverse_flag (< 0).
# Never mutates the underlying quantities.

suppressPackageStartupMessages({ library(data.table) })

classify_concordance <- function(metrics, threshold) {
  d <- copy(metrics[concordance_evaluable == TRUE])
  d[, concordance_threshold := threshold]
  d[, concordance_class := fifelse(stat_value >= threshold, "concordant", "discordant")]
  d[, classification_uncertain := is.finite(rho_ci_low) & is.finite(rho_ci_high) &
       rho_ci_low < threshold & rho_ci_high >= threshold]
  d[, inverse_flag := is.finite(stat_value) & stat_value < 0]
  d[]
}

# Class-consistency against an optional validation comparison (e.g. a second gradient),
# matched by pg_canonical. Returns the primary table with a `validation_support` column:
# supported / not_supported / not_assessable.
add_validation_support <- function(primary, validation, threshold) {
  if (is.null(validation) || !nrow(validation)) {
    primary[, validation_support := "not_assessable"]; return(primary[])
  }
  v <- validation[concordance_evaluable == TRUE,
                  .(pg_canonical, class_v = fifelse(stat_value >= threshold, "concordant", "discordant"))]
  m <- merge(primary, v, by = "pg_canonical", all.x = TRUE)
  m[, validation_support := fifelse(is.na(class_v), "not_assessable",
                             fifelse(class_v == concordance_class, "supported", "not_supported"))]
  m[, class_v := NULL]
  m[]
}

# Reconciliation crosswalk: compare the current primary classification against a previously
# frozen disposition table (read from disk). Joins by pg_canonical and reports every protein
# whose class changed, carrying both statistic values for transparency. `frozen_path` points to
# a prior run's final_protein_disposition.tsv. Returns NULL if the file is unreadable.
crosswalk_vs_frozen <- function(final, frozen_path) {
  if (is.null(frozen_path) || !file.exists(frozen_path)) return(NULL)
  fr <- tryCatch(fread(frozen_path, sep = "\t"), error = function(e) NULL)
  if (is.null(fr) || !all(c("pg_canonical","concordance_class") %in% names(fr))) return(NULL)
  fr_stat <- intersect(c("spearman_rho","stat_value"), names(fr))[1]
  frz <- fr[, .(pg_canonical, frozen_class = concordance_class,
                frozen_stat = if (!is.na(fr_stat)) get(fr_stat) else NA_real_,
                frozen_threshold = if ("concordance_threshold" %in% names(fr)) concordance_threshold else NA_real_)]
  cur <- final[, .(pg_canonical, new_class = concordance_class, new_stat = stat_value,
                   new_kendall_tau = if ("kendall_tau" %in% names(final)) kendall_tau else NA_real_,
                   new_spearman_rho = if ("spearman_rho" %in% names(final)) spearman_rho else NA_real_,
                   new_threshold = concordance_threshold)]
  m <- merge(cur, frz, by = "pg_canonical", all = TRUE)
  m[, flipped := !is.na(new_class) & !is.na(frozen_class) & new_class != frozen_class]
  setorder(m, -flipped, new_class)
  m[]
}

classification_sensitivity <- function(metrics, threshold, min_paired_grid,
                                       min_distinct = 5) {
  rbindlist(lapply(min_paired_grid, function(mp) {
    d <- metrics[n_paired >= mp & distinct_a >= min_distinct & distinct_b >= min_distinct &
                   is.finite(sd_a) & is.finite(sd_b) & sd_a > 0 & sd_b > 0]
    data.table(min_paired = mp, n_evaluable = nrow(d),
               n_concordant = sum(d$stat_value >= threshold, na.rm = TRUE),
               n_discordant = sum(d$stat_value <  threshold, na.rm = TRUE))
  }))
}
