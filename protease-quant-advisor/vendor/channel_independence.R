# channel_independence.R — the control this whole phase rests on.
#
# The premise of Phase 3 is that trypsin and HTA are two DIFFERENT views of the same injection, so
# choosing between them is a real choice and borrowing across them is real information. Both halves
# of that premise fail where the two searches actually recovered the SAME peptides: then the channels
# are not independent evidence, the recommendation is near-moot, and cross-protease imputation is
# closer to copying a number than to borrowing an independent measurement.
#
# So we measure it rather than assert it:
#   shared_peptide_frac   |peptides seen in both channels| / |union|, over qualifying peptides
#                         (a *_present twin is emitted over all detected peptides). HIGH => the two
#                         roll-ups are built from overlapping inputs; flagged in the output.
#   det_phi, det_jaccard  per-sample detection agreement (vendored detection_concordance)
#   complementarity       fraction of samples where EXACTLY ONE channel is observed. This is exactly
#                         the pool Part B can recover, so it sizes the imputation prize before we
#                         build the imputer.
#
# Also reported dataset-wide: how many DIA-NN Precursor.Ids literally collide between the two
# searches. The frozen Stage-1 grouped its per-precursor metrics by precursor_id alone, so colliding
# ids had their two channels' metrics merged; the count says how much that could have mattered.
#
# Depends on: data.table; vendored detection.R (detection_concordance, .phi_jaccard).

suppressPackageStartupMessages({ library(data.table) })

# Per-protein peptide-sequence overlap between the two channels.
peptide_overlap <- function(prec, qual_ids, labelA = "trypsin", labelB = "hta") {
  d <- prec[is.finite(log2q)]
  d[, is_qual := precursor_id %in% qual_ids]
  ov <- function(sub) {
    a <- unique(sub[protease == labelA, stripped_sequence])
    b <- unique(sub[protease == labelB, stripped_sequence])
    u <- union(a, b); i <- intersect(a, b)
    list(n_pep_a = length(a), n_pep_b = length(b), n_pep_shared = length(i),
         shared_peptide_frac = if (length(u)) length(i) / length(u) else NA_real_)
  }
  all_ov  <- d[, ov(.SD), by = pg_canonical,
               .SDcols = c("protease","stripped_sequence")]
  setnames(all_ov, c("n_pep_a","n_pep_b","n_pep_shared","shared_peptide_frac"),
           c("n_pep_present_a","n_pep_present_b","n_pep_present_shared","shared_peptide_frac_present"))
  qual_ov <- d[is_qual == TRUE, ov(.SD), by = pg_canonical,
               .SDcols = c("protease","stripped_sequence")]
  setnames(qual_ov, c("n_pep_a","n_pep_b","n_pep_shared"),
           c("n_pep_qual_a","n_pep_qual_b","n_pep_qual_shared"))
  merge(all_ov, qual_ov, by = "pg_canonical", all = TRUE)
}

# Per-protein per-sample detection agreement + complementarity, from the two channel long tables
# (rollup_to_channel_long output: one row per detected pg x sample).
detection_agreement <- function(A_long, B_long, pgs, samples) {
  dc <- detection_concordance(A_long[pg_canonical %in% pgs], B_long[pg_canonical %in% pgs], pgs, samples)
  ns <- length(samples)
  a <- A_long[pg_canonical %in% pgs & is.finite(log2q), .(pg_canonical, sample_id, a = TRUE)]
  b <- B_long[pg_canonical %in% pgs & is.finite(log2q), .(pg_canonical, sample_id, b = TRUE)]
  m <- merge(a, b, by = c("pg_canonical","sample_id"), all = TRUE)
  m[is.na(a), a := FALSE][is.na(b), b := FALSE]
  comp <- m[, .(n_samples_both = sum(a & b),
                n_samples_only_a = sum(a & !b),
                n_samples_only_b = sum(!a & b)), by = pg_canonical]
  comp[, `:=`(complementarity = (n_samples_only_a + n_samples_only_b) / ns,
              n_samples_universe = ns)]
  merge(dc, comp, by = "pg_canonical", all = TRUE)
}

# Dataset-level: how many precursor ids are literally the same string in both searches.
precursor_id_collisions <- function(prec, labelA = "trypsin", labelB = "hta") {
  a <- unique(prec[protease == labelA, precursor_id])
  b <- unique(prec[protease == labelB, precursor_id])
  i <- intersect(a, b)
  data.table(n_precursor_ids_a = length(a), n_precursor_ids_b = length(b),
             n_precursor_ids_colliding = length(i),
             frac_of_a_colliding = if (length(a)) length(i)/length(a) else NA_real_,
             frac_of_b_colliding = if (length(b)) length(i)/length(b) else NA_real_)
}

# Assemble the full independence table for a condition.
channel_independence <- function(prec, qual_ids, A_long, B_long, pgs, samples,
                                 labelA = "trypsin", labelB = "hta", high_overlap = 0.50) {
  po <- peptide_overlap(prec[pg_canonical %in% pgs], qual_ids, labelA, labelB)
  da <- detection_agreement(A_long, B_long, pgs, samples)
  out <- merge(po, da, by = "pg_canonical", all = TRUE)
  out[, shared_peptide_frac := fifelse(
    is.finite(n_pep_qual_shared) & (n_pep_qual_a + n_pep_qual_b - n_pep_qual_shared) > 0,
    n_pep_qual_shared / (n_pep_qual_a + n_pep_qual_b - n_pep_qual_shared), NA_real_)]
  out[, flag_low_independence := is.finite(shared_peptide_frac) & shared_peptide_frac >= high_overlap]
  setorder(out, pg_canonical)
  out[]
}
