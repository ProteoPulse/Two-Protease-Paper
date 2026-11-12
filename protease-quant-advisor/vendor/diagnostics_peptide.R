# diagnostics_peptide.R — OPTIONAL peptide-level audit for discordant proteins (from CP05).
# Activated only when precursor (pr) matrices are supplied for both channels.
# SURROGATE DIAGNOSTIC ONLY: never replaces the reported protein quantity or the frozen class.
#
#   1. Load pr matrices; proteotypic filter; collapse charge/precursor variants -> unique
#      Stripped.Sequence per (protein, sample) via median log2 (record n precursors).
#   2. core peptide = observed in >= core_frac of samples, >= min_obs paired, >= 5 distinct, var>0.
#   3. Sibling coherence: z each core peptide; sibling consensus = median of the others;
#      peptide-to-sibling Spearman; candidate inconsistent iff rho<=rmax AND remaining coherence>=cohmin.
#   4. Leave-one-sequence-out surrogate: surrogate profile = median of core peptide z-profiles;
#      cross-channel surrogate rho; remove each sequence, recompute; candidate_single_peptide_influence
#      iff removal crosses the frozen threshold AND raises rho >= min_delta_rho.
#
# Depends on: utils_ids.R (canonicalize_pg), utils_stats.R (to_log2).

suppressPackageStartupMessages({ library(data.table) })

load_peptides_long <- function(pr_path, sample_cols, already_log = NA, id_fn = identity) {
  pr <- fread(pr_path, sep = "\t", na.strings = c("","NA","NaN"), check.names = FALSE)
  keyc <- intersect(c("Protein.Group","Stripped.Sequence","Precursor.Id","Proteotypic"), names(pr))
  runs <- intersect(sample_cols, names(pr))
  long <- melt(pr[, c(keyc, runs), with = FALSE], id.vars = keyc,
               measure.vars = runs, variable.name = "original_column", value.name = "raw")
  long[, original_column := as.character(original_column)]
  long[, sample_id := vapply(original_column, function(x) as.character(id_fn(x)), character(1))]
  if ("Proteotypic" %in% names(long)) long <- long[Proteotypic == 1 | Proteotypic == TRUE]
  long[, log2q := to_log2(raw, already_log = already_log)]
  long <- long[is.finite(log2q)]
  long[, .(log2q = median(log2q), n_prec = .N),
       by = .(pg_canonical = canonicalize_pg(Protein.Group), stripped = Stripped.Sequence, sample_id)]
}

# sample x peptide (log2) matrix of CORE peptides for one protein/channel
core_matrix <- function(pep, pg, core_frac, min_obs) {
  d <- pep[pg_canonical == pg]
  if (!nrow(d)) return(NULL)
  n_samp <- d[, uniqueN(sample_id)]
  w <- dcast(d, sample_id ~ stripped, value.var = "log2q")
  m <- as.matrix(w[, -1, drop = FALSE]); rownames(m) <- w$sample_id
  keep <- apply(m, 2, function(v) {
    o <- is.finite(v); sum(o) >= max(min_obs, core_frac * n_samp) &&
      length(unique(v[o])) >= 5 && stats::var(v[o]) > 0 })
  m <- m[, keep, drop = FALSE]
  if (ncol(m) == 0) return(NULL)
  m
}
zcols <- function(m) apply(m, 2, function(v){ o<-is.finite(v); if(sum(o)<2||sd(v[o])==0) return(rep(NA,length(v))); (v-mean(v[o]))/sd(v[o]) })
median_pairwise_coherence <- function(m) {
  if (ncol(m) < 2) return(NA_real_)
  cc <- suppressWarnings(cor(m, method = "spearman", use = "pairwise.complete.obs"))
  median(cc[upper.tri(cc)], na.rm = TRUE)
}
surrogate_profile <- function(mz) apply(mz, 1, function(r) if (all(is.na(r))) NA_real_ else median(r, na.rm = TRUE))
xrho <- function(sa, sb, min_paired = 15) { ok <- is.finite(sa) & is.finite(sb)
  if (sum(ok) < min_paired) return(NA_real_); suppressWarnings(cor(sa[ok], sb[ok], method = "spearman")) }

# Run the full peptide audit over a set of discordant protein groups.
peptide_diagnostics <- function(pepA, pepB, discordant_pgs, threshold,
                                core_frac = 0.50, min_obs = 20,
                                rmax = 0.20, cohmin = 0.50, min_delta_rho = 0.20,
                                labelA = "trypsin", labelB = "hta",
                                max_sibling_core = 60L) {
  # max_sibling_core: skip the per-peptide sibling-candidate scan for proteins with more core
  # peptides than this. That scan recomputes a pairwise-coherence matrix inside a per-peptide
  # loop (O(k^3) in #core peptides), which explodes for hyper-abundant proteins (e.g. transferrin
  # ~286 HTA core peptides). Single-peptide-candidate analysis is meaningless for such
  # many-peptide proteins anyway, and this output feeds neither the tier labels nor any figure.
  support <- list(); coherence <- list(); candidates <- list(); loso <- list()
  for (pg in discordant_pgs) {
    mt <- core_matrix(pepA, pg, core_frac, min_obs); mk <- core_matrix(pepB, pg, core_frac, min_obs)
    nt <- if (is.null(mt)) 0L else ncol(mt); nk <- if (is.null(mk)) 0L else ncol(mk)
    st <- function(n) fifelse(n==0,"no_core",fifelse(n==1,"one_core",fifelse(n==2,"two_core","three_or_more_core")))
    support[[pg]] <- data.table(pg_canonical=pg, n_core_A=nt, n_core_B=nk,
                                support_A=st(nt), support_B=st(nk))
    coherence[[pg]] <- data.table(pg_canonical=pg,
      coherence_A = if (nt>=2) median_pairwise_coherence(mt) else NA_real_,
      coherence_B = if (nk>=2) median_pairwise_coherence(mk) else NA_real_)
    sib_one <- function(m, chan) {
      if (is.null(m) || ncol(m) < 3 || ncol(m) > max_sibling_core) return(NULL)
      mz <- zcols(m)
      rbindlist(lapply(seq_len(ncol(m)), function(j) {
        sib <- apply(mz[,-j,drop=FALSE], 1, function(r) if(all(is.na(r))) NA_real_ else median(r,na.rm=TRUE))
        data.table(pg_canonical=pg, channel=chan, stripped=colnames(m)[j],
                   sibling_rho=xrho(mz[,j], sib),
                   remaining_coherence=median_pairwise_coherence(m[,-j,drop=FALSE]))
      }))
    }
    cand <- rbind(sib_one(mt,labelA), sib_one(mk,labelB))
    if (!is.null(cand) && nrow(cand)) {
      cand[, candidate_inconsistent := is.finite(sibling_rho) & sibling_rho <= rmax &
             is.finite(remaining_coherence) & remaining_coherence >= cohmin]
      candidates[[pg]] <- cand
    }
    if (!is.null(mt) && !is.null(mk) && ncol(mt) >= 1 && ncol(mk) >= 1) {
      surT <- surrogate_profile(zcols(mt)); names(surT) <- rownames(mt)
      surK <- surrogate_profile(zcols(mk)); names(surK) <- rownames(mk)
      common <- intersect(names(surT), names(surK))
      base_rho <- xrho(surT[common], surK[common])
      loso_one <- function(m, chan, other_sur) {
        if (ncol(m) < 3) return(NULL)
        rbindlist(lapply(seq_len(ncol(m)), function(j) {
          s <- surrogate_profile(zcols(m[,-j,drop=FALSE])); names(s) <- rownames(m)
          cm <- intersect(names(s), names(other_sur)); rw <- xrho(s[cm], other_sur[cm])
          data.table(pg_canonical=pg, channel=chan, stripped=colnames(m)[j],
                     baseline_surrogate_rho=base_rho, rho_without=rw, delta_rho=rw-base_rho,
                     crosses_threshold=is.finite(rw) & rw>=threshold)
        }))
      }
      lo <- rbind(loso_one(mt,labelA,surK), loso_one(mk,labelB,surT))
      if (!is.null(lo) && nrow(lo)) {
        lo[, candidate_single_peptide_influence := crosses_threshold &
             is.finite(delta_rho) & delta_rho >= min_delta_rho]
        loso[[pg]] <- lo
      }
    }
  }
  list(support = rbindlist(support, fill=TRUE),
       coherence = rbindlist(coherence, fill=TRUE),
       candidates = rbindlist(candidates, fill=TRUE),
       loso = rbindlist(loso, fill=TRUE))
}
