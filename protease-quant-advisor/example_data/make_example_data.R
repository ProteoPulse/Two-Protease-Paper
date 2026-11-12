#!/usr/bin/env Rscript
# make_example_data.R — a small SYNTHETIC two-protease dataset, so the Protease Quant Advisor and its
# tests run on a fresh clone without the MassIVE data.
#
# Nothing here is real. Protein accessions (SYN00001...), genes, sequences and intensities are all
# simulated. The design mimics the reproducibility cohort of the manuscript at a much smaller scale:
#
#   16 biological samples x 3 preparation days x 1 injection = 48 runs, named {day}_{sample}_{rep}
#   (e.g. 2_07_1), with full Windows-style raw-file paths as column headers, as DIA-NN writes them.
#
#   220 proteins: 140 seen by both searches, 50 trypsin-only, 30 HTA-only.
#   Each channel measures the shared truth plus its own per-protein offset, precursor effects and
#   noise; HTA is noisier and sparser (fewer precursors), as in the real data. 15 shared proteins are
#   planted as DISCORDANT (the HTA signal does not track the truth) so the concordance stage has
#   something to find. Values go missing not-at-random (low intensities drop out more), which leaves
#   gaps in the chosen channel; a planted group of proteins also loses trypsin's protein-level value
#   in a quarter of the samples while HTA still covers them, so the imputation stage has work to do,
#   and another group is better seen by HTA, so both channels win some calls.
#
# Usage:  Rscript example_data/make_example_data.R      (rewrites the four .tsv files next to it)

suppressPackageStartupMessages({ library(data.table) })
set.seed(20261110)

.here <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else "example_data"
}
OUT <- .here()

n_samples <- 16L; days <- 1:3; rep_id <- 1L
runs <- CJ(day = days, sample = seq_len(n_samples))
runs[, run := sprintf("%d_%02d_%d", day, sample, rep_id)]
run_paths <- sprintf("C:\\Data\\TwoProtease\\Example\\%s.raw", runs$run)

n_shared <- 140L; n_try_only <- 50L; n_hta_only <- 30L
n_prot <- n_shared + n_try_only + n_hta_only
prot <- data.table(
  idx = seq_len(n_prot),
  acc = sprintf("SYN%05d", seq_len(n_prot)),
  gene = sprintf("SYNG%03d", seq_len(n_prot)),
  class = c(rep("shared", n_shared), rep("trypsin_only", n_try_only), rep("hta_only", n_hta_only)))
prot[, mu := rnorm(.N, 21, 1.8)]                         # protein abundance, log2
prot[, bio_sd := runif(.N, 0.4, 1.2)]                    # between-sample biological spread
prot[, discordant := FALSE][class == "shared", discordant := seq_len(.N) %in% sample(.N, 15L)]
# per-protein, per-protease offsets (different ionisation of different peptides)
prot[, off_try := rnorm(.N, 0, 0.8)][, off_hta := rnorm(.N, -1.5, 1.2)]
# Two planted groups among the concordant shared proteins:
#   try_weak : trypsin sees them poorly, HTA well -> HTA is the more complete channel and gets chosen
#   gappy    : trypsin is chosen, but its PROTEIN-level value is missing for a quarter of the samples
#              (as DIA-NN's protein q-value filter does) while HTA covers them -> imputation can fill
ok_shared <- prot[class == "shared" & !discordant, idx]
weak <- sample(ok_shared, 25L); gappy <- sample(setdiff(ok_shared, weak), 30L)
prot[, `:=`(try_weak = idx %in% weak, gappy = idx %in% gappy)]
prot[try_weak == TRUE, `:=`(off_try = off_try - 3.5, off_hta = 0.5)]
prot[gappy == TRUE, `:=`(mu = mu + 1, off_hta = 0)]

# true sample-level abundance; a small day effect; independent run noise per channel later
truth <- CJ(idx = prot$idx, sample = seq_len(n_samples))
truth[prot, on = "idx", `:=`(mu = i.mu, bio_sd = i.bio_sd)]
truth[, bio := rnorm(.N, 0, bio_sd)]
truth <- truth[runs, on = "sample", allow.cartesian = TRUE]
truth[, day_eff := rnorm(.N, 0, 0.10)]
truth[, value := mu + bio + day_eff]

AA <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
rand_seq <- function(n, tryptic) vapply(seq_len(n), function(i) {
  core <- paste(sample(setdiff(AA, c("K", "R")), sample(6:16, 1), replace = TRUE), collapse = "")
  if (tryptic) paste0(core, sample(c("K", "R"), 1)) else paste0(core, sample(AA, 1))
}, character(1))

make_channel <- function(channel) {
  is_try <- channel == "trypsin"
  keep <- if (is_try) prot[class != "hta_only"] else prot[class != "trypsin_only"]
  # precursors per protein: trypsin richer than HTA
  keep[, n_prec := if (is_try) pmin(rpois(.N, 5) + 1L, 12L) else pmin(rpois(.N, 2) + 1L, 6L)]
  pr <- keep[rep(seq_len(.N), n_prec), .(idx, acc, gene, discordant,
                                         off = if (is_try) off_try else off_hta)]
  pr[, prec_eff := rnorm(.N, 0, 1.0)]
  pr[, seq := rand_seq(.N, is_try)]
  pr[, charge := sample(2:3, .N, replace = TRUE, prob = c(0.7, 0.3))]
  pr[, prec_id := paste0(seq, charge)]
  pr <- unique(pr, by = "prec_id")
  pr[, pid := .I]

  noise_sd <- if (is_try) 0.20 else 0.35
  L <- pr[, .(pid, idx, discordant, off, prec_eff)][truth[, .(idx, sample, day, run, value)],
            on = "idx", allow.cartesian = TRUE, nomatch = NULL]
  # a discordant protein's HTA signal follows its own random sample profile, not the truth
  if (!is_try) {
    rnd <- unique(L[discordant == TRUE, .(idx, sample)])[, scramble := rnorm(.N, 0, 1.0)]
    L[rnd, on = .(idx, sample), value := value - (value - mean(value)) + i.scramble]
  }
  L[, x := value + off + prec_eff + rnorm(.N, 0, noise_sd)]
  # missing not at random: the lower the value, the likelier it drops out
  lod <- if (is_try) 18.5 else 19.5
  L[, p_miss := plogis((lod - x) / 0.6) + if (is_try) 0.02 else 0.06]
  L[runif(.N) < p_miss, x := NA_real_]

  # precursor matrix (linear intensities, as DIA-NN writes them)
  W <- dcast(L, pid ~ run, value.var = "x")
  setcolorder(W, c("pid", runs$run))
  for (r in runs$run) set(W, j = r, value = round(2^W[[r]], 1))
  ann <- pr[W$pid, .(Protein.Group = acc, Protein.Ids = acc, Protein.Names = paste0(gene, "_SYNTH"),
                     Genes = gene, First.Protein.Description = "Synthetic example protein",
                     Proteotypic = 1L, Stripped.Sequence = seq, Modified.Sequence = seq,
                     Precursor.Charge = charge, Precursor.Id = prec_id)]
  prm <- cbind(ann, W[, runs$run, with = FALSE])
  setnames(prm, runs$run, run_paths)

  # protein-group matrix: sum of the protein's detected precursors per run, plus a few protein-level
  # dropouts (DIA-NN's protein q-value filtering removes some protein x run cells independently)
  P <- L[is.finite(x), .(q = sum(2^x)), by = .(idx, run)]
  P <- P[runif(.N) > 0.04]
  if (is_try) {
    gap <- prot[gappy == TRUE, .(sample = sample(n_samples, 4L)), by = idx]
    gap <- gap[runs, on = "sample", allow.cartesian = TRUE, nomatch = NULL]
    P <- P[!gap, on = .(idx, run)]
  }
  PW <- dcast(P, idx ~ run, value.var = "q")
  for (r in setdiff(runs$run, names(PW))) PW[, (r) := NA_real_]
  setcolorder(PW, c("idx", runs$run))
  for (r in runs$run) set(PW, j = r, value = round(PW[[r]], 1))
  pa <- prot[PW$idx, .(Protein.Group = acc, Protein.Names = paste0(gene, "_SYNTH"), Genes = gene,
                       First.Protein.Description = "Synthetic example protein")]
  pgm <- cbind(pa, PW[, runs$run, with = FALSE])
  setnames(pgm, runs$run, run_paths)
  list(pr = prm, pg = pgm)
}

for (ch in c("trypsin", "hta")) {
  m <- make_channel(ch)
  fwrite(m$pr, file.path(OUT, sprintf("%s_report.pr_matrix.tsv", ch)), sep = "\t", na = "")
  fwrite(m$pg, file.path(OUT, sprintf("%s_report.pg_matrix.tsv", ch)), sep = "\t", na = "")
  cat(sprintf("%-8s precursors %5d | proteins %4d | missing pr %.0f%% | missing pg %.0f%%\n", ch,
              nrow(m$pr), nrow(m$pg),
              100 * mean(is.na(as.matrix(m$pr[, run_paths, with = FALSE]))),
              100 * mean(is.na(as.matrix(m$pg[, run_paths, with = FALSE])))))
}
cat("wrote 4 files to", OUT, "\n")
