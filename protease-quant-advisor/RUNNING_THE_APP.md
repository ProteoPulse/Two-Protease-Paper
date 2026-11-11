# Running the Protease Quant Advisor

A practical guide: requirements, how to run it (GUI and command line), what to feed it, what comes
out, every setting, and the nuances that decide whether a number should be trusted. For the scientific
rationale behind the defaults see `README.md` and the manuscript's Methods.

---

## 1. What it does, in one line

Given two DIA-NN searches of the **same injections** (one per protease), it returns **one protein
quantity per protein per sample** — choosing, protein by protein, which search to trust — plus the
evidence for each call and an optional cross-protease gap-fill. It never averages the two searches and
never invents a new quantity.

## 2. Requirements

- **R ≥ 4.5**
- Packages: `data.table`, `ggplot2`, `yaml`, `base64enc` (engine); `shiny`, `DT` (GUI only).

```r
install.packages(c("data.table","ggplot2","yaml","shiny","DT","base64enc"))
```

No compilation, no internet, no external services. A discovery-scale run (≈ 300 samples) completes in
a few minutes on a laptop; memory scales with the pr_matrix size.

## 3. Inputs — exactly four files

One pair per search, all standard DIA-NN outputs:

| input | what it is | why it is needed |
|---|---|---|
| `pr_a` | trypsin `report.pr_matrix.tsv` | precursor evidence (qualification, coherence, agreement) |
| `pr_b` | HTA/Krakatoa `report.pr_matrix.tsv` | same, other search |
| `pg_a` | trypsin `report.pg_matrix.tsv` | the protein quantities being adjudicated |
| `pg_b` | HTA/Krakatoa `report.pg_matrix.tsv` | same, other search |

Requirements on the files:
- Both searches must be of the **same raw files**, so the run (column) names match between `pr` and
  `pg` and between the two searches. The tool intersects on sample columns; a run only one search saw
  is dropped from the comparison.
- Column headers are raw-file paths; the tool reduces them to a bare run name (basename, extension
  removed) before matching, so Windows/Unix path differences do not matter.
- The `pg_matrix` may have any filename (e.g. `trypsin_library_search_report.pg_matrix.tsv`); pass the
  actual path.

## 4. Running it

### GUI

```bash
Rscript -e 'shiny::runApp("app", launch.browser = TRUE)'
```

Load the four matrices, set the run filters, press **Run advisor**. (Stage 6 and its settings are
command-line only.) Each stage has a tab with its
decision plot and a downloadable table; the full table set downloads as a zip and the report as one
self-contained HTML file.

### Command line (same engine, same numbers)

```bash
Rscript run_advisor.R \
  pr_a=<trypsin report.pr_matrix.tsv>  pr_b=<hta report.pr_matrix.tsv> \
  pg_a=<trypsin report.pg_matrix.tsv>  pg_b=<hta report.pg_matrix.tsv> \
  run_include=Plasma_Mix run_exclude=20min \
  selection_mode=completeness impute=true offset_model=protein_median \
  tag=mix30 out=results/mix30
```

Worked examples for every condition are in `run_all_conditions.sh`; scenarios with expected results are
in `tests/SCENARIOS.md`.

## 5. Settings

All defaults live in `config_default.yml`, each annotated with the reason it was chosen. The settings
listed in `run_advisor.R`'s usage header can be set on the command line as `key=value`; anything else is
changed with `config=my.yml`, a YAML file holding only the keys to change, layered over the defaults.

**Run selection / ingest**
| key | default | meaning |
|---|---|---|
| `run_include` | "" (all) | regex on the normalized run name; keep only matching runs |
| `run_exclude` | "" | regex; drop matching runs (e.g. `20min`) |
| `run_sample_regex` | "" | capture group mapping a run name to a biological sample id; blank = one run per sample. **Set this when you have technical replicates**, so replicate-based criteria (CV, ICC) switch on. |

**Channel selection (Stage 3)**
| key | default | meaning |
|---|---|---|
| `selection_mode` | `completeness` | `completeness` = the manuscript method: the channel that quantified the protein in more samples (higher mean precursor detection rate) wins, held per protein. `multi_criteria` = the weighted-vote generalization (completeness/coherence/split-half wt 2, depth wt 1, and with replicates technical CV wt 2 + ICC wt 1). Both are always computed; the unused rule's answer is reported as `vote_channel` and `selection_disagrees_with_vote` flags where they differ. |

Completeness-mode confidence tiers on the size of the detection-rate gap
(`config_default.yml → recommend`: `completeness_high_gap` 0.20, `completeness_moderate_gap` 0.05); the
manuscript defines high/moderate/low only for the weighted vote, so these are an app convenience.

**Imputation (Stage 5)**
| key | default | meaning |
|---|---|---|
| `impute` | false | turn gap-filling on. Off by default because on well-measured data it fires on few cells, by design |
| `offset_model` | `protein_median` | `protein_median` (the protein's median offset, as in the manuscript) or `sample_knn` (offset from the k nearest co-observed samples; ~0.04 log2 more accurate on discovery, near-identical on replicate data) |
| `min_paired` / `n_perm` | 15 / 1000 | agreement (Stage 4) paired-sample minimum and permutation count |

Deeper imputation gates (`config_default.yml → impute`): donor completeness ≥ 0.60, ≥ 8 co-observed
runs, CCC ≥ 0.50 or τ ≥ 0.30, offset MAD ≤ 1.0 log2, discordant proteins hard-blocked, borrowed
offsets off. These are set by benchmark; change them only deliberately.

**Reproducibility / well-quantified counts (Stage 6)**
| key | default | meaning |
|---|---|---|
| `run_day_regex` | "" | regex whose first group captures a preparation-day label from the run name. With replicates **and** a day map → corrected inter-day CV; with replicates but no map → between-replicate SD; with no replicates → cohort completeness (CV is suppressed, not fabricated). For the reproducibility cohort (`Day_Sample_Rep` run names) use `run_day_regex='^([0-9]+)_.*$'`. |
| `cv_threshold` | 20 | % CV below which a protein × sample counts as well quantified (CV modes) |
| `completeness_min` | 0.50 | fraction of runs a protein must be quantified in (completeness mode) |

**Other**: `fasta=` enables sequence-coverage criteria; `probe_in_score=true` lets the roll-up QC probe
vote (off by default — it correlates 0.91 with peptide count and changes no consequential calls);
`tag=` / `out=` name the run and its output directory.

## 6. Outputs

Written to `out/` (default `results/<tag>/`):

| file | contents |
|---|---|
| **`final_abundance_matrix.tsv`** | **the deliverable** — one value per protein per sample, with source channel, confidence, and measured-vs-imputed provenance |
| `recommendation.tsv` | per shared protein: chosen channel, confidence, margin, scores, criteria used, consequential flag |
| `precursor_evidence.tsv` | every evidence metric, per protein and per protease |
| `concordance.tsv` | per-protein agreement statistic, threshold, class |
| `precursor_qualification.tsv`, `qualification_sensitivity.tsv` | per-precursor QC and the threshold sweep |
| `detection_summary.tsv`, `ingest_summary.tsv` | shared vs single-channel counts; what was read after filtering |
| `imputation_ledger.tsv`, `imputation_refused.tsv`, `imputation_per_protein.tsv` | every filled cell (donor, offset, gates, tier), every refused gap with its reason, per-protein gate outcomes |
| `reproducibility_wellquantified` (Stage 6) | well-quantified counts per protease and union, with the metric used |
| `report.html` | self-contained stage-by-stage report with plots |
| `INTERPRETATION_AND_METHODS.md` | a methods/interpretation document generated from this run's own numbers |
| `advisor_state.rds` | the complete run state; re-render reports without recomputing |

**Reading the deliverable.** Key columns of `final_abundance_matrix.tsv`:
`recommended_channel` (constant per protein), `source_channel` (= recommended_channel, the scale the
value lives on), `value_source` (`measured` or `imputed_cross_protease`), and for imputed cells
`donor_channel`, `donor_value`, `offset_used`. To use only measured data, filter
`value_source == "measured"`; to keep only confident calls, filter `recommendation_confidence`; to see
only the proteins where the choice mattered, filter `consequential == TRUE`.

## 7. Nuances that decide whether to trust a number

1. **The decision is per protein, held across all samples** — not per sample (which would put the
   between-search offset inside a protein's profile and manufacture fold-changes) and not one channel
   for the whole study. This is enforced: the run halts if any protein ever carries two source
   channels, and `build_index.R` re-checks it (`proteins_with_two_channels` must be 0).
2. **The two proteases are never averaged.** Combined results are the union or a single chosen channel.
3. **The pr_matrix is pre-filtered.** DIA-NN's 1% q-value filter is inherited as a floor and cannot be
   loosened; the tool states this in the run header. Its own QC (detection, replicate CV, coherence)
   is computed from intensities and adds on top.
4. **Most gap-filling requests are refused, by design** — a donor that does not track the recipient is
   not a surrogate. `imputation_refused.tsv` is usually the more informative table; a gap left empty is
   a correct decision, not a shortfall.
5. **Imputation is a refinement, not the coverage story.** It substantially improves a sparse single
   channel but adds little to the union of both proteases; report it as such.
6. **Stage 6 picks its metric from the data.** Do not read an inter-day CV where there are no
   replicates; with one run per sample the stage reports completeness and suppresses CV rather than
   fabricating one. Supply `run_sample_regex` (replicates) and `run_day_regex` (days) to unlock the
   precision metrics.

## 8. Validation

`Rscript tests/run_all_tests.R` runs the whole suite. `tests/test_example_end_to_end.R` runs every
stage on the bundled synthetic data and checks the tool's invariants, and `tests/test_shiny_app.R`
drives the GUI headlessly on the same data; neither needs study data, and both run in CI. `tests/test_manuscript_defaults.R`
asserts the shipped defaults match the manuscript and that completeness selection reproduces an
independent "more complete channel" recomputation; `tests/test_ingest.R` checks the reader against known
ground truth (54 runs per search and 392 shared proteins on Mix30);
`tests/test_reproducibility_stage.R` asserts Stage 6 reproduces the published Panel A counts in
inter-day-CV mode and falls back to completeness without replicates. The real-data tests need
`TWO_PROTEASE_DATA` set to the MassIVE search results and skip without it. `tests/SCENARIOS.md` lists
copy-paste scenarios with expected numbers. `vendor/MANIFEST.yml` holds md5s of the reused modules so
drift from the source frameworks stays detectable.
