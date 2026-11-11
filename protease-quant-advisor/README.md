# Protease Quant Advisor

A standalone R application for **two-protease DIA proteomics**. You searched the same injections twice —
once against a trypsin library, once against a second protease (here HTA / Krakatoa) — so every protein
that both searches report has **two** protein-group quantities for the same sample. This tool decides,
protein by protein, **which of the two to trust**, and returns one clean value per protein per sample.

It does **not** invent a new quantity. It never averages the two searches. It adjudicates between the
two numbers you already have, shows the evidence behind every call, and can optionally fill gaps in the
chosen channel from the other protease — always reversibly, and never for proteins where the two
searches disagree.

The defaults reproduce the method in the accompanying manuscript (Alba, Bharadwaj *et al.*, bioRxiv
2026, doi:10.64898/2026.10.05.756883); every default that was set by an experiment
rather than by convention is annotated with that experiment in `config_default.yml`.

---

## Contents

- [Install](#install)
- [Quickstart](#quickstart) · [GUI](#gui) · [command line](#command-line)
- [Inputs — exactly four files](#inputs--exactly-four-files)
- [What it does, stage by stage](#what-it-does-stage-by-stage)
- [Channel selection: `completeness` vs `multi_criteria`](#channel-selection-completeness-vs-multi_criteria)
- [The grain of the decision](#the-grain-of-the-decision-and-why-it-is-enforced)
- [Cross-protease gap filling](#cross-protease-gap-filling-optional)
- [Outputs](#outputs)
- [Settings reference](#settings-reference)
- [Verify it works](#verify-it-works)
- [Repository layout](#repository-layout)
- [Data and code availability](#data-and-code-availability)

---

## Install

Requires **R ≥ 4.5** (tested with R 4.5.2). Install the packages once:

```bash
Rscript install.R          # engine + GUI packages
Rscript install.R --engine # engine only (skip the Shiny GUI)
```

That installs `data.table`, `ggplot2`, `yaml`, `base64enc` (engine) and `shiny`, `DT` (GUI). No
compilation, no internet at run time, no external services. See `requirements.txt` for the list.

**Try it on the bundled synthetic data** (16 samples × 3 days; ~30 s):

```bash
Rscript run_advisor.R \
  pr_a=example_data/trypsin_report.pr_matrix.tsv pr_b=example_data/hta_report.pr_matrix.tsv \
  pg_a=example_data/trypsin_report.pg_matrix.tsv pg_b=example_data/hta_report.pg_matrix.tsv \
  run_sample_regex='^[0-9]+_([0-9]+)_[0-9]+$' run_day_regex='^([0-9]+)_.*$' \
  impute=true tag=example
```

The example matrices are simulated (`example_data/make_example_data.R` regenerates them); they exist
only so the tool and its tests run without downloading the study data.

---

## Quickstart

Two ways to run it, over the **same engine** (`engine/pipeline.R::run_advisor()`), so the GUI and the
command line can never drift apart.

### GUI

```bash
Rscript -e 'shiny::runApp("app", launch.browser = TRUE)'
```

Load the four matrices, set the run filters, choose the channel-selection mode, press **Run advisor**.
Stages 0–5 each get a tab with their decision plot and a downloadable table; the whole table set downloads
as a zip and the report as one self-contained HTML file. Stage 6 (well-quantified counts) and its
settings (`run_day_regex`, `cv_threshold`, `completeness_min`) are command-line only.

### Command line

```bash
Rscript run_advisor.R \
  pr_a=<trypsin report.pr_matrix.tsv>  pr_b=<hta report.pr_matrix.tsv> \
  pg_a=<trypsin report.pg_matrix.tsv>  pg_b=<hta report.pg_matrix.tsv> \
  run_include=Plasma_Mix run_exclude=20min \
  selection_mode=completeness impute=true \
  tag=mix30 out=results/mix30
```

Worked commands for every condition are in `run_all_conditions.sh`; ready-to-run scenarios with expected
results are in `tests/SCENARIOS.md`.

---

## Inputs — exactly four files

One **pair** per search, all standard DIA-NN outputs:

| input | what it is | why it is needed |
|---|---|---|
| `pr_a` | trypsin `report.pr_matrix.tsv` | precursor evidence (qualification, coherence, agreement) |
| `pr_b` | HTA/Krakatoa `report.pr_matrix.tsv` | same, other search |
| `pg_a` | trypsin `report.pg_matrix.tsv` | the protein quantities being adjudicated |
| `pg_b` | HTA/Krakatoa `report.pg_matrix.tsv` | same, other search |

- Both searches must be of the **same raw files**, so the run (column) names line up. The tool
  intersects on sample columns; a run only one search saw is dropped from the comparison.
- Column headers are raw-file paths; they are reduced to a bare run name (basename, extension removed)
  before matching, so Windows/Unix path differences do not matter.
- A `pg_matrix` may have any filename (e.g. `trypsin_library_search_report.pg_matrix.tsv`); pass the
  actual path.
- **The `pr_matrix` is already filtered** at DIA-NN's 1% run-specific and global q-value, and it does
  not retain per-run q-values. That 1% filter is inherited as a floor and cannot be loosened; the tool
  states this in the run header rather than silently re-applying it. Its own qualification (detection
  rate, replicate CV, sibling coherence) is computed from intensities and adds on top.

---

## What it does, stage by stage

| stage | question it answers | key output |
|---|---|---|
| 0 · ingest | which runs, which samples, how much did each search see? | `ingest_summary`, `detection_summary` |
| 1 · qualify | which precursors are trustworthy enough to count as evidence? | `precursor_qualification`, `qualification_sensitivity` |
| 2 · evidence | how well does each protease support this protein? | `precursor_evidence` |
| 4 · agreement | do the two channels agree about this protein? (permutation-null Kendall τ) | `concordance` |
| 3 · recommend | **trypsin or HTA, and how confident?** | `recommendation` |
| 5 · gap filling | can a gap in the chosen channel be filled from the other? *(optional)* | `imputation_ledger`, `imputation_refused` |
| 6 · reproducibility | how many proteins are *well quantified*, per protease and as the union? | `reproducibility_wellquantified` |

The headline deliverable is **`final_abundance_matrix.tsv`**: one value per protein per sample, with the
channel it came from, the confidence of that choice, and whether it was measured or imputed.

---

## Channel selection: `completeness` vs `multi_criteria`

For every protein both searches saw, the tool picks one channel to supply that protein's values,
**held constant across all samples**. Two modes:

- **`completeness` (default — the manuscript method).** The channel that quantified the protein in the
  larger fraction of samples (higher **mean precursor detection rate**) wins. This is the rule used
  throughout the manuscript and its figures: a channel that sees a protein in more runs gives the more
  usable per-sample profile and fewer gaps to fill.
- **`multi_criteria` (a generalization the app offers).** A weighted head-to-head over more evidence:
  completeness, sibling coherence and split-half reliability (weight 2 each), peptide depth (weight 1),
  and — when replicates are present — technical CV (weight 2) and replicate ICC (weight 1). It agrees
  with completeness on the clear cases and matters only for the small set of proteins where two channels
  are comparably complete but differ in precision.

**Both are always computed.** Whichever mode you pick, the *other* rule's answer is reported per protein
as `vote_channel` / `vote_confidence`, and `selection_disagrees_with_vote` flags exactly where the two
would differ — so switching modes never hides anything. Set the mode with `selection_mode=` on the CLI,
the "Channel selection" control in the GUI, or `recommend.selection_mode` in `config_default.yml`.

Every shared protein gets a call — there is no abstention. Thin evidence surfaces as
`recommendation_confidence = low`, a statement about certainty, not a refusal.

---

## The grain of the decision, and why it is enforced

The recommendation is made **once per protein and held across every sample.** Not per sample, and not
once for the whole experiment.

| level | used | why |
|---|---|---|
| per protein **per sample** | no | A protein taking its value from trypsin in one injection and HTA in the next carries the between-search offset inside its own profile and manufactures fold-changes that are pure search artefact. |
| per protein, **fixed across samples** | **yes** | Each shared protein is adjudicated once on its own evidence; that channel then supplies all of its values, so every within-protein comparison stays inside one measurement system. |
| one channel for the **whole experiment** | no | Collapses to "use trypsin" and discards the proteins the other search genuinely measures better. |

This is enforced, not assumed: `run_advisor()` **halts** if any protein ever carries two source
channels, and `build_index.R` re-checks it from each written matrix (`proteins_with_two_channels` must
be 0).

---

## Cross-protease gap filling (optional)

When the chosen channel has no value for a sample, the other protease's measurement of the **same
protein in the same injection** is translated onto the chosen channel's scale by that protein's own
log2 offset — so the written number still lives in the chosen channel's frame and `source_channel`
stays constant. Imputation is **not** an exception to the per-protein grain.

- **Offset model (default `protein_median`, as described in the manuscript).** The offset is the
  protein's median log2(chosen − other) over the samples where both proteases measured it (≥ 8
  co-observed runs required). `sample_knn` is offered as an alternative: it takes the median over the
  *k* = 10 co-observed samples nearest the target run in donor-profile space, so the translation adapts
  to local sample context. In our benchmark it was slightly more accurate on the discovery cohort
  (≈ 0.04 log2; paired *p* = 3.5 × 10⁻⁷ / 1.0 × 10⁻⁸) and near-identical on replicate data. A per-protein *slope* is deliberately not offered (best under
  MCAR, among the worst under MNAR); borrowed offsets from similar proteins are off (benchmarked worse
  than leaving the value missing).
- **Strict, and gated.** Only proteins the chosen channel measured at least once are filled. A fill must
  clear: donor completeness ≥ 0.60, ≥ 8 co-observed runs, CCC ≥ 0.50 **or** τ ≥ 0.30, offset MAD ≤ 1.0
  log2 — and any protein whose channels **disagree** is hard-blocked. Every filled cell and every
  refused gap is written to a ledger with its reason. On well-measured data most requests are refused,
  by design: a donor that does not track the recipient is not a surrogate, and `imputation_refused.tsv`
  is usually the more informative table.

Off by default; turn it on with `impute=true`.

---

## Outputs

Written to `out/` (default `results/<tag>/`):

| file | contents |
|---|---|
| **`final_abundance_matrix.tsv`** | **the deliverable** — one value per protein per sample, with source channel, confidence, and measured-vs-imputed provenance |
| `recommendation.tsv` | per shared protein: chosen channel, confidence, margin; the other mode's `vote_channel`/`vote_confidence`; `selection_disagrees_with_vote`; `consequential` flag |
| `precursor_evidence.tsv` | every evidence metric, per protein and per protease |
| `concordance.tsv` | per-protein agreement statistic, permutation-null threshold, class |
| `precursor_qualification.tsv`, `qualification_sensitivity.tsv` | per-precursor QC and the threshold sweep |
| `detection_summary.tsv`, `ingest_summary.tsv` | shared vs single-channel counts; what was read after filtering |
| `imputation_ledger.tsv`, `imputation_refused.tsv`, `imputation_per_protein.tsv` | every filled cell (donor, offset, gates), every refused gap with its reason, per-protein gate outcomes |
| `reproducibility_wellquantified.tsv` | well-quantified counts per protease and as the union, with the metric the stage chose |
| `report.html` | self-contained stage-by-stage report with plots |
| `INTERPRETATION_AND_METHODS.md` | a methods/interpretation document generated from *this run's own* numbers, so the prose cannot drift from the results |
| `advisor_state.rds` | the complete run state; re-render reports without recomputing (`render_reports.R`) |

**Reading the deliverable.** Key columns of `final_abundance_matrix.tsv`: `recommended_channel`
(constant per protein), `source_channel` (= recommended_channel, the scale the value lives on),
`value_source` (`measured` or `imputed_cross_protease`), and for imputed cells `donor_channel`,
`donor_value`, `offset_used`. To use only measured data, filter `value_source == "measured"`; to keep
only confident calls, filter on `recommendation_confidence`; to see only the proteins where the choice
mattered, filter `consequential == TRUE`.

---

## Settings reference

Every default lives in `config_default.yml`, annotated with the reason it was chosen. The keys below
can be set on the command line as `key=value`. Anything else (QC and imputation-gate thresholds, channel
labels, seed) is changed with `config=my.yml`: a YAML file holding only the keys you want to change,
layered over the defaults.

| key | default | meaning |
|---|---|---|
| `run_include` / `run_exclude` | "" | regex on the normalized run name to keep / drop runs |
| `run_sample_regex` | "" | capture group mapping a run to a biological sample id — **set this when you have technical replicates**, so CV/ICC criteria and Stage 6 precision metrics switch on |
| `selection_mode` | `completeness` | `completeness` (manuscript method) or `multi_criteria` (weighted vote) |
| `impute` | `false` | turn cross-protease gap filling on |
| `offset_model` | `protein_median` | `protein_median` or `sample_knn` |
| `min_paired` / `n_perm` | 15 / 1000 | agreement (Stage 4) paired-sample minimum and permutation count |
| `run_day_regex` | "" | capture a preparation-day label for Stage 6 corrected inter-day CV (reproducibility cohort: `^([0-9]+)_.*$`) |
| `cv_threshold` / `completeness_min` | 20 / 0.50 | Stage 6 well-quantified thresholds (CV mode / completeness mode) |
| `probe_in_score` | `false` | let the roll-up QC probe vote (off — it correlates 0.91 with peptide count and changes no consequential calls) |
| `fasta` | "" | a FASTA unlocks sequence-coverage criteria |

`RUNNING_THE_APP.md` is the deeper operational guide (every setting, every output, the nuances that
decide whether a number should be trusted).

---

## Verify it works

```bash
Rscript tests/run_all_tests.R
```

runs the whole suite in isolated processes and summarises pass / fail / skip. Individually:

- **`tests/test_example_end_to_end.R`** — runs the whole pipeline on the bundled synthetic data (replicate
  design, both offset models, and a discovery design) and checks the tool's invariants: one source
  channel per protein, the chosen channel is the more complete one, measured values are untouched, filled
  values are exactly donor + offset, discordant proteins are never filled, Stage 6 picks the right
  metric, and every output is written. Needs no study data; this is what CI runs.
- **`tests/test_shiny_app.R`** — drives the GUI's server headlessly (`shiny::testServer`) on the same
  synthetic data: loads the four files, presses *Run advisor*, renders every table and plot, and
  downloads the HTML report. Skips if `shiny`/`DT` are not installed.
- **`tests/test_manuscript_defaults.R`** — asserts `config_default.yml` carries the manuscript
  defaults, and that a completeness-mode run's `recommended_channel` matches an **independent**
  recomputation of "the more complete channel" (not the pipeline's own code path), the one-channel-
  per-protein invariant holds, and the selection reporting is complete.
- **`tests/test_ingest.R`** — reading the pivoted Mix30 matrices yields 54 runs per search and **392
  shared proteins**, matching the trypsin∩HTA overlap reported in the manuscript.
- **`tests/test_reproducibility_stage.R`** — Stage 6 reproduces the published Fig 4A per-sample counts
  in inter-day-CV mode and falls back to completeness on data without replicates.

The last two need the study data: download the DIA-NN search results from MassIVE MSV000103457 and set
`TWO_PROTEASE_DATA` to that folder (layout in `../figure-scripts/DATA_LAYOUT.md`). Without it they
**skip** rather than fail; checks 2–4 of `test_manuscript_defaults.R` likewise skip until a
completeness-mode run exists under `results/` (or under `$ADVISOR_RESULTS`). `tests/SCENARIOS.md` has copy-paste scenarios with expected numbers. `vendor/MANIFEST.yml`
holds md5s of the reused modules so drift from the source frameworks stays detectable.

---

## Repository layout

```
config_default.yml     every default, with the reason where one was measured
install.R              one-shot dependency installer      requirements.txt   the dependency list
run_advisor.R          CLI                                 run_all_conditions.sh  the study's three runs
render_reports.R       re-render report/methods from a saved run   build_index.R  cross-run index page
example_data/          small SYNTHETIC dataset + the script that generates it
app/app.R              Shiny GUI — no statistics live here
engine/
  ingest_matrices.R    pivoted pr_matrix / pg_matrix -> long contract
  pipeline.R           the whole analysis as one function; app and CLI both call it
  impute_pg.R          Stage 5, per-protein median offset (or sample-kNN), refusals recorded
  reproducibility.R    Stage 6, picks the metric for the data
  outputs.R  plots.R  methods_doc.R   tables, figures, the self-contained report, the generated methods
  load_engine.R        single source-order entry point     utils.R  shared helpers
vendor/                snapshot of reused modules + MANIFEST.yml (md5s); recommend_rollup.R is app-modified
tests/                 the test suite + SCENARIOS.md
results/               run outputs (git-ignored; regenerate with run_advisor.R)
```

`vendor/` is a **snapshot, not a fork**: modules copied from our earlier analysis frameworks, with md5s in
`MANIFEST.yml`. The one functionally diverged file, `recommend_rollup.R`, is flagged `modified` there (it
gained `selection_mode`).

---

## Data and code availability

The raw data and DIA-NN search results are on MassIVE (**MSV000103457**). The scripts behind the
manuscript figures live next to this folder in the same repository (`../figure-scripts/`). If you use
the tool, please cite the accompanying manuscript:

> Alba MM, Bharadwaj A, Becker L, Canfield EJ, Katz JE, Van Eyk JE, Kreimer S. *Joining Two Perspectives
> of the Proteome with the Trypsin+HTA Protease Workflow.* bioRxiv (2026).
> doi:10.64898/2026.10.05.756883
