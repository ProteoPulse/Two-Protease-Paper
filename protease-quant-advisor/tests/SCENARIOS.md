# Test scenarios — Protease Quant Advisor

Concrete things to run, and what a correct result looks like. Scenario 0 needs nothing but this
repository; the others need the study data. Paths below assume you are **in the app
folder** (`protease-quant-advisor/`) and that `TWO_PROTEASE_DATA` points at the DIA-NN search results
from MassIVE MSV000103457, laid out as in `../figure-scripts/DATA_LAYOUT.md`. For brevity, define:

```bash
S="$TWO_PROTEASE_DATA"
TRY_PR=$S/Trypsin_Library_Search/report.pr_matrix.tsv
HTA_PR=$S/Krakatoa_Library_Search/report.pr_matrix.tsv
TRY_PG=$S/Trypsin_Library_Search/trypsin_library_search_report.pg_matrix.tsv
HTA_PG=$S/Krakatoa_Library_Search/report.pg_matrix.tsv
```

Run the automated suite at any time with `Rscript tests/run_all_tests.R`.

---

## Scenario 0 — the bundled synthetic data (no download needed)

```bash
Rscript run_advisor.R \
  pr_a=example_data/trypsin_report.pr_matrix.tsv pr_b=example_data/hta_report.pr_matrix.tsv \
  pg_a=example_data/trypsin_report.pg_matrix.tsv pg_b=example_data/hta_report.pg_matrix.tsv \
  run_sample_regex='^[0-9]+_([0-9]+)_[0-9]+$' run_day_regex='^([0-9]+)_.*$' impute=true tag=example
```

**Expect** (~30 s): 16 samples, 129 shared proteins (112 trypsin / 17 HTA), 8 discordant (the planted
ones), 101 cells imputed in 27 proteins, Stage 6 in inter-day-CV mode. `tests/test_example_end_to_end.R`
checks all of this automatically.

## Scenario A — Mix30 discovery, the manuscript default (completeness)

```bash
Rscript run_advisor.R pr_a=$TRY_PR pr_b=$HTA_PR pg_a=$TRY_PG pg_b=$HTA_PG \
  run_include=Plasma_Mix run_exclude=20min impute=true \
  selection_mode=completeness tag=mix30_completeness out=results/mix30_completeness
```

**Expect:** 392 shared proteins adjudicated (matches the trypsin∩HTA overlap in the manuscript), split
**322 trypsin / 70 HTA**, 152 consequential; concordance 82 concordant / 152 discordant / 137 not
assessable; imputation fills **31 cells in 3 proteins** (2,926 gaps refused, each with its reason); a final
matrix of 2,393 proteins × 54 samples (71,561 rows); every protein carrying exactly one `source_channel`
(`INTERPRETATION_AND_METHODS.md` shows `proteins carrying more than one = 0`); the run note stating
"Channel selection = COMPLETENESS". Open `results/mix30_completeness/report.html`.

**Check:** in `tables/recommendation.tsv`, `recommended_channel` should equal the channel with the
higher `protein_detection_rate` in `tables/precursor_evidence.tsv` for every protein. This is exactly
what `tests/test_manuscript_defaults.R` verifies automatically.

## Scenario B — same data, the weighted-vote generalization (multi_criteria)

```bash
Rscript run_advisor.R pr_a=$TRY_PR pr_b=$HTA_PR pg_a=$TRY_PG pg_b=$HTA_PG \
  run_include=Plasma_Mix run_exclude=20min impute=true \
  selection_mode=multi_criteria tag=mix30_multicriteria out=results/mix30_multicriteria
```

**Expect:** the same 392 shared proteins but a different split (the archived multi-criteria Mix30 run
is 299 trypsin / 93 HTA). Compare the two runs:

```bash
Rscript -e 'library(data.table);
  a=fread("results/mix30_completeness/tables/recommendation.tsv");
  b=fread("results/mix30_multicriteria/tables/recommendation.tsv");
  m=merge(a[,.(pg_canonical,completeness=recommended_channel)],
          b[,.(pg_canonical,vote=recommended_channel)],by="pg_canonical");
  cat("disagree:",m[completeness!=vote,.N],"of",nrow(m),"\n")'
```

**Interpret:** the disagreement count is the "small set where two channels are comparably complete but
differ in precision" the manuscript describes. In the completeness run, the same proteins are flagged
in-place by `selection_disagrees_with_vote == TRUE`.

## Scenario C — reproducibility cohort, unlocking Stage 6 (inter-day CV)

```bash
RS=$S/Reproducibility_Search
Rscript run_advisor.R \
  pr_a="$RS/Trypsin Library/report.pr_matrix.tsv" pr_b="$RS/Krakatoa Library/report.pr_matrix.tsv" \
  pg_a="$RS/Trypsin Library/report.pg_matrix.tsv"  pg_b="$RS/Krakatoa Library/report.pg_matrix.tsv" \
  run_sample_regex='^[0-9]+_([0-9]+)_[0-9]+$' run_day_regex='^([0-9]+)_.*$' \
  impute=true tag=repro_v2 out=results/repro_v2
```

**Expect:** 8 biological samples detected (the `run_sample_regex` collapses replicate injections);
`min_paired = 15` is unreachable on 8 samples, so the run prints a **clamp note** and proceeds — this
is correct, not a bug. Stage 6 runs in **corrected inter-day CV** mode (replicates + a day map are
present) and reports well-quantified counts per protease and as the union in
`tables/reproducibility_wellquantified.tsv`: **1,043 trypsin / 894 Krakatoa / 1,510 union** protein ×
sample pairs under the corrected metric, and 1,000 / 731 / 1,322 under the published day-mean CV (whose
per-sample counts reproduce Fig 4A exactly). 94 shared proteins, 68 trypsin / 26 HTA. With replicates present, `technical_cv` and `replicate_icc` become
live evidence criteria.

## Scenario D — imputation on vs off (read the refusals)

Run Scenario A once with `impute=true` and once with `impute=false`. **Expect:** the recommendation and
the measured values are identical; only `final_abundance_matrix.tsv` differs by the handful of imputed
cells (Mix30 fills very few). Open `tables/imputation_refused.tsv`: most gaps are **refused** on
purpose — a donor whose values do not track the recipient is not a surrogate. `imputation_ledger.tsv`
lists every cell that *was* filled with its donor, offset and the gate values it passed. A gap left
empty is a decision, not a shortfall.

## Scenario E — the "no shared samples" guardrail

Feed run filters that keep disjoint runs from the two searches (e.g. `run_include` matching only
trypsin-side run names). **Expect:** the pipeline stops immediately with
*"the two matrices share no sample columns after run filtering"* rather than producing an empty or
misleading result. This confirms the tool refuses to compare searches that never saw the same
injections.

## Scenario F — the GUI end to end

```bash
Rscript -e 'shiny::runApp("app", launch.browser = TRUE)'
```

Load the four matrices, set `keep runs matching = Plasma_Mix`, `drop runs matching = 20min`, leave
**Channel selection = completeness (manuscript method)**, press **Run advisor**. **Expect:** the same
392-shared headline as Scenario A, a tab per stage with its decision plot, and identical numbers to the
CLI (the GUI and CLI call the same `run_advisor()`), downloadable as a zip of tables and a single HTML
report. Switch Channel selection to **multi-criteria** and re-run to see the split move, matching
Scenario B.
