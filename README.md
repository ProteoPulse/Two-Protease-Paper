# Two-Protease-Paper

Code for **"Joining Two Perspectives of the Proteome with the Trypsin+HTA Protease Workflow"**
(Alba MM, Bharadwaj A, Becker L, Canfield EJ, Katz JE, Van Eyk JE, Kreimer S. *bioRxiv* 2026,
[doi:10.64898/2026.10.05.756883](https://doi.org/10.64898/2026.10.05.756883)).

The repository has two parts:

| folder | what it is |
|---|---|
| [`protease-quant-advisor/`](protease-quant-advisor/) | **Protease Quant Advisor**: an R tool (command line + Shiny GUI) that joins the DIA-NN outputs of a trypsin-library search and an HTA (Krakatoa)-library search of the same injections. For each protein seen by both, it picks the protease to quantify it with, checks that the two proteases agree, and fills gaps in the chosen protease from the other one (cross-protease inference). Ships with a small synthetic dataset and tests. |
| [`figure-scripts/`](figure-scripts/) | The R scripts behind the manuscript's figures and numbers, organised by analysis step. |

**Data.** Raw files, DIA-NN search results and spectral libraries are on MassIVE,
**[MSV000103457](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000103457)**. Nothing in
this repository is real patient data; `protease-quant-advisor/example_data/` is simulated.

---

## Quick start: the Quant Advisor

Requires R ≥ 4.5.

```bash
cd protease-quant-advisor
Rscript install.R                       # data.table, ggplot2, yaml, base64enc (+ shiny, DT for the GUI)
Rscript tests/run_all_tests.R           # runs the pipeline end to end on the bundled synthetic data

# your own data: two DIA-NN searches of the same runs
Rscript run_advisor.R pr_a=trypsin/report.pr_matrix.tsv pr_b=hta/report.pr_matrix.tsv \
                      pg_a=trypsin/report.pg_matrix.tsv pg_b=hta/report.pg_matrix.tsv \
                      impute=true tag=my_run

# or the GUI
Rscript -e 'shiny::runApp("app", launch.browser = TRUE)'
```

[`protease-quant-advisor/README.md`](protease-quant-advisor/README.md) covers the method, inputs, outputs
and every setting. [`RUNNING_THE_APP.md`](protease-quant-advisor/RUNNING_THE_APP.md) is the
operational guide.

---

## Reproducing the figures

1. Download the DIA-NN search results from MassIVE MSV000103457 and lay them out as in
   [`figure-scripts/DATA_LAYOUT.md`](figure-scripts/DATA_LAYOUT.md).
2. Point the scripts at them:

   ```bash
   export TWO_PROTEASE_DATA=/path/to/search_results    # required
   export TWO_PROTEASE_OUT=$PWD/output                 # optional; default <repo>/output
   export TWO_PROTEASE_FASTA=/path/to/uniprotkb_Human_AND_reviewed_true_AND_m_2025_10_28.fasta  # 03_, 04_
   ```
3. Run the scripts with `Rscript` from any directory, e.g.
   `Rscript figure-scripts/05_cross_protease/completeness_experiment.R`. Outputs go to
   `$TWO_PROTEASE_OUT/<step>/{figures,tables}`.

### Figure → script

| manuscript | content | script |
|---|---|---|
| Methods: libraries | merge DDA + 1 m/z DIA libraries per protease; PEAKS → DIA-NN | `01_libraries/Combine_Libraries.R`, `01_libraries/PEAKS_to_DIA_NN_Library_Converter.R` |
| Fig 2A–B | identification matrices, completeness, peptides per protein | `02_identification/Identification_Analysis.R` |
| Fig 1B, Fig S5 | search overlap counts; mixture vs trypsin R² | `02_identification/Identification_Analysis.R` (overlap heatmaps, R² scatters) |
| Fig 3A | PTM sites per sample and PTM-site overlap | `03_ptm_coverage/PTM_Analysis.R` |
| Fig 3B | sequence coverage and PTM maps | `03_ptm_coverage/sequence_coverage_comparison.Rmd` |
| Fig 4A (measured) | inter-day CV violins, well-quantified counts per sample | `04_reproducibility/Reproducibility_Data_Processing.R` |
| Fig 4A (with inference) | the same panel with cross-protease inference on | `05_cross_protease/make_panel_figures.R` (+ `make_manuscript_figures.R` for the standalone variants; shared code in `cv_common.R`) |
| Fig 4B | CEP95 and AGRN profiles: cross-protease inference vs mean and minimum imputation | `05_cross_protease/profile_mean_imputation.R` |
| Fig 4C | data completeness, measured vs recovered by inference (hatched) | `05_cross_protease/completeness_experiment.R` → `completeness_experiment_shaded.R` |
| Fig S7 | masking test: 30% of values hidden; cross-protease vs minimum / median / mean / kNN | `06_supplementary/masking_validation.R` |
| text | well-quantified gain per protease, discovery completeness, offset-model check | `05_cross_protease/panel_c_imputation.R`, `06_supplementary/{discovery_imputation,offset_model_check,imputation_gain_and_validity}.R` |
| verification | independent re-derivation of the Mix 30 vs Trypsin 2 completeness/coverage numbers; the PTM inventory submitted to MassIVE | `07_verification/*` |

Figures 1A, S1–S4 are schematics or instrument diagrams and have no code.

### Run order

Most scripts are independent. Where one script reads another's tables, run them in this order:

- `05_cross_protease/completeness_experiment.R` → `completeness_experiment_shaded.R`
- `06_supplementary/discovery_imputation.R` and `masking_validation.R` → `imputation_gain_and_validity.R`
- `05_cross_protease/profile_mean_imputation.R` optionally cross-checks its channel choice against a
  full Quant Advisor run on Mix 30. Run `protease-quant-advisor/run_all_conditions.sh` first, or point
  `ADVISOR_STATE` at an `advisor_state.rds`. It skips the check if neither is available.
- `07_verification/discovery_lens.R` needs the Human Protein Atlas table (`HPA_TSV`; see `DATA_LAYOUT.md`).

### Collaborator scripts (`01_`–`04_`)

These are released as used for the manuscript. Only file paths, line endings and a provenance header
were changed. They were written for interactive use on Windows:

- Folder and file pickers are bypassed by the environment variables above. Text prompts take their
  default answers (`ASSUME_DEFAULTS <- TRUE`).
- `01_libraries/*` are function libraries. Call `combine_diann_libraries(files = ..., output_file = ...)`
  and `convert_peaks_to_diann(...)` with explicit paths, or run them interactively.
- `PTM_Analysis.R` fetches UniMod names from unimod.org when its local table lacks one.
- `sequence_coverage_comparison.Rmd` loads several packages it does not use (including a mouse
  annotation package). Install them, or delete those `library()` lines, before knitting.

### Packages

- Our scripts (`05_`–`07_`): data.table, ggplot2, patchwork, scales, ragg, duckdb (07 only), grid.
- Collaborator scripts additionally need: cowplot, viridisLite, ggnewscale, Biostrings (Bioconductor),
  stringr, xml2, rmarkdown/knitr, plus the packages listed at the top of the Rmd.

---

## Citation

Please cite the preprint (see [`CITATION.cff`](CITATION.cff)). Code is released under the
[MIT License](LICENSE).
