# Expected data layout

The scripts read the DIA-NN search results deposited at MassIVE (**MSV000103457**). Download them and
arrange (or symlink) them like this. `TWO_PROTEASE_DATA` points at the top folder:

```
$TWO_PROTEASE_DATA/
├── Trypsin_Library_Search/            discovery runs (324 raw files) searched with the trypsin library
│   ├── report.pr_matrix.tsv
│   ├── trypsin_library_search_report.pg_matrix.tsv     (the trypsin-library report.pg_matrix.tsv)
│   ├── report.parquet                                   (07_verification only)
│   └── lib.parquet                                      (07_verification/modification_inventory.R)
├── Krakatoa_Library_Search/           the same 324 raw files searched with the Krakatoa (HTA) library
│   ├── report.pr_matrix.tsv
│   ├── report.pg_matrix.tsv
│   ├── report.parquet
│   └── lib.parquet
├── Reproducibility_Search/            8 plasma samples x 3 days x 8 preparations = 192 runs
│   ├── Trypsin Library/
│   │   ├── report.pr_matrix.tsv  report.pg_matrix.tsv  report.parquet  lib.parquet
│   └── Krakatoa Library/
│       ├── report.pr_matrix.tsv  report.pg_matrix.tsv  report.parquet  lib.parquet
└── reference/
    └── proteinatlas.tsv               Human Protein Atlas table (07_verification/discovery_lens.R only;
                                       download proteinatlas.tsv.zip from proteinatlas.org, or set HPA_TSV)
```

Notes:
- **The trypsin discovery protein matrix is the one exception to `report.*` naming.** Our scripts read it
  as `trypsin_library_search_report.pg_matrix.tsv`. If your download names it `report.pg_matrix.tsv`,
  copy or symlink it to that name. The collaborator scripts (`02_`–`04_`) read `report.pg_matrix.tsv`
  directly, so keep both names.
- **Discovery run names** encode the preparation group: `..._DIA_Plasma_Mix_<donor>` (Mix 30 min),
  `..._Mix_<donor>_20min` (Mix 20 min), `..._Try_...` / `..._Kra_...` and the August batch
  `250807_..._Plasma_...` (Trypsin 1 / Krakatoa 1). The scripts select runs by these patterns.
- **Reproducibility run names** are `{day}_{sample}_{replicate}`, e.g. `2_07_5`.
- **FASTA.** Sequence coverage and PTM site numbering (`03_ptm_coverage`, `04_reproducibility`) use the
  UniProtKB human reviewed proteome downloaded 2025-10-28. It is bundled as
  `reference/uniprotkb_human_reviewed_2025_10_28.fasta.gz` and is the default when the scripts are run
  from the repository root; `TWO_PROTEASE_FASTA` overrides it.
  The DIA-NN searches themselves used the FASTA files deposited with the MassIVE dataset.
- **Libraries.** The combined spectral libraries are not in this repository. They are deposited with
  the MassIVE dataset.
