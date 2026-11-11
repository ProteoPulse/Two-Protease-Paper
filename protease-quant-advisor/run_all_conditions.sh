#!/usr/bin/env bash
# run_all_conditions.sh — the three datasets in this study, through the same tool with the same
# settings. Run sequentially, not in parallel: Stage 2 holds the whole precursor table in memory.
set -uo pipefail
cd "$(dirname "$0")"
# Point TWO_PROTEASE_DATA at the DIA-NN search results (MassIVE MSV000103457; layout in
# ../figure-scripts/DATA_LAYOUT.md).
S="${TWO_PROTEASE_DATA:?set TWO_PROTEASE_DATA to the folder holding the DIA-NN search results}"
mkdir -p results
DISC_TRY="$S/Trypsin_Library_Search"; DISC_HTA="$S/Krakatoa_Library_Search"
REP_TRY="$S/Reproducibility_Search/Trypsin Library"; REP_HTA="$S/Reproducibility_Search/Krakatoa Library"

run () {  # tag, pr_a, pr_b, pg_a, pg_b, then extra key=value args
  local tag="$1"; shift
  echo "════════ $tag ════════  $(date '+%H:%M:%S')"
  Rscript run_advisor.R pr_a="$1" pr_b="$2" pg_a="$3" pg_b="$4" \
    "${@:5}" impute=true tag="$tag" out="results/$tag" 2>&1 | tee "results/${tag}_console.log" | tail -3
  echo "──────── $tag exit=$? ────────"
}

# Discovery, 30-minute gradient: Plasma_Mix runs that are not the 20-minute set.
run mix30 "$DISC_TRY/report.pr_matrix.tsv" "$DISC_HTA/report.pr_matrix.tsv" \
          "$DISC_TRY/trypsin_library_search_report.pg_matrix.tsv" "$DISC_HTA/report.pg_matrix.tsv" \
          run_include="Plasma_Mix" run_exclude="20min"

# Discovery, 20-minute gradient.
run mix20 "$DISC_TRY/report.pr_matrix.tsv" "$DISC_HTA/report.pr_matrix.tsv" \
          "$DISC_TRY/trypsin_library_search_report.pg_matrix.tsv" "$DISC_HTA/report.pg_matrix.tsv" \
          run_include="Plasma_Mix.*20min"

# Reproducibility: 192 runs named {day}_{sample}_{replicate}; field 2 is the biological sample.
# This is the only dataset with replicates, so technical CV and replicate ICC become live criteria.
run repro "$REP_TRY/report.pr_matrix.tsv" "$REP_HTA/report.pr_matrix.tsv" \
          "$REP_TRY/report.pg_matrix.tsv" "$REP_HTA/report.pg_matrix.tsv" \
          run_sample_regex='^[0-9]+_([0-9]+)_[0-9]+$'

echo "ALL DONE $(date '+%H:%M:%S')"
