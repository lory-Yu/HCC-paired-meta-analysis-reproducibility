#!/usr/bin/env bash
set -euo pipefail

project_root=${1:?project root required}
analysis_root="$project_root/results/analysis"
required=(
  "$analysis_root/annotation/annotation_summary.tsv"
  "$analysis_root/differential/cohort_differential_summary.tsv"
  "$analysis_root/meta/meta_input_long.tsv.gz"
  "$analysis_root/meta/meta_all_genes.tsv.gz"
  "$analysis_root/meta/meta_primary_hits.tsv"
  "$analysis_root/meta/meta_top20_primary_hits.tsv"
  "$analysis_root/meta/meta_analysis_failures.tsv"
  "$analysis_root/meta/meta_analysis_warnings.tsv"
  "$analysis_root/meta/meta_analysis_retries.tsv"
  "$analysis_root/figures/meta_volcano.png"
  "$analysis_root/figures/top12_forest_plot.png"
  "$analysis_root/analysis_summary.tsv"
  "$analysis_root/run_summary.md"
)
for path in "${required[@]}"; do
  [[ -s "$path" ]] || { echo "Missing analysis artifact: $path" >&2; exit 1; }
done

awk -F '\t' 'NR>1 {n++; pairs+=$3} END {exit n!=6 || pairs!=167}' \
  "$analysis_root/differential/cohort_differential_summary.tsv"

for cohort in GSE121248 GSE45114 GSE57555 GSE57957 GSE76427 GSE84402; do
  [[ -s "$analysis_root/differential/${cohort}_paired_limma.tsv.gz" ]] || exit 1
  gzip -t "$analysis_root/differential/${cohort}_paired_limma.tsv.gz"
done
gzip -t "$analysis_root/meta/meta_all_genes.tsv.gz"
awk -F '\t' '$1=="genes_meta_requested" {requested=$2} $1=="genes_meta_analyzed" {done=$2} $1=="genes_meta_failed" {failed=$2} END {exit requested!=(done+failed)}' \
  "$analysis_root/analysis_summary.tsv"
echo "Differential/meta-analysis verification passed."
