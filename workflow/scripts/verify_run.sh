#!/usr/bin/env bash
set -euo pipefail

project_root=${1:?project root required}
required=(
  "$project_root/results/provenance/download_manifest.tsv"
  "$project_root/results/provenance/raw_file_manifest.tsv"
  "$project_root/results/provenance/conda-packages.json"
  "$project_root/results/provenance/R-sessionInfo.txt"
  "$project_root/results/qc/input_validation_summary.tsv"
  "$project_root/results/qc/qc_summary.tsv"
  "$project_root/results/qc/qc_decisions.tsv"
  "$project_root/results/qc/qc_exclusions.tsv"
  "$project_root/results/qc/qc_review.md"
  "$project_root/workflow/environment-linux-64.lock"
  "$project_root/results/provenance/output_checksums.sha256"
)

for path in "${required[@]}"; do
  [[ -s "$path" ]] || { echo "Missing required artifact: $path" >&2; exit 1; }
done

awk -F '\t' 'NR>1 && $7!="pass" {bad++} END {exit bad>0}' \
  "$project_root/results/qc/input_validation_summary.tsv"

awk -F '\t' 'NR>1 {n++; status=$2; gsub(/"/, "", status); if (status!="completed") bad++} END {exit n!=6 || bad>0}' \
  "$project_root/results/qc/qc_summary.tsv"

sha256sum -c <(awk 'NR>1 {print $6"  "$4}' \
  "$project_root/results/provenance/download_manifest.tsv") >/dev/null

for archive in "$project_root"/data/downloads/GSE*/GSE*_RAW.tar; do
  tar -tf "$archive" >/dev/null
done

echo "Workflow verification passed."
