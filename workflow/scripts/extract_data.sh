#!/usr/bin/env bash
set -euo pipefail

project_root=${1:?project root required}
download_root="$project_root/data/downloads"
work_root="$project_root/data/work"
manifest="$project_root/00_protocol/sample_manifest_v1.0.tsv"
mkdir -p "$work_root"

for cohort in GSE121248 GSE45114 GSE57555 GSE57957 GSE76427 GSE84402; do
  archive="$download_root/$cohort/${cohort}_RAW.tar"
  cohort_dir="$work_root/$cohort/raw"
  marker="$cohort_dir/.archive_extracted"
  mkdir -p "$cohort_dir"
  if [[ ! -f "$marker" ]]; then
    tar -xf "$archive" -C "$cohort_dir"
    touch "$marker"
  fi
done

# Decompress only analysis-eligible sample files. Original archives and compressed
# members remain untouched.
awk -F '\t' 'NR>1 && $13=="include" {print $1"\t"$2}' "$manifest" | while IFS=$'\t' read -r cohort gsm; do
  while IFS= read -r -d '' compressed; do
    uncompressed=${compressed%.gz}
    if [[ ! -s "$uncompressed" ]]; then
      gzip -cd "$compressed" > "$uncompressed.part"
      mv "$uncompressed.part" "$uncompressed"
    fi
  done < <(find "$work_root/$cohort/raw" -type f -name "${gsm}*.gz" -print0)
done

echo "Raw archives extracted under $work_root"
