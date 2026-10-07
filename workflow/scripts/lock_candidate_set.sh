#!/usr/bin/env bash
set -euo pipefail

project_root=${1:?project root required}
source_table="$project_root/results/analysis/meta/meta_primary_hits.tsv"
lock_dir="$project_root/results/analysis/candidate_lock"
locked_table="$lock_dir/candidate_set_v1.0.tsv"
metadata="$lock_dir/candidate_set_v1.0.txt"

[[ -s "$source_table" ]] || { echo "Missing discovery candidate source: $source_table" >&2; exit 1; }
mkdir -p "$lock_dir"

if [[ -e "$locked_table" ]]; then
  cmp -s "$source_table" "$locked_table" || {
    echo "Candidate lock differs from current discovery result; refusing to overwrite v1.0" >&2
    exit 1
  }
else
  cp "$source_table" "$locked_table"
fi

source_sha=$(sha256sum "$source_table" | awk '{print $1}')
candidate_count=$(awk 'END {print NR-1}' "$locked_table")
{
  echo "candidate_lock_version=v1.0"
  echo "source=results/analysis/meta/meta_primary_hits.tsv"
  echo "selection_basis=frozen discovery meta-analysis rule only"
  echo "candidate_count=$candidate_count"
  echo "sha256=$source_sha"
  echo "external_validation_must_not_modify_discovery_set=true"
} > "$metadata"

echo "Candidate set v1.0 locked: $candidate_count genes ($source_sha)"
