#!/usr/bin/env bash
set -euo pipefail

project_root=${1:?project root required}
config="$project_root/workflow/config/downloads.tsv"
download_root="$project_root/data/downloads"
provenance_root="$project_root/results/provenance"
mkdir -p "$download_root" "$provenance_root"

download_one() {
  local cohort=$1
  local url=$2
  local local_name=$3
  destination_dir="$download_root/$cohort"
  destination="$destination_dir/$local_name"
  mkdir -p "$destination_dir"
  if [[ -s "$destination" ]]; then
    echo "Present: $cohort/$local_name"
    return 0
  fi
  echo "Downloading: $cohort/$local_name"
  curl --location --fail --silent --show-error --retry 5 --retry-delay 5 --continue-at - \
    --output "$destination.part" "$url"
  mv "$destination.part" "$destination"
  echo "Completed: $cohort/$local_name"
}

pids=()
while IFS=$'\t' read -r cohort role url local_name required; do
  download_one "$cohort" "$url" "$local_name" &
  pids+=("$!")
  if (( ${#pids[@]} >= 4 )); then
    wait "${pids[0]}"
    pids=("${pids[@]:1}")
  fi
done < <(tail -n +2 "$config")
for pid in "${pids[@]}"; do
  wait "$pid"
done

manifest="$provenance_root/download_manifest.tsv"
printf 'cohort\trole\turl\tlocal_path\tsize_bytes\tsha256\taccessed_utc\n' > "$manifest"
accessed_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
tail -n +2 "$config" | while IFS=$'\t' read -r cohort role url local_name required; do
  local_path="$download_root/$cohort/$local_name"
  if [[ ! -s "$local_path" ]]; then
    if [[ "$required" == "yes" ]]; then
      echo "Required download missing: $local_path" >&2
      exit 1
    fi
    continue
  fi
  size_bytes=$(stat -c '%s' "$local_path")
  sha256=$(sha256sum "$local_path" | awk '{print $1}')
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$cohort" "$role" "$url" "$local_path" "$size_bytes" "$sha256" "$accessed_utc" >> "$manifest"
done

echo "Wrote $manifest"
