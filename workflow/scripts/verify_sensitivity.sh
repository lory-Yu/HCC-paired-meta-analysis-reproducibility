#!/usr/bin/env bash
set -euo pipefail

project_root=${1:?project root required}
root="$project_root/results/analysis/sensitivity"
required=(
  "$project_root/results/analysis/candidate_lock/candidate_set_v1.0.tsv"
  "$project_root/results/analysis/candidate_lock/candidate_set_v1.0.txt"
  "$root/fixed_effect_all_genes.tsv.gz"
  "$root/fixed_effect_candidate_comparison.tsv"
  "$root/leave_one_cohort_out/leave_one_cohort_out_all_genes.tsv.gz"
  "$root/leave_one_cohort_out/leave_one_cohort_out_failures.tsv"
  "$root/leave_one_cohort_out/candidate_leave_one_out_summary.tsv"
  "$root/all_samples/all_samples_cohort_summary.tsv"
  "$root/all_samples/all_samples_meta_all_genes.tsv.gz"
  "$root/all_samples/all_samples_candidate_comparison.tsv"
  "$root/figures/candidate_effect_concordance.png"
  "$root/figures/leave_one_out_candidate_stability.png"
  "$root/figures/leave_one_out_pass_distribution.png"
  "$root/sensitivity_summary.tsv"
  "$root/run_summary.md"
  "$root/prespecified_subsets/subset_sensitivity_summary.tsv"
  "$root/prespecified_subsets/cohort_exclusion_summary.tsv"
  "$root/prespecified_subsets/subset_sensitivity_review.md"
  "$root/prespecified_subsets/affymetrix_only_all_genes.tsv.gz"
  "$root/prespecified_subsets/illumina_only_all_genes.tsv.gz"
  "$root/prespecified_subsets/adjacent_non_tumour_only_all_genes.tsv.gz"
  "$root/candidate_downgrade_audit.tsv"
  "$root/sensitivity_review_zh.md"
  "$project_root/results/analysis/candidate_lock/candidate_set_v2_high_confidence.tsv"
  "$project_root/results/analysis/candidate_lock/candidate_set_v2_intermediate.tsv"
  "$project_root/results/analysis/candidate_lock/candidate_set_v2_metadata.txt"
)
for path in "${required[@]}"; do
  [[ -s "$path" ]] || { echo "Missing sensitivity artifact: $path" >&2; exit 1; }
done

locked=$(awk 'END {print NR-1}' "$project_root/results/analysis/candidate_lock/candidate_set_v1.0.tsv")
[[ "$locked" -eq 446 ]] || { echo "Unexpected locked candidate count: $locked" >&2; exit 1; }

audit_n=$(awk 'END {print NR-1}' "$root/candidate_downgrade_audit.tsv")
[[ "$audit_n" -eq 446 ]] || { echo "Unexpected downgrade audit count: $audit_n" >&2; exit 1; }

high_n=$(awk 'END {print NR-1}' "$project_root/results/analysis/candidate_lock/candidate_set_v2_high_confidence.tsv")
int_n=$(awk 'END {print NR-1}' "$project_root/results/analysis/candidate_lock/candidate_set_v2_intermediate.tsv")
expl_n=$(awk 'END {print NR-1}' "$root/exploratory/candidate_set_v2_exploratory_only.tsv")
[[ $((high_n + int_n + expl_n)) -eq 446 ]] || {
  echo "Confidence strata do not sum to 446: high=$high_n intermediate=$int_n exploratory=$expl_n" >&2
  exit 1
}

samples=$(awk -F '\t' '$1=="all_samples_retained" {print $2}' "$root/sensitivity_summary.tsv")
[[ "$samples" -eq 438 ]] || { echo "Unexpected retained all-sample count: $samples" >&2; exit 1; }

awk -F '\t' 'NR>1 {n++; samples+=$2} END {exit n!=6 || samples!=438}' \
  "$root/all_samples/all_samples_cohort_summary.tsv"
gzip -t "$root/fixed_effect_all_genes.tsv.gz"
gzip -t "$root/leave_one_cohort_out/leave_one_cohort_out_all_genes.tsv.gz"
gzip -t "$root/all_samples/all_samples_meta_all_genes.tsv.gz"
for scenario in affymetrix_only illumina_only adjacent_non_tumour_only; do
  gzip -t "$root/prespecified_subsets/${scenario}_all_genes.tsv.gz"
  awk -F '\t' -v s="$scenario" '
    NR==1 {next}
    {scenario_name=$1; gsub(/^"|"$/, "", scenario_name)}
    scenario_name==s {
      found=1
      if ($5 != 446 || $6 != 446 || $9 != "FALSE") exit 1
    }
    END {if (!found) exit 1}
  ' "$root/prespecified_subsets/subset_sensitivity_summary.tsv"
done
echo "Sensitivity analysis verification passed."
