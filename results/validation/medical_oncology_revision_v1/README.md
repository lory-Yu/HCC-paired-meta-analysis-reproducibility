# Medical Oncology revision audit outputs

This directory contains the version-locked identifier and denominator audits added for release `v0.2.0-submission`.

## Primary analysis set

The primary external-validation analysis uses only unambiguously harmonized current gene symbols:

- 173 prespecified discovery-robust candidates;
- 167 estimable in ICGC-LIRI-JP after identifier harmonization;
- 163 jointly estimable in TCGA-LIHC and ICGC-LIRI-JP;
- 157 passing candidate-family FDR in both cohorts;
- 155 passing genome-wide FDR in both cohorts.

The four ambiguous historical symbols are documented in `historical_symbol_ambiguous_four_decision.tsv` and excluded from the primary denominator. The 34-symbol and current-symbol-only variants are explicitly labelled sensitivity analyses rather than interchangeable primary estimates.

## Secondary single-cell sensitivity

`singlecell_donor_unpaired_countsum_pseudobulk_results.tsv` and its companion files report an unpaired donor-level count-sum sensitivity analysis of public GSE202642 data. One hundred candidate–cell-type contrasts were testable and none passed FDR < 0.05. This negative result is retained because it constrains interpretation; it is not a cell-level independent-replicate analysis and is not treated as mechanistic validation.

## File guide

- `external_validation_primary_unambiguous_harmonized_gene_level.tsv`: primary gene-level TCGA/ICGC validation table.
- `external_validation_denominator_decision_summary.tsv`: locked counts and analysis-set labels.
- `historical_symbol_34_ambiguity_audit.tsv`: audit of all historical-symbol cases.
- `historical_symbol_ambiguous_four_decision.tsv`: explicit decisions for the four unresolved mappings.
- `genomewide_FDR_sensitivity_gene_level.tsv`: genome-wide multiplicity sensitivity.
- `genes_pass_candidate_fail_genomewide_FDR.tsv`: genes whose classification changes under the broader family.
- `singlecell_donor_design_decision.tsv`: prespecified interpretation boundary for the donor-level sensitivity.
- `singlecell_donor_tissue_celltype_cellcounts.tsv`: donor/tissue/cell-type support counts.
- `singlecell_donor_unpaired_countsum_pseudobulk_results.tsv`: valid-JSON donor-level results.
- `singlecell_donor_unpaired_countsum_summary.tsv`: compact test-count and FDR summary.

All files are derived from public datasets. No restricted clinical records or third-party raw expression files are redistributed.
