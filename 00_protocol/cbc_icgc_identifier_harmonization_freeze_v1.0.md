# ICGC-LIRI-JP identifier harmonization sensitivity analysis v1.0

Date: 2026-10-06

## Purpose

This is a post hoc technical sensitivity analysis prompted by an audit of ICGC non-estimability. It tests whether discovery candidates recorded under current HGNC symbols can be linked to historical symbols in the frozen ICGC-LIRI-JP count matrix. It does not change the six-cohort discovery analysis, either candidate lock, the paired ICGC model, the expression filter, or any significance threshold.

## Locked inputs

- Discovery robust set: `candidate_set_v2_high_confidence.tsv`, 173 genes, SHA-256 `55110e282bf9d26619de70d6e8c1e5972d050766b034a6fa06009021ae9fd3cd`.
- Intermediate set: `candidate_set_v2_intermediate.tsv`, 119 genes, SHA-256 `23812a47c8df31301aca8a9061ecd7da56dd41afa5f214a3fc6e4a585494aa59`.
- Frozen ICGC full paired count matrix: `results/validation/rnaseq_normalized_v2/icgc_liri_jp_paired_full_counts.tsv.gz`.
- Existing complete-transcriptome ICGC model: `results/validation/cbc_external_validation_v1/icgc_complete_filtered_transcriptome.tsv.gz`.
- Frozen historical-symbol map: `00_protocol/cbc_icgc_historical_symbol_map_v1.tsv`.

## Mapping rule

1. Match the current candidate symbol directly when it exists in the frozen matrix.
2. Otherwise use only a historical symbol listed in the frozen map when the candidate Entrez ID and the `org.Hs.eg.db` alias relation agree and the historical symbol exists in the matrix.
3. Each candidate may map to at most one matrix row. Ambiguous or absent mappings remain not estimable.
4. Mapping is determined without inspecting expression effects or P values.

## Statistical rule

The paired limma-voom model is not re-fit. The harmonized candidate table reuses the already fitted ICGC complete-transcriptome row corresponding to the locked historical symbol. Genome-wide FDR therefore remains the BH value from the same 19,353-gene filtered universe. Candidate-family BH-FDR is recalculated within each locked tier after harmonization. TCGA results are unchanged.

## Interpretation

The symbol-only analysis remains available as the original external-validation result. The harmonized result is explicitly labelled post hoc and technical. It may show that non-estimability was caused by identifier versioning, but it must not be presented as an independent cohort or as biological validation beyond the underlying ICGC data.

