# HCC/CBC public analysis specification

## Scope

This public specification covers only the paired multi-cohort hepatocellular carcinoma transcriptomic analysis described in the accompanying manuscript. It excludes pharmacology, herbal-material, treatment, prevention and causal-mechanism claims.

## Frozen discovery inputs

- Cohorts: GSE121248, GSE45114, GSE57555, GSE57957, GSE76427 and GSE84402.
- Primary contrast: tumour minus patient-matched non-tumour tissue.
- Intent-to-analyse samples: 439; quality-controlled samples: 438.
- Primary complete pairs: 167.
- The sole sample-level exclusion was GSE57957/GSM1398656 under the prespecified two-metric QC rule; its pair was excluded from the primary paired model.
- Platform-specific preprocessing, sample status and pairing are recorded in `cohort_manifest_v1.0.tsv`, `sample_manifest_v1.0.tsv` and the workflow scripts.

## Discovery inference

- Each cohort was analysed separately using a patient-adjusted limma model.
- Gene effects represented in at least four cohorts were combined with a REML random-effects model and Hartung–Knapp inference.
- Benjamini–Hochberg adjustment used the complete successfully fitted discovery universe of 17,900 genes.
- A discovery hit required FDR <0.05, absolute pooled log2 fold change >=0.5849625 and directional agreement in at least 75% of contributing cohorts.
- The discovery analysis yielded 446 recurrent changes. Prespecified robustness/downgrade rules produced 173 discovery-robust, 119 intermediate and 154 exploratory genes.

## Robustness analyses

- Leave-one-cohort-out analyses were recomputed over the complete estimable universe in every round.
- A fixed-effect model was retained as a comparison.
- A blocked analysis used all 438 quality-controlled samples.
- Affymetrix-only, Illumina-only, adjacent-non-tumour-only and declared cohort-exclusion analyses were sensitivity analyses and did not enlarge the candidate set.

## External paired RNA-seq validation

- TCGA-LIHC used 50 paired patients.
- ICGC-LIRI-JP Release 28 used 199 paired donors.
- Candidate-family FDR was calculated within each locked evidence tier among estimable genes and was the primary locked-list replication criterion.
- Complete-transcriptome FDR was calculated across all finite fitted P values in each cohort-specific filtered transcriptome and was labelled post hoc sensitivity analysis.
- Identifier harmonization in ICGC used a frozen Entrez/alias mapping without inspecting recovered effects; the original symbol-only results were retained.
- Among 173 discovery-robust genes, 167 were jointly estimable, 166 retained direction in both external cohorts, 161 met the dual candidate-family criterion and 159 met the dual complete-transcriptome sensitivity criterion.

## Single-cell endpoint

GSE202642 was used only for descriptive cellular localization. It contained seven tumour and four adjacent-tissue samples after sample accounting. Cells were not treated as independent biological replicates, and the exploratory tumour-versus-adjacent comparison was not used as differential validation.

## Evidence boundaries

The public outputs support recurrent tumour-associated tissue-level expression changes and a strictly defined external replication subset. They do not establish diagnostic performance, clinical prognosis, treatment benefit, recurrence prevention, protein-level validation or causal mechanism.
