# Release notes: v0.2.0-submission

Release date: 2026-11-22

This minor release aligns the public research compendium with the Medical Oncology submission manuscript and makes the revised external-validation denominator fully auditable.

## Added

- a primary, unambiguously harmonized 173-gene external-validation table;
- explicit 173 → 167 → 163 denominator accounting;
- separate candidate-family and genome-wide FDR classifications (157 and 155 dual-cohort passes, respectively);
- a 34-symbol historical-identifier audit and explicit exclusion decisions for four ambiguous symbols;
- current-symbol-only and inclusive historical-symbol sensitivity outputs;
- a donor-level, unpaired count-sum single-cell sensitivity analysis with its negative-result summary;
- revised Fig. 1, Fig. 2 and Fig. 5 assets and figure source data;
- portable rebuild scripts, release-specific verification and privacy auditing.

## Interpretation boundary

The donor-level single-cell analysis is a secondary sensitivity analysis and does not establish cell-intrinsic regulation or mechanism. The compendium contains computational evidence only; it contains no wet-laboratory validation and makes no treatment-efficacy claim.

## Deliberately excluded

- journal correspondence, decision letters, cover letters and submission forms;
- manuscript submission files in Markdown, Word or PDF form containing author contact details or journal-specific text;
- private access instructions or internal review reports;
- local absolute paths, temporary files, credentials and access tokens;
- third-party raw GEO, TCGA, ICGC or single-cell data.

These exclusions keep the public archive focused on the scientific record and prevent disclosure of submission-process or workstation information.
