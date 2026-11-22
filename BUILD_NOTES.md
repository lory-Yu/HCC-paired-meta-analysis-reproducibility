# Build notes

This directory was assembled from the HCC project workspace on 2026-10-07 and updated on 2026-11-22 for `v0.2.0-submission`.

## Selection rule

Included files must support at least one manuscript result, figure, statistical denominator, identifier audit or reproducibility claim. This release adds the primary unambiguous identifier analysis, explicit external-validation denominators, multiplicity sensitivity and donor-level single-cell sensitivity needed by the Medical Oncology submission.

## Deliberate exclusions

- third-party raw GEO, TCGA, ICGC and single-cell input files;
- local caches, R libraries, credentials and tokens;
- private reviewer access material;
- Word/PDF submission files and cover letters;
- editorial correspondence, decision letters and internal response reports;
- full GEO SOFT exports used by the metadata audit;
- projects and endpoints outside the submitted HCC transcriptomic article;
- early mixed-scope protocol documents superseded for this release by `00_protocol/HCC_CBC_PUBLIC_ANALYSIS_SPEC.md`.

## Path normalization

Author-machine absolute paths in copied manifests were converted to project-relative paths. Repository URLs, accession identifiers, original file sizes and source-file SHA-256 values were retained.

## Validation

The public Makefile verifies package checksums, discovery/meta-analysis outputs, sensitivity outputs, both external-validation multiplicity families, identifier harmonization, the Medical Oncology revision denominators and the release-increment privacy boundary. Verification results are recorded in `VERIFICATION_REPORT.txt`.

From `v0.1.1-submission` the Makefile also rebuilds the figures (`cbc-figures`, with a matplotlib 3.11.2 guard), Supplementary Table S1 (`cbc-supplementary-s1`, with internal count assertions) and the GEO metadata audit (`cbc-geo-metadata-audit`, requires `GEO_SOFT_DIR`). Submission-document builders and generated Word/PDF files are excluded from `v0.2.0-submission`; this archive is limited to data, code, figures and research documentation.
