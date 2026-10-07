# Build notes

This directory was assembled from the HCC project workspace on 2026-10-07 as a public reproducibility-core release and updated on 2026-10-08 for `v0.1.1-submission`.

## Selection rule

Included files must support at least one manuscript result, figure, statistical denominator, identifier audit or reproducibility claim. The release includes frozen HCC/CBC protocols, original workflow code, environment specifications, complete derived discovery/sensitivity tables, derived external-validation tables, descriptive single-cell localization outputs, figure source data and provenance manifests.

## Deliberate exclusions

- third-party raw GEO, TCGA, ICGC and single-cell input files;
- local caches, R libraries, credentials and tokens;
- private reviewer access material;
- Word/PDF submission files and cover letters;
- full GEO SOFT exports used by the metadata audit;
- projects and endpoints outside the submitted HCC transcriptomic article;
- early mixed-scope protocol documents superseded for this release by `00_protocol/HCC_CBC_PUBLIC_ANALYSIS_SPEC.md`.

## Path normalization

Author-machine absolute paths in copied manifests were converted to project-relative paths. Repository URLs, accession identifiers, original file sizes and source-file SHA-256 values were retained.

## Validation

The public Makefile verifies package checksums, discovery/meta-analysis outputs, sensitivity outputs, both external-validation multiplicity families, identifier harmonization and presence of the targeted revision audits. Verification results are recorded in `VERIFICATION_REPORT.txt`.

From `v0.1.1-submission` the Makefile also rebuilds the figures (`cbc-figures`, with a matplotlib 3.11.2 guard), Supplementary Table S1 (`cbc-supplementary-s1`, with internal count assertions) and the GEO metadata audit (`cbc-geo-metadata-audit`, requires `GEO_SOFT_DIR`). The Markdown-to-DOCX and delivery-checksum builders are included so that the manuscript deliverables can be reconstructed, while the private cover-letter source and generated Word/PDF files remain excluded.
