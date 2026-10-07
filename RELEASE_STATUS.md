# Release status and public-upload gate

## Current state

This directory is the locally assembled source of the checksum-verified public reproducibility-core release. The public repository is `https://github.com/lory-Yu/HCC-paired-meta-analysis-reproducibility`; tag `v0.1.0-submission` was archived by Zenodo as `https://doi.org/10.5281/zenodo.23210222` on 7 October 2026.

## Included

- HCC/CBC protocols, change-control records and identifier-harmonization rules;
- original analysis and verification code plus the locked environment;
- derived meta-analysis, sensitivity, TCGA/ICGC and single-cell tables;
- figure source data, figures, manifests and provenance checksums.

## Excluded

- raw GEO, TCGA-LIHC and ICGC-LIRI-JP files;
- credentials, access tokens, personal data and temporary editor files;
- all projects and experimental templates outside the submitted HCC article;
- the private reviewer ZIP and its access instructions.

## Author approval and remaining publication metadata

1. The authors approved public release on 7 October 2026: original code under MIT and original derived outputs under CC BY 4.0.
2. Verified ORCID identifiers were not supplied and are therefore left blank rather than inferred.
3. GitHub and Zenodo publication completed; the version DOI is `10.5281/zenodo.23210222` and the concept DOI is `10.5281/zenodo.23210221`.
4. DOI and repository links were backfilled only after the Zenodo record became publicly resolvable.

## Integrity

The root `SHA256SUMS.sha256` file is generated after all release files are finalized. It is the integrity check for this release, not the original workspace checksum list.
