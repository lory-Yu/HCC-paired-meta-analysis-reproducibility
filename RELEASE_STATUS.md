# Release status and public-upload gate

## Current state

This directory is the locally assembled source of the checksum-verified public reproducibility-core release. The public repository is `https://github.com/lory-Yu/HCC-paired-meta-analysis-reproducibility`. Version `v0.1.1-submission` adds Supplementary Table S1, the GEO metadata audit and revised manuscript text and figures; its version DOI is recorded here after Zenodo publishes it.

## Included

- HCC/CBC protocols, change-control records and identifier-harmonization rules;
- original analysis and verification code plus the locked environment;
- derived meta-analysis, sensitivity, TCGA/ICGC and single-cell tables;
- GEO sample-characteristics field audit and gene-level Supplementary Table S1;
- figure source data, figures, manifests and provenance checksums.

## Excluded

- raw GEO, TCGA-LIHC and ICGC-LIRI-JP files;
- full GEO SOFT exports (about 540 MB; contain third-party expression tables);
- credentials, access tokens, personal data and temporary editor files;
- all projects and experimental templates outside the submitted HCC article;
- the private reviewer ZIP and its access instructions.

## Author approval and remaining publication metadata

1. The authors approved public release on 7 October 2026: original code under MIT and original derived outputs under CC BY 4.0.
2. Verified ORCID identifiers were not supplied and are therefore left blank rather than inferred.
3. The concept DOI is `10.5281/zenodo.23210221`; the `v0.1.1-submission` version DOI is pending publication.
4. Version-specific DOI links are backfilled only after the Zenodo record becomes publicly resolvable.

## Integrity

The root `SHA256SUMS.sha256` file is generated after all release files are finalized. It is the integrity check for this release, not the original workspace checksum list.
