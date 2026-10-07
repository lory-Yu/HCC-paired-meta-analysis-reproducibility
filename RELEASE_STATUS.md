# Release status and public-upload gate

## Current state

This directory is a locally assembled, checksum-verified public reproducibility-core release. It has not been uploaded to GitHub, Zenodo or another public repository by this workflow, and therefore has no public URL or DOI.

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
3. Upload the exact checksum-verified archive, record its final SHA-256, and obtain the real DOI.
4. Replace the DOI placeholder in the public README and manuscript Data/Code availability only after publication of the repository record.

## Integrity

The root `SHA256SUMS.sha256` file is generated after all release files are finalized. It is the integrity check for this release, not the original workspace checksum list.
