# Release status and public-upload gate

## Current state

This directory is the locally assembled source of release `v0.2.0-submission`. The public repository is `https://github.com/lory-Yu/HCC-paired-meta-analysis-reproducibility`. The tagged GitHub snapshot was archived by Zenodo as `https://doi.org/10.5281/zenodo.24433095` on 22 November 2026. The stable concept DOI is `https://doi.org/10.5281/zenodo.23210221`.

## Included

- HCC/CBC protocols, change-control records and identifier-harmonization rules;
- original analysis and verification code plus the locked environment;
- derived meta-analysis, sensitivity, TCGA/ICGC and single-cell tables;
- the locked unambiguous identifier analysis with 173 → 167 → 163 denominator accounting;
- candidate-family and genome-wide FDR sensitivity results (157 and 155 dual-cohort passes);
- the donor-level single-cell sensitivity (100 testable contrasts; zero FDR passes);
- GEO sample-characteristics field audit and gene-level Supplementary Table S1;
- figure source data, figures, manifests and provenance checksums.

## Excluded

- raw GEO, TCGA-LIHC and ICGC-LIRI-JP files;
- full GEO SOFT exports (about 540 MB; contain third-party expression tables);
- credentials, access tokens, personal data and temporary editor files;
- all projects and experimental templates outside the submitted HCC article;
- journal correspondence, decision letters, cover letters and private access instructions;
- manuscript Word/PDF files containing author contact details;
- the unfiltered revision ZIP, which contains local paths and submission-process material.

## Author approval and remaining publication metadata

1. The authors approved public release on 7 October 2026: original code under MIT and original derived outputs under CC BY 4.0.
2. Verified ORCID identifiers were not supplied and are therefore left blank rather than inferred.
3. The prior public version is `v0.1.1-submission` (`10.5281/zenodo.23225515`); it remains immutable.
4. The v0.2.0 version DOI and repository links were backfilled after the Zenodo record became publicly resolvable and its downloaded archive passed the checksum, statistical and privacy audits.

## Integrity

The root `SHA256SUMS.sha256` file is generated after all release files are finalized. It is the integrity check for this release, not the original workspace checksum list.
