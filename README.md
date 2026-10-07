# HCC reproducibility core release

This pre-submission research compendium accompanies *Paired-cohort meta-analysis with dual RNA-seq replication of recurrent transcriptional alterations in hepatocellular carcinoma* by Simiao Yu and Zhouping Wang.

## Version

`0.1.0-submission` — a public reproducibility-core release prepared before journal submission. It is not the final accepted-manuscript archive.

## Contents

- `00_protocol/`: frozen cohort manifests, multiplicity rules, identifier audits and change control.
- `workflow/`: executable analysis and verification scripts plus environment files.
- `results/analysis/`: probe annotation, paired cohort effects, complete meta-analysis tables and prespecified sensitivity analyses.
- `results/validation/`: derived TCGA/ICGC validation tables, identifier-harmonization audits, Allain comparison and descriptive single-cell localization.
- `results/provenance/`: download manifests, environment snapshots and provenance checksums.
- `manuscript/source_data/`: source data for all manuscript figures and tables.
- `manuscript/figures/`: editable and rendered figure assets used in the manuscript.

## Primary data

No third-party raw data are redistributed. GEO, TCGA-LIHC and ICGC-LIRI-JP accession information, versions, download locations and checksums are included in the manifests. Users must obtain those files from the original repositories and comply with their terms. The included matrices and tables are derived research outputs; they are not a replacement for the source repositories.

## Reproduction

Create the environment from `workflow/environment.yml` or the explicit Linux lock, reacquire the public source data using the recorded manifests, and run the documented scripts or the public Makefile. Verification targets do not require the raw files. First verify the release itself:

```bash
sha256sum -c SHA256SUMS.sha256
```

The public Makefile exposes integrity and verification targets without the author's local absolute paths. Full reruns require the public source data and a local R/Python environment. The original discovery protocol also contained an abandoned pharmacology extension; that mixed-scope file is intentionally excluded. `00_protocol/HCC_CBC_PUBLIC_ANALYSIS_SPEC.md` is the public HCC-only specification and points to the frozen cohort and candidate files used here.

## Licenses

- Original analysis and verification code: MIT License (`LICENSE_CODE.txt`).
- Original derived tables, figure source data and documentation: CC BY 4.0 (`LICENSE_DERIVED_DATA.txt`).
- Third-party source data are not redistributed or relicensed and remain governed by their source repositories.

## Citation

The public code repository is `https://github.com/lory-Yu/HCC-paired-meta-analysis-reproducibility`. The archived submission release is available at `https://doi.org/10.5281/zenodo.23210222`; the concept DOI for all versions is `https://doi.org/10.5281/zenodo.23210221`.

Yu S, Wang Z. Paired-cohort meta-analysis with dual RNA-seq replication of recurrent transcriptional alterations in hepatocellular carcinoma: reproducibility core. Version v0.1.0-submission. Zenodo; 2026. doi:10.5281/zenodo.23210222.

If the manuscript changes materially, create a new version rather than overwriting this release.

## Contact

Zhouping Wang, `wangzp@jiangnan.edu.cn`.
