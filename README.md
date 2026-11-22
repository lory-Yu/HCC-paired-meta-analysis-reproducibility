# HCC reproducibility core release

This pre-submission research compendium accompanies *Paired cross-platform transcriptomic analysis with dual-cohort validation identifies recurrent tissue-level alterations in hepatocellular carcinoma* by Simiao Yu and Zhouping Wang.

## Version

`v0.2.0-submission` — a public research-compendium release prepared for the Medical Oncology submission. It adds a locked unambiguous identifier analysis, explicit 173 → 167 → 163 denominator accounting, candidate-family and genome-wide FDR sensitivity, a donor-level single-cell sensitivity analysis, revised figures and release-specific verification. It is not an accepted-manuscript archive.

## Contents

- `00_protocol/`: frozen cohort manifests, multiplicity rules, identifier audits and change control.
- `workflow/`: executable analysis and verification scripts plus environment files.
- `results/analysis/`: probe annotation, paired cohort effects, complete meta-analysis tables and prespecified sensitivity analyses.
- `results/validation/`: derived TCGA/ICGC validation tables, identifier-harmonization audits, Allain comparison, descriptive single-cell localization and the `medical_oncology_revision_v1` audit set.
- `results/audit/`: GEO sample-characteristics field audit (script output; rerunning needs re-downloaded GEO SOFT exports).
- `results/provenance/`: download manifests, environment snapshots and provenance checksums.
- `manuscript/source_data/`: source data for all manuscript figures and tables.
- `manuscript/figures/`: editable and rendered figure assets, including the revised Fig. 1, Fig. 2 and Fig. 5.
- `manuscript/supplementary_table_S1_gene_level_CBC.{tsv,xlsx}`: gene-level discovery and external-validation results for all 446 discovery hits.

## Primary data

No third-party raw data are redistributed. GEO, TCGA-LIHC and ICGC-LIRI-JP accession information, versions, download locations and checksums are included in the manifests. Users must obtain those files from the original repositories and comply with their terms. The included matrices and tables are derived research outputs; they are not a replacement for the source repositories.

## Reproduction

Create the environment from `workflow/environment.yml` or the explicit Linux lock, reacquire the public source data using the recorded manifests, and run the documented scripts or the public Makefile. Verification targets do not require the raw files. First verify the release itself:

```bash
sha256sum -c SHA256SUMS.sha256
```

The public Makefile exposes integrity and verification targets without the author's local absolute paths. Full reruns require the public source data and a local R/Python environment. The original discovery protocol also contained an abandoned pharmacology extension; that mixed-scope file is intentionally excluded. `00_protocol/HCC_CBC_PUBLIC_ANALYSIS_SPEC.md` is the public HCC-only specification and points to the frozen cohort and candidate files used here.

Figures, Supplementary Table S1 and the GEO metadata audit can be rebuilt with `make cbc-figures` (requires matplotlib 3.11.2), `make cbc-supplementary-s1` and `GEO_SOFT_DIR=<dir> make cbc-geo-metadata-audit`. Rebuilding changes PDF/SVG/XLSX timestamps, so run these in a copy if you also check `SHA256SUMS.sha256`.

The v0.2.0 denominator and donor-sensitivity assertions can be run with `make verify-medical-oncology-revision`; the release-increment privacy check is `make privacy-audit`. `make verify` runs both together with the existing integrity and analysis checks.

## Licenses

- Original analysis and verification code: MIT License (`LICENSE_CODE.txt`).
- Original derived tables, figure source data and documentation: CC BY 4.0 (`LICENSE_DERIVED_DATA.txt`).
- Third-party source data are not redistributed or relicensed and remain governed by their source repositories.

## Citation

The public code repository is `https://github.com/lory-Yu/HCC-paired-meta-analysis-reproducibility`. Release `v0.2.0-submission` will receive a version-specific DOI when Zenodo completes the GitHub release archive; the stable concept DOI for all versions is `https://doi.org/10.5281/zenodo.23210221`.

Yu S, Wang Z. Paired cross-platform transcriptomic analysis with dual-cohort validation identifies recurrent tissue-level alterations in hepatocellular carcinoma: research compendium. Version v0.2.0-submission. Zenodo; 2026. Concept DOI: 10.5281/zenodo.23210221.

If the manuscript changes materially, create a new version rather than overwriting this release.

## Contact

Zhouping Wang, `wangzp@jiangnan.edu.cn`.
