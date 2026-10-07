# Implementation audit before differential analysis

Date: 2026-09-16 (Asia/Shanghai)  
Scope: download, input validation, environment construction, and six-cohort QC only. No differential-expression or candidate-gene result was computed or inspected.

## Recorded implementation corrections

1. The GSE76427 non-normalized matrix labels adjacent non-tumour samples as `ANTT`, not `NT`. The parser was corrected to the deposited label; the frozen GSM manifest was unchanged.
2. GEO does not publish `annot/*.annot.gz` files for GPL5918 or GPL16699. Their official GEO full platform-text endpoints were used instead; cohort and sample sources were unchanged.
3. The conda `preprocessCore` binary could not create pthreads on this host. The same Bioconductor version (1.68.0) was rebuilt from official source with `--disable-threading`. The algorithm and frozen GPL570 RMA method were unchanged. Source URL, SHA-256, and build marker are retained.
4. GSE45114 legacy GPR headers include non-UTF-8 scanner paths and inconsistent quoting. All feature blocks were verified as 23,232 rows by 82 columns, and aligned by the stable Block/Column/Row coordinates. Only the feature blocks are parsed; no sample was excluded for header encoding.
5. A preliminary QC implementation allowed post-quantile median/IQR robust Z scores to create flags. Because quantile normalization deliberately equalizes distributions and leaves only tiny rank/tie residuals, these are not valid independent exclusion diagnostics. Before any differential analysis, distribution values were retained for description but their flags were disabled for RMA-, Agilent-quantile-, and Illumina-quantile-normalized data. GPR M-value distribution flags remain eligible. This correction is conservative: it prevents a normalization artefact from excluding a sample.

All code corrections were based on file structure or computational diagnostics, not disease-group differences, candidate genes, P values, or downstream effect estimates.
