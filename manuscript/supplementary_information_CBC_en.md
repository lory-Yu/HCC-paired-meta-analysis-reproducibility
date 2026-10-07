# Supplementary Information

## Supplementary Table S1. Gene-level discovery and external validation results

Supplied as a separate spreadsheet, `supplementary_table_S1_gene_level_CBC.xlsx`; its README sheet defines every column. The table lists all 446 discovery hits (173 discovery robust, 119 intermediate and 154 exploratory) with pooled tumour-minus-non-tumour log2 fold changes, 95% confidence and prediction intervals, Hartung–Knapp P values, discovery FDR and heterogeneity statistics. For the 292 discovery robust and intermediate genes carried into external validation, it also gives the TCGA-LIHC and ICGC-LIRI-JP paired effects with 95% confidence intervals, nominal P values, candidate-family and genome-wide FDR, directional calls, reasons for non-estimability, and the ICGC matrix symbol with its mapping basis. Exploratory genes were not tested externally. Every count in main-text Table 1 is reproduced from these rows.

## Supplementary Table S2. Non-estimability after identifier harmonization

### TCGA-LIHC

Five discovery robust genes were present before filtering but did not meet the complete-matrix expression filter: **IGF2BP3, SRXN1, COX7B2, FAM133A and PAGE4**. The expression filter was not relaxed after these genes were identified. A count-level audit confirmed that all five failed the unchanged cohort-wide `filterByExpr` rule; for example, IGF2BP3 had a median of 3 counts in normal samples and 21.5 counts in tumour samples, with CPM ≥1 in only 17/100 libraries. These genes therefore remain not estimable in TCGA rather than being labelled directionally discordant.

### ICGC-LIRI-JP

Thirty-four candidates absent under their current symbols were recovered under frozen historical symbols. **REXO5** remained absent after harmonization. **COX7B2** was present but removed by the expression filter. The current symbol, historical matrix symbol, Entrez ID and mapping basis for every candidate are listed in Supplementary Table S1.

## Supplementary Table S3. Descriptive external concordance statistics

| Scope | Cohort | Direction | Exact 95% CI | Pearson r | Spearman rho | OLS slope (95% CI) |
|---|---|---:|---:|---:|---:|---:|
| All estimable | TCGA-LIHC | 167/168 | 0.967–1.000 | 0.914 | 0.870 | 1.73 (1.61–1.85) |
| Identifier-harmonized | ICGC-LIRI-JP | 171/171 | 0.979–1.000 | 0.897 | 0.881 | 1.57 (1.46–1.69) |
| Jointly estimable | TCGA-LIHC | 166/167 | 0.967–1.000 | 0.913 | 0.871 | 1.73 (1.61–1.85) |
| Jointly estimable | ICGC-LIRI-JP | 167/167 | 0.978–1.000 | 0.924 | 0.878 | 1.50 (1.41–1.60) |
| Jointly estimable, both cohorts | Both cohorts | 166/167 | 0.967–1.000 | not applicable | not applicable | not applicable |

Direction summaries report proportions and exact Clopper–Pearson 95% confidence intervals. They are descriptive because genes are correlated and the candidate set was selected before this audit; genes are not independent biological replicates. No binomial null test is reported. Ordinary least-squares slopes compare platform-specific effect scales and are not calibration coefficients or causal effects.

## Supplementary Note 1. Descriptive single-cell localization

GSE202642 supplied 115 732 cells from 11 tissue samples. After quality control (minimum 500 counts, 200–6000 detected features and no more than 20% mitochondrial reads), 97 255 cells remained from seven tumour and four adjacent-tissue samples. Ten broad classes were retained: B cell, cholangiocyte, endothelial, fibroblast/cancer-associated fibroblast, hepatocyte, Kupffer/macrophage, NK cell, neutrophil, plasma cell and T cell.

Of the 292 locked discovery robust and intermediate genes, 291 were detected. Among the 173 discovery robust genes, the dominant broad classes were cholangiocyte (57), hepatocyte (43), Kupffer/macrophage (21), fibroblast/CAF (18), endothelial (7), B cell (7), neutrophil (7), T cell (6), plasma cell (4) and NK cell (2). These assignments are descriptive and are based on broad-class mean expression. An exploratory sample-level tumour-versus-adjacent comparison was performed, but the seven-versus-four design was insufficient for reliable inference. Its direction counts were therefore not retained as evidence or reported as differential validation; cells were not treated as independent samples.

![Supplementary Figure S1. Single-cell localization](figures/FigS1_single_cell_localization_CBC.png)

**Supplementary Figure S1.** (a) Broad cell classes among 97 255 quality-controlled cells. (b) Dominant localization of discovery robust genes. (c) Example expression localization heat map. Localization is descriptive and is not an independent differential-expression validation. Heat-map abbreviations are: B, B cell; Chol, cholangiocyte; Endo, endothelial; CAF, fibroblast/cancer-associated fibroblast; Hep, hepatocyte; Kup, Kupffer/macrophage; NK, natural killer cell; Neut, neutrophil; Plas, plasma cell; and T, T cell.

## Supplementary Table S4. Descriptive comparison with the Allain et al. 935-gene signature

| Current analysis set | Current denominator | Allain denominator | Overlap | Same direction among overlapping genes |
|---|---:|---:|---:|---:|
| Discovery robust | 173 | 935 | 39 | 39/39 |
| Jointly estimable in TCGA and ICGC | 167 | 935 | 38 | 38/38 |

The comparison used the official Allain et al. Supplementary Table 2 archive (dataset DOI: 10.1158/0008-5472.22412988.v1). It is descriptive and is not counted as independent validation because study-level overlap between public microarray compendia cannot be excluded. No overlap-enrichment or direction-null P value was calculated.
