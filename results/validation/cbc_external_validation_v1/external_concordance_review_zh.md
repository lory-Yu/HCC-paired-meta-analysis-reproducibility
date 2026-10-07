# External concordance audit (CBC v1)

This is a derived, descriptive audit of the locked post-rejection external-validation tables. It does not alter the discovery analysis, candidate locks, or manuscript files.

## Statistical definitions

- Direction concordance is the number of external log2FC estimates having the same sign as the frozen discovery log2FC. Exact binomial tests use a 0.5 null; intervals are two-sided 95% Clopper--Pearson intervals.
- Effect-size concordance is assessed with Pearson correlation, Spearman correlation, and ordinary least-squares regression of external log2FC on discovery log2FC.
- The `all_estimable` scope uses every discovery-robust gene estimable in that external cohort. The `joint_133` scope uses only genes estimable in both TCGA-LIHC and ICGC-LIRI-JP.
- The gene-level binomial p values are descriptive: genes are correlated and the locked candidate set was selected before this audit. They are not evidence that genes are independent biological replicates.

## Direction results

| scope | cohort | metric | success/total | proportion | CP 95% CI | exact two-sided p |
|---|---|---|---:|---:|---:|---:|
| all_estimable | TCGA-LIHC | same_direction | 167/168 | 0.9940 | 0.9673–0.9998 | 9.034e-49 |
| all_estimable | ICGC-LIRI-JP | same_direction | 137/137 | 1.0000 | 0.9734–1.0000 | 1.148e-41 |
| joint_133 | TCGA-LIHC | same_direction | 132/133 | 0.9925 | 0.9588–0.9998 | 2.461e-38 |
| joint_133 | ICGC-LIRI-JP | same_direction | 133/133 | 1.0000 | 0.9726–1.0000 | 1.837e-40 |
| joint_133 | TCGA-and-ICGC | dual_same_direction | 132/133 | 0.9925 | 0.9588–0.9998 | 2.461e-38 |
| joint_133 | TCGA-and-ICGC | dual_candidate_family_FDR_and_same_direction | 127/133 | 0.9549 | 0.9044–0.9833 | 1.321e-30 |
| joint_133 | TCGA-and-ICGC | dual_genomewide_FDR_and_same_direction | 125/133 | 0.9398 | 0.8849–0.9737 | 3.840e-28 |
| identifier_harmonized_all_estimable | ICGC-LIRI-JP | same_direction | 171/171 | 1.0000 | 0.9787–1.0000 | 6.682e-52 |
| identifier_harmonized_joint_167 | TCGA-LIHC | same_direction | 166/167 | 0.9940 | 0.9671–0.9998 | 1.796e-48 |
| identifier_harmonized_joint_167 | ICGC-LIRI-JP | same_direction | 167/167 | 1.0000 | 0.9782–1.0000 | 1.069e-50 |
| identifier_harmonized_joint_167 | TCGA-and-ICGC | dual_same_direction | 166/167 | 0.9940 | 0.9671–0.9998 | 1.796e-48 |
| identifier_harmonized_joint_167 | TCGA-and-ICGC | dual_candidate_family_FDR_and_same_direction | 161/167 | 0.9641 | 0.9234–0.9867 | 3.054e-40 |
| identifier_harmonized_joint_167 | TCGA-and-ICGC | dual_genomewide_FDR_and_same_direction | 159/167 | 0.9521 | 0.9078–0.9791 | 1.424e-37 |

## Effect-size results

| scope | cohort | n | Pearson r | Spearman rho | OLS slope (95% CI) | intercept (95% CI) | R² |
|---|---|---:|---:|---:|---:|---:|---:|
| all_estimable | TCGA-LIHC | 168 | 0.9137 | 0.8699 | 1.7328 (1.6146–1.8509) | -0.2452 (-0.3600–-0.1304) | 0.8348 |
| all_estimable | ICGC-LIRI-JP | 137 | 0.8914 | 0.8821 | 1.6286 (1.4876–1.7695) | -0.0432 (-0.1825–0.0961) | 0.7946 |
| joint_133 | TCGA-LIHC | 133 | 0.9119 | 0.8721 | 1.7901 (1.6509–1.9294) | -0.2682 (-0.4055–-0.1308) | 0.8316 |
| joint_133 | ICGC-LIRI-JP | 133 | 0.9203 | 0.8775 | 1.5404 (1.4272–1.6536) | -0.1332 (-0.2449–-0.0215) | 0.8469 |
| identifier_harmonized_all_estimable | ICGC-LIRI-JP | 171 | 0.8970 | 0.8810 | 1.5735 (1.4557–1.6913) | -0.0559 (-0.1708–0.0589) | 0.8045 |
| identifier_harmonized_joint_167 | TCGA-LIHC | 167 | 0.9132 | 0.8712 | 1.7312 (1.6125–1.8499) | -0.2472 (-0.3627–-0.1316) | 0.8340 |
| identifier_harmonized_joint_167 | ICGC-LIRI-JP | 167 | 0.9245 | 0.8776 | 1.5014 (1.4062–1.5965) | -0.1278 (-0.2204–-0.0352) | 0.8547 |

The slopes are descriptive scale comparisons, not calibration coefficients or causal effects. Values above one indicate that the external effect estimates are larger in magnitude on average than the discovery estimates under the respective platform and model specifications.

## Provenance

- Script: `workflow/scripts/compute_cbc_external_concordance_v1.py`
- Input SHA-256: `tcga_locked_candidate_multiplicity.tsv` = `76665fbf1ad03852e689684fefb18702b1ea954784d56d2e83ccffe3bc27e71c`
- Input SHA-256: `icgc_locked_candidate_multiplicity.tsv` = `f188a6767a378f0ea89df7795957eacb1c9756de7b0621130fde14fb3216332c`
- Input SHA-256: `tcga_icgc_locked_candidate_cross_validation.tsv` = `36cd16d071d72d0793618813aa410e4e64e6f22e17c18b0b83d98ae400727d9c`
- Input SHA-256: `icgc_locked_candidate_identifier_harmonized_sensitivity.tsv` = `c2d904b62b7952dd4e020cf91f5c21a592db334d36a840c226191ce6290b9481`
- Input SHA-256: `tcga_icgc_identifier_harmonized_sensitivity.tsv` = `e9c8aa745e91c940e3ede3ad00c3da6d2f491bab62be027a2ac4897e3b5e7da8`
