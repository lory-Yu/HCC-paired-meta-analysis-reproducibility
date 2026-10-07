# CBC 首投包导航

本目录是针对 *Computational Biology and Chemistry* 的首投包。它不是对 Functional & Integrative Genomics 原稿的简单换刊，而是根据编辑提出的四点意见重构了外部验证统计口径，并删减了不能形成主证据链的阴性下游分析。

## 先看什么

1. [revision_audit_zh.md](revision_audit_zh.md)：中文修改说明。
2. [manuscript_CBC_submission.docx](manuscript_CBC_submission.docx)：可编辑英文主稿。
3. [supplementary_information_CBC.docx](supplementary_information_CBC.docx)：补充表、单细胞定位及端点处置。
4. [cover_letter_CBC.docx](cover_letter_CBC.docx) 和 [highlights_CBC.txt](highlights_CBC.txt)：投稿附加文件。
5. `figures/`：两个主图、一个补充图和图形摘要的 PDF/SVG/PNG/TIFF；图件经过文字大小、碰撞和版面审计。
6. `results/validation/cbc_targeted_revision_audits_v1/`：Allain 2016 定量对照和 GEO–TCGA/ICGC 公开样本标识符重叠审计。
7. 补充表 S6：按共同分母重建可估性、标识符调和、方向一致性和两套 FDR 的诊断审计；它不引入新的候选筛选。

## 统计主线

- 发现阶段锁定：439 个 intent-to-analyse 样本；QC 后 438 个；167 对；17,900 个基因；446 个发现命中；173 个 discovery robust、119 个 intermediate、154 个 exploratory。
- TCGA-LIHC：50 对；完整过滤转录组 15,657 个基因；discovery robust 168/173 可估计，161/168 通过全转录组 BH-FDR 且方向一致。
- ICGC-LIRI-JP：199 对；完整过滤转录组 19,353 个基因。对历史符号进行效应盲法统一后，171/173 可估计，170/171 通过全转录组 BH-FDR 且方向一致。
- 两库共同可估计：167/173；锁定候选家族的主要双库标准：161/167；全转录组 BH-FDR 敏感性标准：159/167。

161/167 不是合并两个数据库后重新计算的联合 FDR，而是两个外部队列分别达到候选家族 BH-FDR<0.05 且方向一致后的交集；159/167 是两个数据库分别进行全转录组 BH 校正后的敏感性结果。原始符号直接匹配得到的 127/133 和 125/133 只作为技术敏感性分析。不可估计基因不计为阴性。

## 被移出主文的内容

GSE202642 只用于描述信号在肝细胞、胆管上皮细胞、Kupffer 巨噬细胞、内皮细胞和成纤维细胞等区室中的定位。由于只有 7 个肿瘤和 4 个癌旁样本，样本级差异检验不进入稿件证据链。蛋白、药理、扰动和其他超出转录组复现范围的端点均不进入本次投稿稿件。

## 数据与复现

统计重构结果位于 `results/validation/cbc_external_validation_v1/`，通过 `make verify-cbc-external-validation` 和 `make verify-cbc-icgc-identifier-sensitivity`。源码包括外部验证、历史符号统一、效应量一致性与作图脚本。`source_data_manifest_CBC.tsv` 记录每个图件源表的 SHA-256、用途和上游路径；主交付文件的校验和见 [delivery_checksums_CBC.sha256](delivery_checksums_CBC.sha256)。

当前未伪写公开仓库或 DOI。如作者决定在首投前公开，需先确认授权和许可证，完成 GitHub/Zenodo 归档后再将真实 URL/DOI 写入 Data availability 和 Code availability。
