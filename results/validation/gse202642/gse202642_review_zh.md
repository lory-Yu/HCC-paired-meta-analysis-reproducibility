Takeaway: GSE202642 完成组织样本级 pseudobulk 定位；97,255 个 QC 细胞、7 肿瘤 + 4 癌旁样本；291/292 候选可检出；high_confidence 主要落在胆管样/肝细胞/Kupffer/CAF。

# GSE202642 单细胞定位审阅

本阶段在 Phase 2（TCGA/ICGC）审阅通过后运行。仅使用锁定候选集 v2（173 high_confidence + 119 intermediate）；未回改发现阈值。统计单位为组织样本（GSM），不是单个细胞。

## 数据与口径

- 队列：GSE202642（HBV 相关 HCC；PMID 36878933）。
- 矩阵：合并 10x MTX；条码后缀 `-1`…`-11` 对应 11 个 GSM（7 肿瘤 + 4 癌旁）。
- 原始细胞：115,732；QC 后：97,255（nCount≥500、nFeature 200–6000、线粒体比例≤20%）。
- 细胞类型：marker 评分 + Louvain 聚类映射为 10 个广义类型。
- 独立单位：`tissue_sample_GSM_not_cell`。
- 肿瘤—癌旁推断门槛：每个组织类别至少 3 个样本单位；不足则记 `insufficient_sample_units`。
- 候选集 SHA-256：high `55110e…`；intermediate `23812a…`（未改）。

## 定位结果

- 292 个锁定候选中 **291** 在矩阵中可检出（1 个缺失）。
- high_confidence（173）主导细胞类型（按最大平均表达、≥50 细胞）：

| 主导类型 | 基因数 |
|---------|--------|
| Cholangiocyte | 57 |
| Hepatocyte | 43 |
| Kupffer_macrophage | 21 |
| Fibroblast_CAF | 18 |
| Endothelial / B_cell / Neutrophil | 各 7 |
| T_cell | 6 |
| Plasma | 4 |
| NK_cell | 2 |

## 肿瘤—癌旁方向（样本级，次要）

在主导细胞类型内、满足 ≥3 肿瘤样本与 ≥3 癌旁样本时：

- high_confidence 可检验：104 / 173；不足单位：68。
- 可检验中与发现 Meta 同向：51；反向：53。
- 主导类型内 Wilcoxon 候选集多重校正后 FDR&lt;0.05：**0**（样本数小，检验力有限）。

因此本阶段**主要支持细胞类型定位**；肿瘤—癌旁方向在单细胞伪批量层仅为探索性，不能替代 TCGA/ICGC bulk 配对验证。

## 机器可读产出

- `candidate_celltype_localization.tsv`
- `candidate_dominant_celltype.tsv`
- `candidate_patient_pseudobulk.tsv.gz`
- `candidate_celltype_tumor_adjacent_pseudobulk.tsv`
- `cell_type_marker_table.tsv`、`cell_metadata.tsv.gz`
- `figures/`（UMAP、feature、主导类型条形图）
- `download_checksums.tsv`（矩阵 SHA-256 `855789bd…`）

## 允许与禁止的表述

- 允许：候选基因主要定位于胆管样上皮、肝细胞、Kupffer/CAF 等；结果基于患者/组织样本级 pseudobulk。
- 禁止：把细胞数写成患者数；把单细胞方向不一致写成推翻 TCGA/ICGC；或声称治疗获益、直接靶点及因果机制。

## 局限

- 仅 11 个组织样本（7+4），肿瘤—癌旁推断普遍受限。
- 细胞类型为 marker 广义注释，非作者官方标签复现。
- 未纳入 DepMap / LINCS；HPA/ChEMBL 另见对应目录。
