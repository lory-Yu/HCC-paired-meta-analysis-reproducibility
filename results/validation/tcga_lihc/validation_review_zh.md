Takeaway: TCGA-LIHC 配对验证只检验已锁定的候选集 v2，不回改发现阈值。

# TCGA-LIHC 独立验证审阅

- 表达矩阵：UCSC Xena GDC Hub `TCGA-LIHC.star_counts.tsv.gz`（xena_gdc_log2_star_count_plus_1）。
- 表型：`TCGA-LIHC.clinical.tsv.gz`，使用官方 `sample_type.samples`。
- 完整肿瘤—正常配对数：50。
- 主验证集：173 个 high_confidence；次级：119 个 intermediate。
- FDR：high_confidence 与 intermediate 分别在各自候选集内做 BH，不把 17,900 基因组重新筛选。
- high_confidence 方向一致：170；候选集内 FDR<0.05：159；validated_primary：158。
- intermediate 方向一致：115；候选集内 FDR<0.05：108。
- 映射后缺失基因数：2。
- 未配对分析和临床关联仅为次级，不能替代配对主验证，也不能写成疗效。
- 本阶段不支持治疗获益、直接靶点结合或因果机制表述。
