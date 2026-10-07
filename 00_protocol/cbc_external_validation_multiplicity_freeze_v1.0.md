# CBC 外部验证多重检验口径重构：冻结后变更控制 v1.0

日期：2026-10-06

## 性质与触发原因

这是 Functional & Integrative Genomics 编辑意见触发的拒稿后统计纠正，不是发现阶段预设分析。编辑指出：发现阶段在全基因宇宙控制 FDR，而外部验证仅在预选候选层内控制 FDR，两者不能作为等价标准直接比较。因此，本次在查看既有验证结果后冻结以下分析口径；旧结果和候选锁均保持只读。

## 锁定输入

- `results/analysis/candidate_lock/candidate_set_v2_high_confidence.tsv`，SHA-256 `55110e282bf9d26619de70d6e8c1e5972d050766b034a6fa06009021ae9fd3cd`，173 个基因。本次稿件称为 **discovery robust set**，不再称为跨队列 high-confidence anchor。
- `results/analysis/candidate_lock/candidate_set_v2_intermediate.tsv`，SHA-256 `23812a47c8df31301aca8a9061ecd7da56dd41afa5f214a3fc6e4a585494aa59`，119 个基因。
- TCGA-LIHC：冻结的 50 对肿瘤—正常样本及 Xena `log2(STAR count + 1)` 文件，按既有 v2 流程反变换为非负整数计数。
- ICGC-LIRI-JP：冻结的 199 对供者及既有 v2 生成的完整 raw-count 配对矩阵。矩阵仅作为锁定输入复用，v2 结果目录不得覆盖。

## 锁定模型

两个外部队列分别执行：

1. 在各自完整可检验转录组上使用 `edgeR::filterByExpr` 过滤低表达基因；
2. TMM 归一化；
3. `limma-voom` 均值—方差建模；
4. 患者/供者固定效应；
5. 对比方向固定为 tumour minus normal；
6. 双侧检验，阈值固定为 0.05，不按结果更改。

TCGA 与 ICGC 的过滤后转录组是各队列自身可检验基因集合；它们不被宣称与发现阶段约 17,900 基因构成完全相同的检验宇宙。

## 三个必须分开的显著性口径

每个队列同时、并列报告：

- `nominal_P`：未经多重校正的模型 P 值；
- `candidate_family_FDR`：分别在锁定的 discovery robust 173 基因和 intermediate 119 基因内部，对可估计候选的 P 值使用 BH 校正；
- `genomewide_FDR`：在该队列完整过滤后基因宇宙的所有 P 值上使用 BH 校正。

三者不得互换。候选家族 FDR 作为定向复核口径；全转录组 FDR 作为与发现阶段全局校正更可比、但检验宇宙仍不完全相同的严格敏感性口径。

## 方向与跨库分母

- `same_direction` 要求外部队列 log2FC 与发现阶段 pooled log2FC 同号。
- 外部验证的 nominal、candidate-family 和 genome-wide 计数均同时要求方向一致和相应 P/FDR < 0.05。
- 跨库汇总仅通过锁定候选的 gene symbol 合并，不对两个完整转录组执行 first-hit 或任意映射。
- 必须分别报告：总锁定数、每库可估计数、两库共同可估计数、共同可估计者中的双库同方向数、双库 nominal 数、双库 candidate-family FDR 数和双库 genome-wide FDR 数。
- 不可估计基因不进入该队列的 BH 分母，但必须保留在候选表和审计分母中。

## 禁止事项

- 不改变六个 GEO 队列、发现模型、样本排除、效应方向、候选阈值或候选分层。
- 不依据本次结果增删候选、调整 0.05 阈值或选择更有利的多重校正口径。
- 不覆盖 `results/validation/rnaseq_normalized_v2/` 或更早验证目录。
- 不把拒稿后重构写成预设分析。
- 不将同方向或显著表达关联写成机制、诊断性能、预后效应或治疗获益。

## 输出位置

所有新结果仅写入 `results/validation/cbc_external_validation_v1/`。
