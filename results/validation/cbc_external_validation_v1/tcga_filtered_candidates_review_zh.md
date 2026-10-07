# TCGA 低表达候选审计

本审计只解释 5 个在 TCGA 原始矩阵中存在、但被预设 `filterByExpr` 规则过滤的发现稳健基因。
未放宽过滤阈值、未重新拟合被过滤基因，也未根据 ICGC 或发现阶段效应恢复候选。

- 审计基因数：5
- filterByExpr 保留数：0
- 结论：这些基因在 TCGA 中继续记为 not estimable；低表达过滤不是方向不一致。

详细的非零样本数、中位计数和 CPM 见 `tcga_filtered_discovery_robust_candidates_audit.tsv`。
