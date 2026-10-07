Takeaway: ICGC-LIRI-JP 使用 Release 28 公开 exp_seq，并完成真正的 donor 配对验证。

# ICGC-LIRI-JP 独立验证审阅

- 项目确认：ICGC-LIRI-JP（开放桶 `icgc25k-open/release_28/data/LIRI-JP/`）。
- 未用其他 ICGC 肝癌项目、TCGA 或 GEO 冒充。
- 配对：是。完整肿瘤—正常 donor 数 = 199。
- 未发生 unpaired fallback。
- 表达单位：log2(raw_read_count + 1)；基因标识为 RefSeq gene_id（符号）。
- 主验证集：173 high_confidence；次级：119 intermediate。
- high_confidence 方向一致：118；候选集内 FDR<0.05：114；validated_primary：113。
- intermediate 方向一致：94；候选集内 FDR<0.05：91。
- 缺失/不可估计基因：0。
- 结果支持跨队列表达关联的可重复性；不支持治疗获益、直接靶点或因果机制表述。
