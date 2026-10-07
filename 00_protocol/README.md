# 冻结分析协议

本目录保存公开复现包所需的冻结科学协议、样本登记和 CBC 修订审计规则。

- `cohort_manifest_v1.0.tsv`：队列级样本与平台清单。
- `sample_manifest_v1.0.tsv`：逐 GSM 纳入、排除、患者和配对清单。
- `change_control.md`：冻结后修订记录。
- `HCC_CBC_PUBLIC_ANALYSIS_SPEC.md`：HCC-only 公开分析规范及冻结文件索引。
- `cbc_external_validation_multiplicity_freeze_v1.0.md`：外部验证两套多重检验口径。
- `cbc_icgc_identifier_harmonization_freeze_v1.0.md` 与 `cbc_icgc_historical_symbol_map_v1.tsv`：ICGC 历史符号调和规则和映射表。

原始 GEO SOFT 导出不随公开复现包分发；其获取方式和元数据字段审计的重跑要求见根目录 README 和相应审计脚本。

现有旧分析结果只用于识别数据完整性问题，不用于改变本协议的主要终点或阈值。
