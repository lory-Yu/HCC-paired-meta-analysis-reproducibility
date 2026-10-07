# 冻结分析协议

本目录保存正式重分析开始前的不可变科学协议和样本登记。

- `analysis_plan_v1.0.md`：供人工审阅的完整方案。
- `analysis_freeze_v1.0.yaml`：供工作流读取的冻结参数。
- `cohort_manifest_v1.0.tsv`：队列级样本与平台清单。
- `sample_manifest_v1.0.tsv`：逐 GSM 纳入、排除、患者和配对清单。
- `source_metadata/`：NCBI GEO 官方 SOFT 导出。
- `scripts/parse_geo_soft.py`：从官方 SOFT 重建样本清单的确定性脚本。
- `change_control.md`：冻结后修订记录。

现有旧分析结果只用于识别数据完整性问题，不用于改变本协议的主要终点或阈值。
