# Zenodo/GitHub 公开交接清单

## 当前完成状态

- 已整理纯 HCC/CBC 可复现核心版；
- 已排除第三方原始数据、无关项目、投稿私密包和本机绝对路径；
- 已准备 Zenodo 元数据、引用文件、公开后稿件声明和正式许可证；
- 已执行发现分析、稳健性、双 FDR、ICGC 标识符调和和定点审计验证；
- 已公开 GitHub 仓库并发布标签 `v0.1.0-submission`；Zenodo 已完成归档。
- GitHub：`https://github.com/lory-Yu/HCC-paired-meta-analysis-reproducibility`
- 版本 DOI：`https://doi.org/10.5281/zenodo.23210222`；概念 DOI：`https://doi.org/10.5281/zenodo.23210221`。

## 已确认的公开授权

1. 原创代码正式采用 MIT License。
2. 原创派生数据、图源数据和文档正式采用 CC BY 4.0。
3. 作者已于 2026 年 10 月 7 日授权公开；第三方原始数据未被再分发或重新许可。
4. 两位作者的 ORCID 未提供，在 Zenodo 留空，不能猜测。
5. Zenodo 联系邮箱使用通讯作者 `wangzp@jiangnan.edu.cn`。

## 建议发布顺序

1. GitHub 公开仓库、标签和 Release 已完成。
2. Zenodo 已通过 GitHub 集成归档同一标签并生成可解析 DOI。
3. DOI 和公开代码 URL 已回填到 `CITATION.cff`、README 及稿件 Data/Code availability。
4. 本次 DOI 回填属于发布后的文本更新；若分析或 release 文件发生实质改变，应发布新版本，而不是覆盖已公开版本。
