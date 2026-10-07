# Zenodo/GitHub 公开交接清单

## 当前完成状态

- 已整理纯 HCC/CBC 可复现核心版；
- 已排除第三方原始数据、无关项目、投稿私密包和本机绝对路径；
- 已准备 Zenodo 元数据、引用文件、公开后稿件声明和正式许可证；
- 已执行发现分析、稳健性、双 FDR、ICGC 标识符调和和定点审计验证；
- 尚未实际上传，因此当前仍无公开 URL 或 DOI。

## 已确认的公开授权

1. 原创代码正式采用 MIT License。
2. 原创派生数据、图源数据和文档正式采用 CC BY 4.0。
3. 作者已于 2026 年 10 月 7 日授权公开；第三方原始数据未被再分发或重新许可。
4. 两位作者的 ORCID 未提供，在 Zenodo 留空，不能猜测。
5. Zenodo 联系邮箱使用通讯作者 `wangzp@jiangnan.edu.cn`。

## 建议发布顺序

1. 在 GitHub/GitLab 创建公开仓库，上传解压后的 release 内容并建立标签 `v0.1.0-submission`。
2. 在 Zenodo 新建 Dataset 记录，上传根目录 ZIP；也可启用 GitHub–Zenodo 集成归档同一标签。
3. 使用 `ZENODO_METADATA_TEMPLATE.yaml` 填写元数据；不要填写尚不存在的文章 DOI。
4. 发布 Zenodo 记录，取得可以解析的真实 DOI。
5. 将 DOI 和公开代码 URL 回填到 `CITATION.cff`、README、稿件 Data availability、Code availability 及参考文献。
6. DOI 回填属于文本更新；若 release 文件本身改变，应发布新版本并重新生成校验和，而不是覆盖已公开文件。
