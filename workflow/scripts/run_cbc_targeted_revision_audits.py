#!/usr/bin/env python3
"""Targeted pre-submission audits requested for the CBC manuscript.

This script performs two bounded, descriptive checks:
1. overlap and direction agreement with the 935-gene Allain et al. signature;
2. public-identifier overlap between discovery GEO samples and the two external cohorts.

Neither check changes the locked discovery or validation results.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/validation/cbc_targeted_revision_audits_v1"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.upper().isin({"TRUE", "T", "1", "YES"})


def allain_overlap() -> tuple[pd.DataFrame, pd.DataFrame]:
    source = ROOT / "data/external/allain2016/allain2016_supplementary_table_2.xlsx"
    allain = pd.read_excel(source, sheet_name="Supplementary Table 2", header=3)
    allain = allain[allain["Regulation"].isin(["Up in HCC", "Down in HCC"])].copy()
    allain["gene_symbol"] = allain["Gene symbol"].astype(str).str.strip().str.upper()
    allain = allain[["gene_symbol", "Regulation"]].rename(columns={"Regulation": "allain_direction"})
    if len(allain) != 935 or allain["gene_symbol"].nunique() != 935:
        raise RuntimeError("The official Allain signature did not resolve to 935 unique symbols.")

    locked = pd.read_csv(
        ROOT / "results/analysis/candidate_lock/candidate_set_v2_high_confidence.tsv", sep="\t"
    )
    locked["gene_symbol"] = locked["gene_symbol"].str.upper()
    locked["analysis_set"] = "Discovery robust"
    locked["analysis_denominator"] = len(locked)
    locked["discovery_effect"] = locked["pooled_log2FC"]

    cross = pd.read_csv(
        ROOT
        / "results/validation/cbc_external_validation_v1/tcga_icgc_identifier_harmonized_sensitivity.tsv",
        sep="\t",
    )
    joint = cross[(cross["evidence_tier"] == "discovery_robust") & as_bool(cross["jointly_estimable"])].copy()
    joint["gene_symbol"] = joint["gene_symbol"].str.upper()
    joint["analysis_set"] = "Jointly estimable in TCGA and ICGC"
    joint["analysis_denominator"] = len(joint)
    joint["discovery_effect"] = joint["discovery_log2FC"]

    detail_parts = []
    summary_rows = []
    for frame in (locked, joint):
        merged = frame[["analysis_set", "analysis_denominator", "gene_symbol", "discovery_effect"]].merge(
            allain, on="gene_symbol", how="inner", validate="one_to_one"
        )
        merged["current_direction"] = merged["discovery_effect"].map(
            lambda value: "Up in HCC" if value > 0 else "Down in HCC"
        )
        merged["direction_agreement"] = merged["current_direction"] == merged["allain_direction"]
        detail_parts.append(merged)
        n_overlap = len(merged)
        n_agree = int(merged["direction_agreement"].sum())
        summary_rows.append(
            {
                "analysis_set": frame["analysis_set"].iloc[0],
                "analysis_denominator": int(frame["analysis_denominator"].iloc[0]),
                "allain_signature_denominator": len(allain),
                "overlap_genes": n_overlap,
                "overlap_fraction_of_analysis_set": n_overlap / int(frame["analysis_denominator"].iloc[0]),
                "direction_agreement_genes": n_agree,
                "direction_agreement_fraction_among_overlap": n_agree / n_overlap if n_overlap else float("nan"),
                "interpretation": "descriptive prior-signature comparison; not an independent validation test",
            }
        )

    detail = pd.concat(detail_parts, ignore_index=True)
    summary = pd.DataFrame(summary_rows)
    allain.to_csv(OUT / "allain2016_signature_935.tsv", sep="\t", index=False)
    detail.to_csv(OUT / "allain2016_overlap_gene_level.tsv", sep="\t", index=False)
    summary.to_csv(OUT / "allain2016_overlap_summary.tsv", sep="\t", index=False)
    return summary, detail


def sample_identifier_audit() -> pd.DataFrame:
    discovery = pd.read_csv(ROOT / "00_protocol/sample_manifest_v1.0.tsv", sep="\t")
    discovery = discovery[discovery["status"].eq("include")].copy()
    tcga = pd.read_csv(ROOT / "results/validation/tcga_lihc/sample_manifest.tsv", sep="\t")
    tcga = tcga[as_bool(tcga["used_in_paired_primary"])].copy()
    icgc = pd.read_csv(ROOT / "results/validation/icgc_liri_jp/sample_manifest.tsv", sep="\t")
    icgc = icgc[as_bool(icgc["used_in_paired_primary"])].copy()

    sets = {
        "GEO discovery": set(discovery["gsm"].dropna().astype(str)),
        "TCGA-LIHC external": set(tcga["sample_barcode"].dropna().astype(str)),
        "ICGC-LIRI-JP external": set(icgc["icgc_sample_id"].dropna().astype(str)),
    }
    rows = []
    for left, right in (("GEO discovery", "TCGA-LIHC external"), ("GEO discovery", "ICGC-LIRI-JP external")):
        overlap = sorted(sets[left] & sets[right])
        rows.append(
            {
                "discovery_resource": left,
                "external_resource": right,
                "discovery_public_sample_ids": len(sets[left]),
                "external_public_sample_ids": len(sets[right]),
                "identical_public_sample_ids": len(overlap),
                "overlapping_ids": ";".join(overlap),
                "audit_scope": "exact released public sample identifiers",
                "limitation": "de-identified cross-repository patient identity cannot be tested from these identifier namespaces",
            }
        )
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "discovery_external_sample_identifier_overlap_audit.tsv", sep="\t", index=False)
    return result


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary, _ = allain_overlap()
    sample = sample_identifier_audit()

    source = ROOT / "data/external/allain2016/allain2016_supplementary_table_2.xlsx"
    manifest = pd.DataFrame(
        [
            {
                "resource": "Allain et al. Supplementary Table 2",
                "source_url": "https://ndownloader.figshare.com/files/39858996",
                "dataset_doi": "10.1158/0008-5472.22412988.v1",
                "local_path": str(source.relative_to(ROOT)),
                "sha256": sha256(source),
                "official_md5": "d6591136186c7a4fb075b0f1c19c6ffc",
            }
        ]
    )
    manifest.to_csv(OUT / "data_manifest.tsv", sep="\t", index=False)

    robust = summary.iloc[0]
    joint = summary.iloc[1]
    report = f"""# CBC 定点修订审计

## Allain 2016 定量对照

- 官方补充表解析出 935 个唯一基因。
- 173 个 discovery-robust 基因中，{int(robust.overlap_genes)} 个与该签名重叠，重叠基因方向一致 {int(robust.direction_agreement_genes)}/{int(robust.overlap_genes)}。
- TCGA/ICGC 联合可估计的 167 个基因中，{int(joint.overlap_genes)} 个重叠，方向一致 {int(joint.direction_agreement_genes)}/{int(joint.overlap_genes)}。
- 这是既往签名的描述性对照，不是第三个独立验证；两个分析的公开微阵列来源可能存在研究层面交集，本次不用重叠比例做显著性或独立复现主张。

## 发现队列与外部队列的样本标识审计

- GEO 的 GSM 标识与 TCGA-LIHC 的 sample barcode 无完全相同项。
- GEO 的 GSM 标识与 ICGC-LIRI-JP 的 sample ID 无完全相同项。
- 这只能确认已释放的公共标识符不重叠；不同数据库的去标识患者身份无法跨库核验。

## 状态

审计只新增描述性输出，未改动冻结候选集、TCGA/ICGC 模型或多重校正结果。
"""
    (OUT / "targeted_revision_audit_zh.md").write_text(report, encoding="utf-8")

    checksums = []
    for path in sorted(OUT.glob("*")):
        if path.name != "checksums.sha256" and path.is_file():
            checksums.append(f"{sha256(path)}  {path.name}")
    (OUT / "checksums.sha256").write_text("\n".join(checksums) + "\n", encoding="utf-8")

    print(summary.to_string(index=False))
    print(sample.to_string(index=False))


if __name__ == "__main__":
    main()
