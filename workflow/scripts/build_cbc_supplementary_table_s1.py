#!/usr/bin/env python3
"""Build CBC Supplementary Table S1: gene-level results for all 446 discovery hits.

Joins the locked discovery tier table with the TCGA-LIHC and identifier-harmonized
ICGC-LIRI-JP validation tables from the CBC source-data folder. No statistic is
recomputed; the script only reshapes and checks that the counts reported in the
manuscript (Table 1) are reproduced from the gene-level rows.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
# Standalone public release keeps outputs under manuscript/; the workspace under results/.
OUTDIR = ROOT / "manuscript" if (ROOT / "manuscript" / "source_data").is_dir() else ROOT / "results/manuscript_computational_biology_chemistry"
SRC = OUTDIR / "source_data"
OUT_XLSX = OUTDIR / "supplementary_table_S1_gene_level_CBC.xlsx"
OUT_TSV = OUTDIR / "supplementary_table_S1_gene_level_CBC.tsv"

TIER = {"high_confidence": "discovery robust", "intermediate_confidence": "intermediate", "exploratory_only": "exploratory"}
TIER_ORDER = {"discovery robust": 0, "intermediate": 1, "exploratory": 2}

DISCOVERY = {
    "gene_id": "entrez_id", "gene_symbol": "gene_symbol", "tier": "discovery_tier",
    "exploratory_reasons": "exploratory_reasons", "k": "discovery_n_cohorts", "total_pairs": "discovery_n_pairs",
    "pooled_log2FC": "discovery_log2FC", "CI_low": "discovery_CI_low", "CI_high": "discovery_CI_high",
    "prediction_low": "discovery_PI_low", "prediction_high": "discovery_PI_high", "p_value": "discovery_P",
    "FDR": "discovery_FDR", "I2": "discovery_I2", "tau2": "discovery_tau2",
    "direction_consistency": "discovery_direction_consistency",
}
EXTERNAL = {
    "estimable": "estimable", "estimation_note": "not_estimable_reason", "validation_log2FC": "log2FC",
    "validation_CI_low": "CI_low", "validation_CI_high": "CI_high", "nominal_P": "P",
    "candidate_family_FDR": "candidate_family_FDR", "genomewide_FDR": "genomewide_FDR",
    "same_direction": "same_direction", "candidate_family_directional": "candidate_family_FDR_lt_0.05_same_direction",
    "genomewide_directional": "genomewide_FDR_lt_0.05_same_direction",
}

COLUMN_NOTES = [
    ("entrez_id", "NCBI Entrez Gene identifier used throughout the analysis."),
    ("gene_symbol", "Current HGNC symbol in the locked discovery annotation."),
    ("discovery_tier", "Locked discovery tier: discovery robust (173), intermediate (119) or exploratory (154). Tiers were fixed before external validation."),
    ("exploratory_reasons", "Locked rule(s) that placed a gene in a lower tier; blank for discovery robust genes."),
    ("discovery_n_cohorts / discovery_n_pairs", "GEO cohorts and matched pairs contributing to the random-effects estimate."),
    ("discovery_log2FC, _CI_low/_CI_high, _PI_low/_PI_high", "Pooled tumour-minus-non-tumour log2 fold change (REML, Hartung-Knapp), 95% confidence and prediction intervals."),
    ("discovery_P / discovery_FDR", "Hartung-Knapp P value and BH FDR across 17 900 fitted genes."),
    ("discovery_I2 / discovery_tau2", "Between-cohort heterogeneity."),
    ("discovery_direction_consistency", "Fraction of contributing cohorts with the pooled direction."),
    ("tcga_* / icgc_*", "Paired limma-voom results in TCGA-LIHC (50 pairs) and ICGC-LIRI-JP (199 pairs). Only the 292 discovery robust and intermediate genes were tested externally; exploratory genes are blank by design."),
    ("*_estimable / *_not_estimable_reason", "FALSE when the gene was absent from the released matrix or removed by the unchanged cohort-wide filterByExpr rule. Non-estimable genes are not counted as failed replication."),
    ("*_candidate_family_FDR", "BH FDR within the locked evidence tier among estimable genes (primary locked-list replication criterion)."),
    ("*_genomewide_FDR", "BH FDR across the complete filtered transcriptome (15 657 TCGA genes; 19 353 ICGC genes); post hoc stricter sensitivity analysis."),
    ("*_same_direction / *_FDR_lt_0.05_same_direction", "External sign equals discovery sign; and additionally the stated FDR <0.05."),
    ("icgc_matrix_symbol / icgc_mapping_status / icgc_mapping_basis", "Symbol used in the frozen ICGC Release 28 matrix and how it was linked (post hoc, effect-blinded identifier harmonization)."),
    ("dual_*", "Criterion met in both external cohorts; defined only for genes jointly estimable in both."),
]


def read(name: str) -> pd.DataFrame:
    return pd.read_csv(SRC / name, sep="\t")


def external(df: pd.DataFrame, prefix: str, extra: list[str]) -> pd.DataFrame:
    keep = df[["gene_id", *EXTERNAL, *extra]].rename(columns=EXTERNAL)
    return keep.rename(columns={c: f"{prefix}_{c}" for c in keep.columns if c != "gene_id"})


def main() -> None:
    disc = read("candidate_downgrade_audit.tsv")
    disc["tier"] = disc["confidence_status"].map(TIER)
    assert disc["tier"].notna().all() and len(disc) == 446
    disc = disc[list(DISCOVERY)].rename(columns=DISCOVERY)

    tcga = external(read("tcga_locked_candidate_multiplicity.tsv"), "tcga", [])
    icgc = external(read("icgc_locked_candidate_identifier_harmonized_sensitivity.tsv"), "icgc",
                    ["icgc_matrix_symbol", "mapping_status", "mapping_basis"])
    icgc = icgc.rename(columns={"icgc_icgc_matrix_symbol": "icgc_matrix_symbol",
                                "icgc_mapping_status": "icgc_mapping_status", "icgc_mapping_basis": "icgc_mapping_basis"})
    assert len(tcga) == len(icgc) == 292

    out = disc.merge(tcga.rename(columns={"gene_id": "entrez_id"}), on="entrez_id", how="left", validate="1:1")
    out = out.merge(icgc.rename(columns={"gene_id": "entrez_id"}), on="entrez_id", how="left", validate="1:1")
    tested = out["tcga_estimable"].notna()
    assert tested.sum() == 292 and set(out.loc[~tested, "discovery_tier"]) == {"exploratory"}

    joint = (out["tcga_estimable"] == True) & (out["icgc_estimable"] == True)  # noqa: E712
    out["dual_jointly_estimable"] = joint.where(tested)
    for c in ("same_direction", "candidate_family_FDR_lt_0.05_same_direction", "genomewide_FDR_lt_0.05_same_direction"):
        out[f"dual_{c}"] = (joint & (out[f"tcga_{c}"] == True) & (out[f"icgc_{c}"] == True)).where(joint)  # noqa: E712

    # Cross-check against the independently written joint table used for Fig. 2.
    ref = read("tcga_icgc_identifier_harmonized_sensitivity.tsv").set_index("gene_id")
    chk = out.set_index("entrez_id").loc[ref.index]
    for a, b in (("dual_jointly_estimable", "jointly_estimable"), ("dual_same_direction", "dual_same_direction"),
                 ("dual_candidate_family_FDR_lt_0.05_same_direction", "dual_candidate_family_directional"),
                 ("dual_genomewide_FDR_lt_0.05_same_direction", "dual_genomewide_directional")):
        assert (chk[a].fillna(False).astype(bool) == ref[b].fillna(False).astype(bool)).all(), a

    # Counts reported in manuscript Table 1 (discovery robust tier).
    r = out[out["discovery_tier"] == "discovery robust"]
    n = lambda s: int((s == True).sum())  # noqa: E731,E712
    got = {
        "tcga_estimable": n(r.tcga_estimable), "tcga_dir": n(r.tcga_same_direction),
        "tcga_cf": n(r["tcga_candidate_family_FDR_lt_0.05_same_direction"]), "tcga_gw": n(r["tcga_genomewide_FDR_lt_0.05_same_direction"]),
        "icgc_estimable": n(r.icgc_estimable), "icgc_dir": n(r.icgc_same_direction),
        "icgc_cf": n(r["icgc_candidate_family_FDR_lt_0.05_same_direction"]), "icgc_gw": n(r["icgc_genomewide_FDR_lt_0.05_same_direction"]),
        "joint": n(r.dual_jointly_estimable), "dual_dir": n(r.dual_same_direction),
        "dual_cf": n(r["dual_candidate_family_FDR_lt_0.05_same_direction"]), "dual_gw": n(r["dual_genomewide_FDR_lt_0.05_same_direction"]),
    }
    want = dict(tcga_estimable=168, tcga_dir=167, tcga_cf=163, tcga_gw=161, icgc_estimable=171, icgc_dir=171,
                icgc_cf=170, icgc_gw=170, joint=167, dual_dir=166, dual_cf=161, dual_gw=159)
    assert got == want, (got, want)

    out["_o"] = out["discovery_tier"].map(TIER_ORDER)
    out = out.sort_values(["_o", "discovery_FDR", "entrez_id"]).drop(columns="_o")
    out.to_csv(OUT_TSV, sep="\t", index=False, na_rep="")

    readme = pd.DataFrame(
        [("Supplementary Table S1", "Gene-level discovery and external validation results for all 446 discovery hits."),
         ("Rows", "446 genes: 173 discovery robust, 119 intermediate, 154 exploratory; sorted by tier then discovery FDR."),
         ("Direction", "All effects are tumour minus matched non-tumour tissue (log2 scale)."),
         ("Blank cells", "Not applicable: exploratory genes were not tested externally; FDR and effect columns are blank for non-estimable genes."),
         ("", "")] + [(c, d) for c, d in COLUMN_NOTES],
        columns=["Field", "Description"],
    )
    with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as xw:
        readme.to_excel(xw, sheet_name="README", index=False)
        out.to_excel(xw, sheet_name="Table_S1", index=False, freeze_panes=(1, 2))
        xw.sheets["README"].column_dimensions["A"].width = 52
        xw.sheets["README"].column_dimensions["B"].width = 120
    print(OUT_XLSX, OUT_TSV, len(out), got)


if __name__ == "__main__":
    main()
