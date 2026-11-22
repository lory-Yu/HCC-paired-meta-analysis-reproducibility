#!/usr/bin/env python3
"""Secondary single-cell donor-level count-sum pseudobulk sensitivity.

Design decision (see feasibility tables):
- Complete paired donor design: NOT supported (only donor 0618 has both tissues).
- Unpaired donor-independent design: POSSIBLE after collapsing technical replicates
  (6 tumour donors vs 3 adjacent donors).
- Partial pairing: only 1 paired donor; retained for description, not tested alone.

This script aggregates raw 10x counts for locked high-confidence genes within each
donor × tissue × broad cell type, then tests tumour vs adjacent donors inside each
gene's dominant cell class by two-sided Mann–Whitney U on log2(CPM+1).
Prespecified filters are fixed before looking at P values.
"""
from __future__ import annotations

import gzip
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/validation/medical_oncology_revision_v1"
DATA = ROOT / "data/external/gse202642"
GSE = ROOT / "results/validation/gse202642"

MIN_CELLS_PER_DONOR_UNIT = 20
MIN_DONORS_PER_GROUP = 3
TEST_FAMILY = "high_confidence_dominant_class_donor_tests"


def donor_id(patient_label: str) -> str:
    if patient_label.startswith(("0618", "0629", "0706", "0707", "0727", "0831")):
        return patient_label.split("_")[0]
    if patient_label in {"X_1b", "B_2b"}:
        return patient_label.split("_")[0]
    return patient_label


def main() -> None:
    meta = pd.read_csv(GSE / "sample_manifest.tsv", sep="\t")
    meta["donor_id"] = meta["patient_label"].map(donor_id)
    cm = pd.read_csv(
        GSE / "cell_metadata.tsv.gz",
        sep="\t",
        usecols=["cell_barcode", "barcode_suffix", "gsm_id", "patient_label", "tissue_class", "cell_type_broad"],
    )
    cm["donor_id"] = cm["patient_label"].map(donor_id)

    # cell counts inventory
    inv = (
        cm.groupby(["donor_id", "tissue_class", "patient_label", "gsm_id", "cell_type_broad"], as_index=False)
        .size()
        .rename(columns={"size": "n_cells"})
    )
    inv.to_csv(OUT / "singlecell_donor_tissue_celltype_cellcounts.tsv", sep="\t", index=False)

    design = pd.DataFrame(
        [
            {"item": "n_tissue_samples", "value": str(meta.shape[0])},
            {"item": "n_tumour_donors_after_collapsing_tech_reps", "value": str(meta.loc[meta.tissue_class.eq("tumor"), "donor_id"].nunique())},
            {"item": "n_adjacent_donors_after_collapsing_tech_reps", "value": str(meta.loc[meta.tissue_class.eq("adjacent"), "donor_id"].nunique())},
            {"item": "n_paired_donors_both_tissues", "value": str(int(meta.groupby("donor_id")["tissue_class"].nunique().ge(2).sum()))},
            {"item": "paired_donor_ids", "value": ",".join(sorted(meta.groupby("donor_id").filter(lambda d: d.tissue_class.nunique() == 2)["donor_id"].unique()))},
            {"item": "tech_replicate_note", "value": "0706 has 2 tumour samples; 0618 has 2 adjacent samples; collapsed by summing counts within donor×tissue×celltype"},
            {"item": "paired_design_supported", "value": "FALSE"},
            {"item": "unpaired_donor_design_supported", "value": "TRUE_WITH_LIMITATIONS"},
            {"item": "min_cells_per_donor_unit", "value": str(MIN_CELLS_PER_DONOR_UNIT)},
            {"item": "min_donors_per_group", "value": str(MIN_DONORS_PER_GROUP)},
            {"item": "model", "value": "MannWhitneyU_on_donor_log2CPM_plus1_within_dominant_celltype; edgeR/limma unavailable in runtime"},
            {"item": "analysis_label", "value": "secondary_donor_unpaired_count_sum_sensitivity"},
        ]
    )
    design.to_csv(OUT / "singlecell_donor_design_decision.tsv", sep="\t", index=False)

    dom = pd.read_csv(GSE / "candidate_dominant_celltype.tsv", sep="\t")
    high = pd.read_csv(ROOT / "results/analysis/candidate_lock/candidate_set_v2_high_confidence.tsv", sep="\t")
    genes = high[["gene_id", "gene_symbol", "pooled_log2FC"]].merge(
        dom.loc[dom.evidence_tier.eq("high_confidence"), ["gene_symbol", "dominant_cell_type", "n_cells"]],
        on="gene_symbol",
        how="left",
    )

    # feature index by symbol (first match)
    with gzip.open(DATA / "GSE202642_features.tsv.gz", "rt") as fh:
        features = [line.rstrip("\n").split("\t") for line in fh]
    symbol_to_row = {}
    for i, row in enumerate(features):
        sym = row[1]
        if sym not in symbol_to_row:
            symbol_to_row[sym] = i
    genes["feature_row"] = genes["gene_symbol"].map(symbol_to_row)
    missing = genes[genes.feature_row.isna()]
    print("genes missing from features", missing.gene_symbol.tolist())

    # barcode -> unit
    # cell_barcode in metadata already includes -suffix
    with gzip.open(DATA / "GSE202642_barcodes.tsv.gz", "rt") as fh:
        barcodes = [line.strip() for line in fh]
    bc_to_j = {b: j for j, b in enumerate(barcodes)}

    unit_info = {}
    for r in cm.itertuples(index=False):
        key = (r.donor_id, r.tissue_class, r.cell_type_broad)
        unit_info.setdefault(key, {"barcodes": [], "n_cells": 0, "gsm_ids": set(), "patient_labels": set()})
        unit_info[key]["barcodes"].append(r.cell_barcode)
        unit_info[key]["n_cells"] += 1
        unit_info[key]["gsm_ids"].add(r.gsm_id)
        unit_info[key]["patient_labels"].add(r.patient_label)

    # Map barcode column index -> unit key for QC cells only
    j_to_unit = {}
    for key, info in unit_info.items():
        for bc in info["barcodes"]:
            j = bc_to_j.get(bc)
            if j is not None:
                j_to_unit[j] = key

    wanted_rows = sorted({int(r) for r in genes.feature_row.dropna().astype(int)})
    wanted_set = set(wanted_rows)
    # accumulate count sums: (gene_row, unit_key) -> count
    sums = defaultdict(float)
    libsize = defaultdict(float)  # total counts per unit across ALL genes for CPM? expensive.
    # For CPM within unit we need total counts per cell aggregated. Approximate using nCount from metadata.
    ncount = pd.read_csv(GSE / "cell_metadata.tsv.gz", sep="\t", usecols=["cell_barcode", "nCount"])
    ncount_map = dict(zip(ncount.cell_barcode, ncount.nCount))
    unit_lib = defaultdict(float)
    for key, info in unit_info.items():
        unit_lib[key] = float(sum(ncount_map.get(bc, 0.0) for bc in info["barcodes"]))

    print("Streaming MTX for", len(wanted_rows), "candidate feature rows ...")
    path = DATA / "GSE202642_matrix.mtx.gz"
    with gzip.open(path, "rt") as fh:
        # skip header comments
        while True:
            line = fh.readline()
            if not line.startswith("%"):
                break
        dims = line.split()
        n_rows, n_cols, n_entries = map(int, dims)
        print("MTX dims", n_rows, n_cols, n_entries)
        for k, line in enumerate(fh, 1):
            r, c, v = line.split()
            r = int(r) - 1
            if r not in wanted_set:
                continue
            c = int(c) - 1
            unit = j_to_unit.get(c)
            if unit is None:
                continue
            sums[(r, unit)] += float(v)
            if k % 20_000_000 == 0:
                print("scanned", k)

    # Build donor-level table for dominant-class tests
    records = []
    for g in genes.itertuples(index=False):
        if pd.isna(g.feature_row) or pd.isna(g.dominant_cell_type):
            records.append(
                {
                    "gene_symbol": g.gene_symbol,
                    "gene_id": g.gene_id,
                    "discovery_log2FC": g.pooled_log2FC,
                    "dominant_cell_type": g.dominant_cell_type,
                    "status": "missing_feature_or_dominant_class",
                    "n_tumor_donors": 0,
                    "n_adjacent_donors": 0,
                    "mean_tumor_log2cpm": np.nan,
                    "mean_adjacent_log2cpm": np.nan,
                    "log2FC_tumor_minus_adjacent": np.nan,
                    "mannwhitney_p": np.nan,
                    "same_direction": np.nan,
                }
            )
            continue
        row = int(g.feature_row)
        ctype = g.dominant_cell_type
        tumor_vals = []
        adj_vals = []
        detail = []
        for tissue, bucket in (("tumor", tumor_vals), ("adjacent", adj_vals)):
            donors = sorted({d for (d, t, ct) in unit_info if t == tissue and ct == ctype})
            for d in donors:
                key = (d, tissue, ctype)
                info = unit_info[key]
                if info["n_cells"] < MIN_CELLS_PER_DONOR_UNIT:
                    continue
                count = sums.get((row, key), 0.0)
                lib = unit_lib[key]
                log2cpm = np.log2((count / lib) * 1e6 + 1.0) if lib > 0 else np.nan
                bucket.append(log2cpm)
                detail.append({"donor_id": d, "tissue_class": tissue, "n_cells": int(info["n_cells"]), "count_sum": count, "libsize": lib, "log2cpm": float(log2cpm)})
        status = "ok"
        p = np.nan
        log2fc = np.nan
        same = np.nan
        if len(tumor_vals) < MIN_DONORS_PER_GROUP or len(adj_vals) < MIN_DONORS_PER_GROUP:
            status = "insufficient_donor_units"
        else:
            # fixed test family before viewing
            stat = stats.mannwhitneyu(tumor_vals, adj_vals, alternative="two-sided")
            p = float(stat.pvalue)
            log2fc = float(np.nanmean(tumor_vals) - np.nanmean(adj_vals))
            same = bool(np.sign(log2fc) == np.sign(g.pooled_log2FC)) if log2fc != 0 else False
        records.append(
            {
                "gene_symbol": g.gene_symbol,
                "gene_id": g.gene_id,
                "discovery_log2FC": g.pooled_log2FC,
                "dominant_cell_type": ctype,
                "status": status,
                "n_tumor_donors": len(tumor_vals),
                "n_adjacent_donors": len(adj_vals),
                "mean_tumor_log2cpm": float(np.nanmean(tumor_vals)) if tumor_vals else np.nan,
                "mean_adjacent_log2cpm": float(np.nanmean(adj_vals)) if adj_vals else np.nan,
                "log2FC_tumor_minus_adjacent": log2fc,
                "mannwhitney_p": p,
                "same_direction": same,
                "donor_detail_json": json.dumps(detail, separators=(",", ":")),
            }
        )

    res = pd.DataFrame(records)
    testable = res[res.status.eq("ok")].copy()
    if len(testable):
        _, fdr, _, _ = multipletests(testable["mannwhitney_p"].to_numpy(float), method="fdr_bh")
        res.loc[testable.index, "fdr"] = fdr
    else:
        res["fdr"] = np.nan
    res["secondary_sensitivity_analysis"] = True
    res["test_family"] = TEST_FAMILY
    res.to_csv(OUT / "singlecell_donor_unpaired_countsum_pseudobulk_results.tsv", sep="\t", index=False)

    summary = pd.DataFrame(
        [
            {"metric": "high_confidence_genes", "value": int(len(genes))},
            {"metric": "donor_unpaired_testable", "value": int((res.status == "ok").sum())},
            {"metric": "same_direction", "value": int((res.same_direction == True).sum())},
            {"metric": "opposite_direction", "value": int((res.same_direction == False).sum())},
            {"metric": "FDR_lt_0.05", "value": int((pd.to_numeric(res.fdr, errors="coerce") < 0.05).fillna(False).sum())},
            {"metric": "insufficient_donor_units", "value": int((res.status == "insufficient_donor_units").sum())},
            {"metric": "paired_design_tested", "value": 0},
            {"metric": "min_cells_per_donor_unit", "value": MIN_CELLS_PER_DONOR_UNIT},
            {"metric": "min_donors_per_group", "value": MIN_DONORS_PER_GROUP},
        ]
    )
    summary.to_csv(OUT / "singlecell_donor_unpaired_countsum_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))
    print("Wrote donor unpaired count-sum sensitivity results")


if __name__ == "__main__":
    main()
