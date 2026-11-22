#!/usr/bin/env python3
"""Rebuild Fig.5 cross-endpoint evidence map with discrete status colours.

Preserves the original gene selection rule, row order and panel-b statistics.
Improves layout, categorical colour legend and export readability.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "manuscript/figures"
MANUSCRIPT_FIG = OUT
SRC = ROOT / "manuscript/source_data"
MM = 1 / 25.4

# Discrete status codes (categorical map, not continuous values)
COVERED = 0
MISSING = 1
QUERY_FAIL = 2
BELOW_THRESH = 3

COLOR_COVERED = "#2878B5"
COLOR_MISSING = "#E6E6E6"
COLOR_QUERY_FAIL = "#E69F00"
COLOR_BELOW = "#9EC8E6"  # only used when source distinguishes non-pass but estimable

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans", "sans-serif"],
        "font.size": 7.0,
        "axes.titlesize": 8.0,
        "axes.labelsize": 7.0,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "legend.fontsize": 6.2,
        "axes.linewidth": 0.7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    }
)


def read_tsv(rel: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / rel, sep="\t", low_memory=False)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.08, 1.04, label, transform=ax.transAxes, ha="right", va="bottom", fontsize=9, fontweight="bold", clip_on=False)


def status_bulk_pass(series: pd.Series) -> np.ndarray:
    """validated_primary → covered; not_estimable → missing; other estimable non-pass → below threshold."""
    out = np.full(len(series), MISSING, dtype=int)
    s = series.astype(str)
    out[s.eq("validated_primary")] = COVERED
    out[s.isin(["direction_consistent", "validated_secondary", "direction_inconsistent"])] = BELOW_THRESH
    out[s.eq("not_estimable")] = MISSING
    return out


def status_bool_covered(series: pd.Series) -> np.ndarray:
    out = np.full(len(series), MISSING, dtype=int)
    flag = series.astype(str).str.upper().isin(["TRUE", "1", "YES"])
    out[flag.to_numpy()] = COVERED
    return out


def status_both_dir(series: pd.Series) -> np.ndarray:
    out = np.full(len(series), MISSING, dtype=int)
    s = series.astype(str).str.upper()
    out[s.eq("TRUE")] = COVERED
    out[s.eq("FALSE")] = BELOW_THRESH
    # NA / empty remain MISSING
    return out


def status_chembl(series: pd.Series) -> np.ndarray:
    out = np.full(len(series), MISSING, dtype=int)
    s = series.astype(str)
    out[s.eq("mapped")] = COVERED
    out[s.eq("no_target")] = MISSING
    out[s.eq("search_failed")] = QUERY_FAIL
    return out


def status_lincs(series: pd.Series) -> np.ndarray:
    """query_* = covered (in SigCom query set); not_in_signature = missing/not covered (not a query failure)."""
    out = np.full(len(series), MISSING, dtype=int)
    s = series.astype(str)
    out[s.str.startswith("query_")] = COVERED
    out[s.eq("not_in_signature")] = MISSING
    return out


def build_selected() -> pd.DataFrame:
    cross = read_tsv("results/validation/phase3_synthesis/cross_endpoint_evidence_table.tsv")
    high = cross[cross["tier"] == "high_confidence"].copy()
    bulk_v2 = read_tsv("results/validation/rnaseq_normalized_v2/tcga_icgc_high_cross_validation_TMM_voom.tsv")
    high = high.merge(
        bulk_v2[["gene_symbol", "tcga_status", "icgc_status", "dual_same_direction"]],
        on="gene_symbol",
        how="left",
        validate="one_to_one",
    )
    high["tcga_validation_status"] = high["tcga_status"]
    high["icgc_validation_status"] = high["icgc_status"]
    high["both_bulk_same_direction"] = high["dual_same_direction"]
    high["both_valid"] = high["tcga_validation_status"].eq("validated_primary") & high["icgc_validation_status"].eq(
        "validated_primary"
    )
    high["score"] = (
        high["both_valid"].astype(int) * 5
        + high["both_bulk_same_direction"].astype(str).str.upper().eq("TRUE").astype(int) * 2
        + high["scrna_detected"].astype(str).str.upper().eq("TRUE").astype(int)
        + high["chembl_target_status"].eq("mapped").astype(int)
    )
    selected = high.sort_values(["score", "discovery_log2FC"], ascending=[False, True]).head(28).copy()
    return high, selected


def fig5() -> None:
    high, selected = build_selected()

    col_specs = [
        ("TCGA", status_bulk_pass(selected["tcga_validation_status"]), "TCGA primary pass"),
        ("ICGC", status_bulk_pass(selected["icgc_validation_status"]), "ICGC primary pass"),
        ("Both dir.", status_both_dir(selected["both_bulk_same_direction"]), "Same direction in both"),
        ("scRNA", status_bool_covered(selected["scrna_detected"]), "Detected in scRNA"),
        ("HPA", status_bool_covered(selected["hpa_matched"]), "HPA annotation match"),
        ("ChEMBL", status_chembl(selected["chembl_target_status"]), "ChEMBL target mapped"),
        ("LINCS", status_lincs(selected["lincs_query_membership"]), "In LINCS query set"),
        ("DepMap", status_bool_covered(selected["depmap_present"]), "DepMap score present"),
    ]
    matrix = np.column_stack([spec[1] for spec in col_specs])
    col_labels = [spec[0] for spec in col_specs]

    # Snapshot for audit (preserve raw source fields + coded status)
    snap_rows = []
    for i, gene in enumerate(selected["gene_symbol"].tolist()):
        row = {
            "row_order": i + 1,
            "gene_symbol": gene,
            "tcga_validation_status": selected.iloc[i]["tcga_validation_status"],
            "icgc_validation_status": selected.iloc[i]["icgc_validation_status"],
            "both_bulk_same_direction": selected.iloc[i]["both_bulk_same_direction"],
            "scrna_detected": selected.iloc[i]["scrna_detected"],
            "hpa_matched": selected.iloc[i]["hpa_matched"],
            "chembl_target_status": selected.iloc[i]["chembl_target_status"],
            "lincs_query_membership": selected.iloc[i]["lincs_query_membership"],
            "depmap_present": selected.iloc[i]["depmap_present"],
        }
        for j, (lab, codes, _) in enumerate(col_specs):
            code = int(codes[i])
            row[f"status_{lab}"] = {COVERED: "covered", MISSING: "missing", QUERY_FAIL: "query_fail", BELOW_THRESH: "below_threshold"}[code]
        snap_rows.append(row)
    SRC.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(snap_rows).to_csv(SRC / "Fig5a_evidence_status_snapshot.tsv", sep="\t", index=False)

    used = set(np.unique(matrix).tolist())
    cmap_colors = [COLOR_COVERED, COLOR_MISSING, COLOR_QUERY_FAIL, COLOR_BELOW]
    cmap = ListedColormap(cmap_colors)

    # Wider/taller canvas for 28 gene labels + legend
    fig_w, fig_h = 190 * MM, 155 * MM
    fig, axs = plt.subplots(
        1,
        2,
        figsize=(fig_w, fig_h),
        gridspec_kw={"width_ratios": [1.55, 1.0], "wspace": 0.28},
    )
    fig.subplots_adjust(left=0.11, right=0.98, top=0.88, bottom=0.16)

    ax = axs[0]
    ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=cmap, vmin=-0.5, vmax=3.5)
    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_yticks(np.arange(len(selected)))
    ax.set_yticklabels(selected["gene_symbol"].tolist())
    ax.set_title("Cross-endpoint evidence map for selected high-confidence genes", loc="left", pad=8)
    ax.tick_params(length=0, pad=2)
    # Very light separators only (avoid dense white grid dominating tiles)
    ax.set_xticks(np.arange(matrix.shape[1] + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(matrix.shape[0] + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="#FBFBFB", linestyle="-", linewidth=0.2)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.45)
        spine.set_color("#A8A8A8")

    legend_items = [
        (COVERED, COLOR_COVERED, "Covered / criterion met"),
        (MISSING, COLOR_MISSING, "Missing / not estimable / not in resource"),
        (QUERY_FAIL, COLOR_QUERY_FAIL, "Query failure (technical)"),
        (BELOW_THRESH, COLOR_BELOW, "Tested, below primary threshold"),
    ]
    handles = [Patch(facecolor=c, edgecolor="#B0B0B0", linewidth=0.4, label=lab) for code, c, lab in legend_items if code in used or code in (COVERED, MISSING, QUERY_FAIL)]
    # Always show the three recommended colours; add below-threshold only if used
    handles = []
    for code, c, lab in legend_items:
        if code == BELOW_THRESH and code not in used:
            continue
        handles.append(Patch(facecolor=c, edgecolor="#B0B0B0", linewidth=0.4, label=lab))
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, frameon=False, handlelength=1.2, columnspacing=1.0)
    panel_label(ax, "a")

    ax = axs[1]
    metrics = [
        ("TCGA validated", int((high["tcga_validation_status"] == "validated_primary").sum()), 173),
        ("ICGC validated", int((high["icgc_validation_status"] == "validated_primary").sum()), 173),
        ("Both validated", int(high["both_valid"].sum()), 173),
        ("Both same direction", int(high["both_bulk_same_direction"].astype(str).str.upper().eq("TRUE").sum()), 173),
        ("scRNA detected", int(high["scrna_detected"].astype(str).str.upper().eq("TRUE").sum()), 173),
        ("HPA matched", int(high["hpa_matched"].astype(str).str.upper().eq("TRUE").sum()), 173),
        ("ChEMBL mapped", int(high["chembl_target_status"].eq("mapped").sum()), 173),
        ("DepMap liver-selective", int(high["depmap_liver_selective_essential_fdr05"].astype(str).str.upper().eq("TRUE").sum()), 173),
    ]
    labels = [m[0] for m in metrics][::-1]
    nums = np.array([m[1] for m in metrics][::-1])
    den = np.array([m[2] for m in metrics][::-1])
    yy = np.arange(len(labels))
    ax.barh(yy, den, color=COLOR_MISSING, height=0.68)
    bar_colors = ["#C95756" if "liver-selective" in x else COLOR_COVERED for x in labels]
    ax.barh(yy, nums, color=bar_colors, height=0.68)
    ax.set_yticks(yy, labels)
    ax.tick_params(axis="y", pad=4)
    ax.set_xlim(0, 185)
    ax.set_xlabel("High-confidence candidates (denominator 173)")
    ax.set_title("Evidence coverage and negative endpoints", loc="left", pad=8)
    for y, n in zip(yy, nums):
        if n > 12:
            ax.text(n - 3, y, str(int(n)), va="center", ha="right", fontsize=6.0, color="white", fontweight="bold")
        else:
            ax.text(3, y, str(int(n)), va="center", ha="left", fontsize=6.0, color="#333333")
    ax.grid(False)
    ax.legend(
        handles=[
            Patch(facecolor=COLOR_COVERED, edgecolor="#B0B0B0", linewidth=0.4, label="Count meeting endpoint"),
            Patch(facecolor="#C95756", edgecolor="#B0B0B0", linewidth=0.4, label="Liver-selective = 0"),
            Patch(facecolor=COLOR_MISSING, edgecolor="#B0B0B0", linewidth=0.4, label="Denominator (173)"),
        ],
        loc="lower right",
        fontsize=5.6,
    )
    panel_label(ax, "b")

    OUT.mkdir(parents=True, exist_ok=True)
    MANUSCRIPT_FIG.mkdir(parents=True, exist_ok=True)
    stem = "Fig5_cross_endpoint_synthesis"
    for folder in (OUT, MANUSCRIPT_FIG):
        fig.savefig(folder / f"{stem}.pdf", bbox_inches="tight")
        fig.savefig(folder / f"{stem}.svg", bbox_inches="tight")
        fig.savefig(folder / f"{stem}.png", dpi=300, bbox_inches="tight")
        fig.savefig(folder / f"{stem}.tiff", dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)

    # Status summary for delivery note
    summary = {
        "n_genes": int(matrix.shape[0]),
        "n_cols": int(matrix.shape[1]),
        "used_codes": sorted(used),
        "lincs_missing_n": int((matrix[:, 6] == MISSING).sum()),
        "depmap_missing_n": int((matrix[:, 7] == MISSING).sum()),
        "query_fail_n": int((matrix == QUERY_FAIL).sum()),
        "below_thresh_n": int((matrix == BELOW_THRESH).sum()),
        "gene_order": selected["gene_symbol"].tolist(),
    }
    pd.Series(summary).to_json(SRC / "Fig5a_rebuild_summary.json")
    print("Wrote Fig5 to", OUT)
    print(summary)


if __name__ == "__main__":
    fig5()
