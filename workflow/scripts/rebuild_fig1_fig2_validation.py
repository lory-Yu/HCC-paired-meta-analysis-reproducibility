#!/usr/bin/env python3
"""Rebuild Fig.1 and Fig.2 with explicit evidence-tier and denominator labels.

Uses the same source tables as workflow/scripts/build_manuscript_figures.py.
Does not change discovery locks or validation statistics—only display labels and
panel titles/legends needed for denominator clarity.
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
REV = ROOT / "results/validation/medical_oncology_revision_v1"
OUT = ROOT / "manuscript/figures"
MANUSCRIPT_FIG = OUT
MM = 1 / 25.4
COLORS = {
    "navy": "#315A7D",
    "blue": "#4C93C3",
    "teal": "#3C9D9B",
    "orange": "#E58B4A",
    "red": "#C95756",
    "grey": "#B7BDC4",
    "dark": "#3F4650",
}

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans", "sans-serif"],
        "font.size": 6.5,
        "axes.titlesize": 7.0,
        "axes.labelsize": 6.8,
        "xtick.labelsize": 5.8,
        "ytick.labelsize": 5.8,
        "legend.fontsize": 5.5,
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

TIER_LABEL = {
    "high_confidence": "discovery-stage robust",
    "intermediate_confidence": "intermediate",
    "exploratory_only": "exploratory",
}


def read_tsv(rel: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / rel, sep="\t", low_memory=False)


def panel_label(ax, label: str) -> None:
    ax.text(-0.11, 1.075, label, transform=ax.transAxes, ha="right", va="top", fontsize=8, fontweight="bold", clip_on=False)


def clean_axis(ax) -> None:
    ax.grid(axis="y", color="#DDE2E7", linewidth=0.45, alpha=0.8)
    ax.set_axisbelow(True)


def export(fig, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{stem}.{ext}", dpi=300 if ext == "png" else None, bbox_inches="tight")
    # also refresh manuscript figure directory used by Word builder
    for ext in ("png", "pdf", "svg"):
        fig.savefig(MANUSCRIPT_FIG / f"{stem}.{ext}", dpi=300 if ext == "png" else None, bbox_inches="tight")
    fig.savefig(MANUSCRIPT_FIG / f"{stem}.tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)


def fig1() -> None:
    cohort = read_tsv("results/analysis/differential/cohort_differential_summary.tsv")
    meta = pd.read_csv(ROOT / "results/analysis/meta/meta_all_genes.tsv.gz", sep="\t", low_memory=False)
    audit = read_tsv("results/analysis/sensitivity/candidate_downgrade_audit.tsv")

    fig, axs = plt.subplots(2, 2, figsize=(183 * MM, 145 * MM), constrained_layout=True)
    ax = axs[0, 0]
    ax.bar(cohort["cohort"], cohort["n_pairs"], color=COLORS["navy"], width=0.65)
    for x, y in enumerate(cohort["n_pairs"]):
        ax.text(x, y - 2.0, str(int(y)), ha="center", va="top", fontsize=5.5, color="white", fontweight="bold")
    ax.set_ylabel("Matched patient pairs")
    ax.set_title("Six prespecified discovery cohorts (167 pairs)", loc="left")
    ax.tick_params(axis="x", rotation=35)
    for tick in ax.get_xticklabels():
        tick.set_rotation_mode("anchor")
        tick.set_ha("right")
    clean_axis(ax)
    ax.grid(False)
    panel_label(ax, "a")

    ax = axs[0, 1]
    x = meta["pooled_log2FC"].to_numpy(float)
    y = -np.log10(np.clip(meta["FDR"].to_numpy(float), 1e-300, 1))
    primary = meta["passes_primary_rule"].astype(str).str.upper().eq("TRUE")
    up = primary & (x > 0)
    down = primary & (x < 0)
    ax.scatter(x[~primary], y[~primary], s=3, color=COLORS["grey"], alpha=0.45, linewidth=0, rasterized=True)
    ax.scatter(x[down], y[down], s=5, color=COLORS["blue"], alpha=0.78, linewidth=0, rasterized=True, label="Primary down")
    ax.scatter(x[up], y[up], s=5, color=COLORS["red"], alpha=0.78, linewidth=0, rasterized=True, label="Primary up")
    ax.axvline(-math.log2(1.5), color=COLORS["dark"], lw=0.7, ls="--")
    ax.axvline(math.log2(1.5), color=COLORS["dark"], lw=0.7, ls="--")
    ax.axhline(-math.log10(0.05), color=COLORS["dark"], lw=0.7, ls="--")
    ax.set_xlabel("Pooled log₂ fold change (tumour − control)")
    ax.set_ylabel("−log₁₀ meta-analysis FDR")
    ax.set_title("Random-effects meta-analysis of 17,900 genes", loc="left")
    ax.legend(loc="lower left")
    panel_label(ax, "b")

    ax = axs[1, 0]
    tier_order = ["high_confidence", "intermediate_confidence", "exploratory_only"]
    tier_colors = {"high_confidence": COLORS["teal"], "intermediate_confidence": COLORS["orange"], "exploratory_only": COLORS["grey"]}
    for tier in tier_order[::-1]:
        d = audit[audit["confidence_status"] == tier]
        ax.scatter(
            d["pooled_log2FC"],
            d["I2"],
            s=11 if tier == "high_confidence" else 8,
            color=tier_colors[tier],
            alpha=0.78,
            linewidth=0,
            label=TIER_LABEL[tier],
            rasterized=True,
        )
    ax.axhline(50, color=COLORS["dark"], lw=0.7, ls="--")
    ax.axhline(75, color=COLORS["dark"], lw=0.7, ls=":")
    ax.axvline(0, color=COLORS["dark"], lw=0.6)
    ax.set_xlabel("Pooled log₂ fold change")
    ax.set_ylabel("Heterogeneity, I² (%)")
    ax.set_title("Discovery-stage tiers among 446 primary hits (n=173/119/154)", loc="left")
    ax.legend(loc="upper right", handletextpad=0.3)
    panel_label(ax, "c")

    ax = axs[1, 1]
    ax.scatter(audit["pooled_log2FC"], audit["fixed_log2FC"], s=7, color=COLORS["teal"], alpha=0.5, linewidth=0, rasterized=True, label="Fixed-effect model")
    ax.scatter(audit["pooled_log2FC"], audit["all_samples_log2FC"], s=7, color=COLORS["orange"], alpha=0.5, linewidth=0, rasterized=True, label="Sensitivity: 438-sample model")
    lim = np.nanmax(np.abs(np.r_[audit["pooled_log2FC"], audit["fixed_log2FC"], audit["all_samples_log2FC"]])) * 1.04
    ax.plot([-lim, lim], [-lim, lim], color=COLORS["dark"], lw=0.7, ls="--")
    ax.axhline(0, color=COLORS["dark"], lw=0.5)
    ax.axvline(0, color=COLORS["dark"], lw=0.5)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("Primary random-effects log₂ fold change")
    ax.set_ylabel("Sensitivity-model log₂ fold change")
    ax.set_title("Effect estimates remained directionally stable", loc="left")
    ax.legend(loc="upper left")
    panel_label(ax, "d")
    export(fig, "Fig1_discovery_and_robustness")


def fig2() -> None:
    # Effects/status: TCGA + intermediate from CBC/rnaseq tables; discovery-robust ICGC from unambiguous primary
    tcga = read_tsv("results/validation/rnaseq_normalized_v2/tcga_candidate_validation_TMM_voom.tsv")
    icgc = read_tsv("results/validation/rnaseq_normalized_v2/icgc_candidate_validation_TMM_voom.tsv")
    joint = read_tsv("results/validation/rnaseq_normalized_v2/tcga_icgc_high_cross_validation_TMM_voom.tsv")
    cross = read_tsv("results/validation/cbc_external_validation_v1/tcga_icgc_locked_candidate_cross_validation.tsv")
    counts = read_tsv("results/validation/cbc_external_validation_v1/multiplicity_counts_with_denominators.tsv")
    primary = pd.read_csv(
        REV / "external_validation_primary_unambiguous_harmonized_gene_level.tsv",
        sep="\t",
        low_memory=False,
    )

    # Save figure source snapshot used this round
    src = ROOT / "manuscript/source_data"
    src.mkdir(parents=True, exist_ok=True)
    counts.to_csv(src / "Fig2_denominator_counts_source.tsv", sep="\t", index=False)
    cross.to_csv(src / "Fig2_cross_validation_source.tsv", sep="\t", index=False)
    primary.to_csv(src / "Fig2_primary_unambiguous_source.tsv", sep="\t", index=False)

    def _flag(s: pd.Series) -> pd.Series:
        return s.astype(str).str.lower().isin(["true", "1", "yes"])

    # ICGC scatter: intermediate from original file; robust from unambiguous primary (includes +30 aliases)
    icgc_inter = icgc[icgc["evidence_tier"].eq("intermediate_confidence")].copy()
    prim_pts = primary.loc[_flag(primary["icgc_ok"]), ["gene_symbol", "discovery_log2FC", "icgc_log2FC"]].copy()
    prim_pts = prim_pts.rename(columns={"icgc_log2FC": "validation_log2FC"})
    prim_pts["evidence_tier"] = "high_confidence"
    icgc_plot = pd.concat([prim_pts, icgc_inter], ignore_index=True, sort=False)

    fig, axs = plt.subplots(2, 2, figsize=(183 * MM, 150 * MM), constrained_layout=True)
    tier_colors = {"high_confidence": COLORS["teal"], "intermediate_confidence": COLORS["orange"]}
    for ax, data, label, n_pairs, letter in [
        (axs[0, 0], tcga, "TCGA-LIHC", 50, "a"),
        (axs[0, 1], icgc_plot, "ICGC-LIRI-JP", 199, "b"),
    ]:
        finite = data["validation_log2FC"].notna()
        for tier, color in tier_colors.items():
            d = data[finite & data["evidence_tier"].eq(tier)]
            ax.scatter(
                d["discovery_log2FC"],
                d["validation_log2FC"],
                s=10,
                color=color,
                alpha=0.7,
                linewidth=0,
                label=TIER_LABEL[tier],
                rasterized=True,
            )
        lo = np.nanmin(np.r_[data.loc[finite, "discovery_log2FC"], data.loc[finite, "validation_log2FC"]])
        hi = np.nanmax(np.r_[data.loc[finite, "discovery_log2FC"], data.loc[finite, "validation_log2FC"]])
        pad = (hi - lo) * 0.05
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color=COLORS["dark"], lw=0.7, ls="--")
        ax.axhline(0, color=COLORS["dark"], lw=0.5)
        ax.axvline(0, color=COLORS["dark"], lw=0.5)
        ax.set_xlabel("Discovery pooled log₂ fold change")
        ax.set_ylabel(f"{label} paired log₂ fold change")
        ax.set_title(f"Independent paired testing ({n_pairs} pairs)", loc="left")
        ax.legend(loc="upper left")
        panel_label(ax, letter)

    ax = axs[1, 0]
    # Robust bars from unambiguous primary; intermediate bars remain CBC current-symbol (no intermediate alias audit this round)
    def status_from_cross(tier_key: str, cohort: str) -> list[int]:
        d = cross[cross["evidence_tier"].eq(tier_key)].copy()
        if cohort == "tcga":
            est = d["tcga_estimable"].astype(str).str.upper().eq("TRUE")
            same = d["tcga_same_direction"].astype(str).str.upper().eq("TRUE")
            cand = d["tcga_candidate_family_directional"].astype(str).str.upper().eq("TRUE")
        else:
            est = d["icgc_estimable"].astype(str).str.upper().eq("TRUE")
            same = d["icgc_same_direction"].astype(str).str.upper().eq("TRUE")
            cand = d["icgc_candidate_family_directional"].astype(str).str.upper().eq("TRUE")
        validated = int((est & cand).sum())
        direction_only = int((est & same & ~cand).sum())
        direction_inconsistent = int((est & ~same).sum())
        not_estimable = int((~est).sum())
        return [validated, direction_only, direction_inconsistent, not_estimable]

    def status_from_primary(cohort: str) -> list[int]:
        if cohort == "tcga":
            est = _flag(primary["tcga_ok"])
            same = _flag(primary["tcga_same"])
            cand = _flag(primary["tcga_cand_dir"])
        else:
            est = _flag(primary["icgc_ok"])
            same = _flag(primary["icgc_same"])
            cand = _flag(primary["icgc_cand_dir"])
        return [
            int((est & cand).sum()),
            int((est & same & ~cand).sum()),
            int((est & ~same).sum()),
            int((~est).sum()),
        ]

    status_order = ["cand.-family pass", "direction only", "direction inconsistent", "not estimable"]
    xlabels = ["TCGA\nrobust", "TCGA\nintermediate", "ICGC\nrobust", "ICGC\nintermediate"]
    arrays = [
        status_from_primary("tcga"),
        status_from_cross("intermediate", "tcga"),
        status_from_primary("icgc"),
        status_from_cross("intermediate", "icgc"),
    ]
    # persist panel c counts
    pd.DataFrame(
        {
            "bar": ["TCGA_robust", "TCGA_intermediate", "ICGC_robust", "ICGC_intermediate"],
            "candidate_family_pass": [a[0] for a in arrays],
            "direction_only": [a[1] for a in arrays],
            "direction_inconsistent": [a[2] for a in arrays],
            "not_estimable": [a[3] for a in arrays],
            "bar_total": [sum(a) for a in arrays],
            "source": ["primary_unambiguous", "cbc_cross", "primary_unambiguous", "cbc_cross"],
        }
    ).to_csv(src / "Fig2c_status_counts_rebuilt.tsv", sep="\t", index=False)

    bottoms = np.zeros(4)
    status_colors = [COLORS["teal"], COLORS["blue"], COLORS["red"], COLORS["grey"]]
    for idx, (status, color) in enumerate(zip(status_order, status_colors)):
        vals = np.array([a[idx] for a in arrays])
        ax.bar(np.arange(4), vals, bottom=bottoms, color=color, width=0.68, label=status)
        bottoms += vals
    ax.set_xticks(np.arange(4), xlabels)
    ax.set_ylabel("Locked discovery-stage candidates")
    ax.set_title("Candidate-family FDR status (tier totals 173 / 119)", loc="left")
    ax.legend(loc="upper right", ncol=2, columnspacing=0.7, handletextpad=0.25)
    # denominator note
    ax.text(
        0.01,
        -0.22,
        "Primary robust (unambiguous): 157/163 jointly estimable (157/173). Sensitivities: symbol-only 127/133; inclusive aliases 161/167. Intermediate bars: current-symbol.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=4.8,
        color=COLORS["dark"],
    )
    clean_axis(ax)
    panel_label(ax, "c")

    ax = axs[1, 1]
    # Example genes that pass revised primary dual candidate-family rule
    top = primary[_flag(primary["dual_cand"])].copy()
    # attach CIs from joint/current-symbol tables when available
    jsub = joint.set_index(joint["gene_symbol"] if "gene_symbol" in joint.columns else joint.iloc[:, 0])
    for col in ["tcga_CI_low", "tcga_CI_high", "icgc_CI_low", "icgc_CI_high"]:
        if col in jsub.columns:
            top[col] = top["gene_symbol"].map(jsub[col])
        else:
            top[col] = np.nan
    top = top.assign(abs_effect=top["discovery_log2FC"].abs()).sort_values("abs_effect", ascending=False).head(10).sort_values("discovery_log2FC")
    y = np.arange(len(top))
    series = [
        ("Discovery", "discovery_log2FC", None, None, COLORS["dark"], -0.20),
        ("TCGA", "tcga_log2FC", "tcga_CI_low", "tcga_CI_high", COLORS["teal"], 0.0),
        ("ICGC", "icgc_log2FC", "icgc_CI_low", "icgc_CI_high", COLORS["orange"], 0.20),
    ]
    for name, col, low, high, color, offset in series:
        xv = top[col].to_numpy(float)
        if low:
            err = np.vstack([xv - top[low].to_numpy(float), top[high].to_numpy(float) - xv])
            ax.errorbar(xv, y + offset, xerr=err, fmt="o", ms=3.0, lw=0.8, capsize=1.5, color=color, label=name)
        else:
            ax.scatter(xv, y + offset, s=11, color=color, label=name, zorder=3)
    ax.axvline(0, color=COLORS["dark"], lw=0.7)
    ax.set_yticks(y, top["gene_symbol"])
    ax.set_xlabel("log₂ fold change (platform-specific scale)")
    ax.set_title("Examples passing dual candidate-family rule", loc="left")
    clean_axis(ax)
    panel_label(ax, "d")
    export(fig, "Fig2_independent_bulk_validation")


if __name__ == "__main__":
    fig1()
    fig2()
    print("Wrote updated Fig1 and Fig2 to", OUT)
