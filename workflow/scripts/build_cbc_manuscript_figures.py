#!/usr/bin/env python3
"""Build the reduced, CBC-focused figure set from locked source tables.

The main figures answer only two questions: how stable is the discovery
effect, and how does it replicate under two explicitly separated FDR
families? Single-cell material is retained as a supplementary localization
audit. The script never changes statistical results.
"""

from __future__ import annotations

import json
import hashlib
import math
import os
import shutil
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "manuscript" if (ROOT / "manuscript" / "source_data").is_dir() else ROOT / "results" / "manuscript_computational_biology_chemistry"
FIG = OUT / "figures"
QA = OUT / "figure_qa"
SRC = OUT / "source_data"
SKILL_ENV = os.environ.get("NATURE_FIGURE_SCRIPTS")
SKILL = Path(SKILL_ENV) if SKILL_ENV else None
HAS_ALIGNMENT_AUDIT = SKILL is not None and SKILL.is_dir()
if HAS_ALIGNMENT_AUDIT:
    sys.path.insert(0, str(SKILL))
    from audit_panel_alignment import require_matplotlib_panel_alignment  # noqa: E402
else:
    def require_matplotlib_panel_alignment(*args, json_out=None, **kwargs):
        """Record that optional external alignment QA was unavailable."""
        if json_out is not None:
            Path(json_out).write_text(
                json.dumps({"status": "skipped", "reason": "NATURE_FIGURE_SCRIPTS not set"}, indent=2) + "\n",
                encoding="utf-8",
            )

MM = 1 / 25.4
COLORS = {
    "navy": "#315A7D",
    "blue": "#4C93C3",
    "teal": "#3C9D9B",
    "orange": "#E58B4A",
    "red": "#C95756",
    "purple": "#7A6BA8",
    "green": "#5A9A68",
    "grey": "#B7BDC4",
    "dark": "#3F4650",
    "light": "#E9EDF1",
}

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans", "sans-serif"],
    "font.size": 6.5,
    "axes.titlesize": 7.2,
    "axes.labelsize": 6.8,
    "xtick.labelsize": 5.8,
    "ytick.labelsize": 5.8,
    "legend.fontsize": 5.8,
    "axes.linewidth": 0.7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
    "figure.facecolor": "white",
})


def tsv(rel: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / rel, sep="\t", low_memory=False)


def label(ax, value: str, y: float = 1.075) -> None:
    ax.text(-0.11, y, value, transform=ax.transAxes, ha="right", va="top",
            fontsize=8, fontweight="bold", clip_on=False)


def grid(ax):
    ax.grid(axis="y", color="#DDE2E7", linewidth=0.45, alpha=0.8)
    ax.set_axisbelow(True)


def export(fig: plt.Figure, stem: str, panel_axes: list[plt.Axes], dpi: int = 600,
           require_labels: bool = True) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    require_matplotlib_panel_alignment(
        fig, json_out=QA / f"{stem}.alignment.json",
        overlay_svg=QA / f"{stem}.alignment.svg", tolerance_pt=1.5,
        gutter_tolerance_pt=1.5, require_panel_labels=require_labels,
        strict=True, exclude_axes=[a for a in fig.axes if a not in panel_axes],
    )
    fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.png", dpi=600, bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.tiff", dpi=dpi, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def copy_source(rel: str, name: str | None = None) -> None:
    SRC.mkdir(parents=True, exist_ok=True)
    src = ROOT / rel
    dst = SRC / (name or src.name)
    if src.suffix == ".gz":
        shutil.copy2(src, dst)
    else:
        shutil.copy2(src, dst)


def fig1() -> list[dict[str, str]]:
    cohorts = tsv("results/analysis/differential/cohort_differential_summary.tsv")
    meta = tsv("results/analysis/meta/meta_all_genes.tsv.gz")
    audit = tsv("results/analysis/sensitivity/candidate_downgrade_audit.tsv")
    fig, axx = plt.subplots(2, 2, figsize=(183 * MM, 145 * MM), constrained_layout=True)
    panels = list(axx.flat)

    ax = axx[0, 0]
    ax.bar(cohorts["cohort"], cohorts["n_pairs"], color=COLORS["navy"], width=0.65)
    for x, y in enumerate(cohorts["n_pairs"]):
        ax.text(x, y - 2, str(int(y)), ha="center", va="top", fontsize=5.5, color="white", fontweight="bold")
    ax.set_ylabel("Matched patient pairs")
    ax.set_title("Six prespecified discovery cohorts (167 pairs)", loc="left")
    ax.tick_params(axis="x")
    for tick in ax.get_xticklabels():
        tick.set_rotation(35); tick.set_rotation_mode("anchor"); tick.set_ha("right")
    ax.grid(False); label(ax, "a")

    ax = axx[0, 1]
    x = meta["pooled_log2FC"].astype(float).to_numpy()
    y = -np.log10(np.clip(meta["FDR"].astype(float).to_numpy(), 1e-300, 1))
    primary = meta["passes_primary_rule"].astype(str).str.upper().eq("TRUE").to_numpy()
    up = primary & (x > 0); down = primary & (x < 0)
    ax.scatter(x[~primary], y[~primary], s=3, color=COLORS["grey"], alpha=.45, linewidth=0, rasterized=True)
    ax.scatter(x[down], y[down], s=5, color=COLORS["blue"], alpha=.78, linewidth=0, rasterized=True, label="Primary down")
    ax.scatter(x[up], y[up], s=5, color=COLORS["red"], alpha=.78, linewidth=0, rasterized=True, label="Primary up")
    ax.axvline(-math.log2(1.5), color=COLORS["dark"], lw=.7, ls="--")
    ax.axvline(math.log2(1.5), color=COLORS["dark"], lw=.7, ls="--")
    ax.axhline(-math.log10(.05), color=COLORS["dark"], lw=.7, ls="--")
    ax.set_xlabel("Pooled log₂ fold change (tumour - control)")
    ax.set_ylabel("-log₁₀ meta-analysis FDR")
    ax.set_title("Random-effects meta-analysis of 17,900 genes", loc="left")
    ax.legend(loc="upper center", bbox_to_anchor=(.5, -.20), ncol=2,
              handletextpad=.3, columnspacing=.8, borderaxespad=0)
    label(ax, "b")

    ax = axx[1, 0]
    order = ["high_confidence", "intermediate_confidence", "exploratory_only"]
    colors = {"high_confidence": COLORS["teal"], "intermediate_confidence": COLORS["orange"], "exploratory_only": COLORS["grey"]}
    names = {"high_confidence": "discovery robust", "intermediate_confidence": "intermediate", "exploratory_only": "exploratory"}
    for tier in order[::-1]:
        d = audit[audit["confidence_status"] == tier]
        ax.scatter(d["pooled_log2FC"], d["I2"], s=11 if tier == "high_confidence" else 8,
                   color=colors[tier], alpha=.78, linewidth=0, label=names[tier], rasterized=True)
    ax.axhline(50, color=COLORS["dark"], lw=.7, ls="--"); ax.axhline(75, color=COLORS["dark"], lw=.7, ls=":"); ax.axvline(0, color=COLORS["dark"], lw=.6)
    ax.set_xlabel("Pooled log₂ fold change"); ax.set_ylabel("Heterogeneity, I² (%)")
    ax.set_title("Evidence tiers for 446 primary hits", loc="left")
    ax.legend(loc="upper center", bbox_to_anchor=(.5, -.20), ncol=3,
              handletextpad=.3, columnspacing=.8, borderaxespad=0)
    label(ax, "c")

    ax = axx[1, 1]
    ax.scatter(audit["pooled_log2FC"], audit["fixed_log2FC"], s=7, color=COLORS["teal"], alpha=.5, linewidth=0, rasterized=True, label="Fixed-effect model")
    ax.scatter(audit["pooled_log2FC"], audit["all_samples_log2FC"], s=7, color=COLORS["orange"], alpha=.5, linewidth=0, rasterized=True, label="Blocked 438-sample model")
    lim = np.nanmax(np.abs(np.r_[audit["pooled_log2FC"], audit["fixed_log2FC"], audit["all_samples_log2FC"]])) * 1.04
    ax.plot([-lim, lim], [-lim, lim], color=COLORS["dark"], lw=.7, ls="--"); ax.axhline(0, color=COLORS["dark"], lw=.5); ax.axvline(0, color=COLORS["dark"], lw=.5)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("Primary random-effects log₂ fold change"); ax.set_ylabel("Sensitivity-model log₂ fold change")
    ax.set_title("Effect estimates remained directionally stable", loc="left"); ax.legend(loc="upper left"); label(ax, "d")
    export(fig, "Fig1_discovery_robustness_CBC", panels)
    return [
        {"figure": "Fig1", "panel": "a", "claim": "Six cohorts contributed 167 matched pairs", "source": "cohort_differential_summary.tsv"},
        {"figure": "Fig1", "panel": "b", "claim": "446 genes met the frozen discovery rule", "source": "meta_all_genes.tsv.gz"},
        {"figure": "Fig1", "panel": "c", "claim": "The 446 hits were stratified by heterogeneity and stability", "source": "candidate_downgrade_audit.tsv"},
        {"figure": "Fig1", "panel": "d", "claim": "Sensitivity estimates retained direction", "source": "candidate_downgrade_audit.tsv"},
    ]


def fig2() -> list[dict[str, str]]:
    tcga = tsv("results/validation/cbc_external_validation_v1/tcga_locked_candidate_multiplicity.tsv")
    icgc = tsv("results/validation/cbc_external_validation_v1/icgc_locked_candidate_identifier_harmonized_sensitivity.tsv")
    cross = tsv("results/validation/cbc_external_validation_v1/tcga_icgc_identifier_harmonized_sensitivity.tsv")
    fig, axx = plt.subplots(2, 2, figsize=(183 * MM, 145 * MM), constrained_layout=True)
    panels = list(axx.flat)
    joint = cross[(cross.evidence_tier.eq("discovery_robust")) & (cross.jointly_estimable)].copy()
    for ax, ycol, name, pairs, letter in [(axx[0,0], "tcga_log2FC", "TCGA-LIHC", 50, "a"), (axx[0,1], "icgc_log2FC", "ICGC-LIRI-JP", 199, "b")]:
        xvals = joint["discovery_log2FC"].astype(float)
        yvals = joint[ycol].astype(float)
        ax.scatter(xvals, yvals, s=11, color=COLORS["teal"], alpha=.72, linewidth=0, rasterized=True)
        slope, intercept = np.polyfit(xvals, yvals, 1)
        corr = np.corrcoef(xvals, yvals)[0, 1]
        vals = np.r_[xvals, yvals]
        lo, hi = np.nanmin(vals), np.nanmax(vals); pad=(hi-lo)*.05
        ax.plot([lo-pad,hi+pad],[lo-pad,hi+pad], color=COLORS["dark"], lw=.7, ls="--")
        ax.axhline(0,color=COLORS["dark"],lw=.5); ax.axvline(0,color=COLORS["dark"],lw=.5)
        ax.set_xlabel("Discovery pooled log₂ fold change"); ax.set_ylabel(f"{name} paired log₂ fold change")
        ax.set_title(f"Jointly estimable set ({pairs} matched pairs; n=167 genes)", loc="left")
        ax.text(.04, .96, f"Pearson r={corr:.3f}\nOLS slope={slope:.2f}", transform=ax.transAxes,
                ha="left", va="top", fontsize=5.7, color=COLORS["dark"])
        label(ax, letter)

    ax = axx[1,0]
    # Explicitly show the two multiplicity families and the unavailable fraction.
    high = cross[cross.evidence_tier.eq("discovery_robust")]
    rows = [
        ("TCGA\n(candidate)", int(high.tcga_candidate_family_directional.sum()), int(high.tcga_estimable.sum())),
        ("TCGA\n(genome-wide)", int(high.tcga_genomewide_directional.sum()), int(high.tcga_estimable.sum())),
        ("ICGC\n(candidate)", int(high.icgc_candidate_family_directional.sum()), int(high.icgc_estimable.sum())),
        ("ICGC\n(genome-wide)", int(high.icgc_genomewide_directional.sum()), int(high.icgc_estimable.sum())),
        ("Both\n(candidate)", int(high.dual_candidate_family_directional.sum()), int(high.jointly_estimable.sum())),
        ("Both\n(genome-wide)", int(high.dual_genomewide_directional.sum()), int(high.jointly_estimable.sum())),
    ]
    x=np.arange(len(rows)); passed=np.array([r[1] for r in rows]); den=np.array([r[2] for r in rows])
    ax.bar(x, den, color=COLORS["light"], width=.72, label="Estimable denominator")
    ax.bar(x, passed, color=COLORS["teal"], width=.72, label="Same direction + FDR < 0.05")
    for i,(n,dn) in enumerate(zip(passed,den)):
        ax.text(i, dn+2.0, f"{n}/{dn}", ha="center", va="bottom", fontsize=5.3)
    ax.set_ylim(0, 188); ax.set_xticks(x); ax.set_xticklabels([]); ax.set_ylabel("Discovery robust genes")
    ax.set_title("Support was stable across declared FDR families", loc="left")
    short_labels=["TCGA-C", "TCGA-G", "ICGC-C", "ICGC-G", "Both-C", "Both-G"]
    for i, text_value in enumerate(short_labels):
        ax.text(i, -0.16, text_value, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=5.2, clip_on=False)
    ax.text(.5, -0.28, "C: locked-candidate FDR; G: complete-transcriptome FDR", transform=ax.transAxes, ha="center", va="top", fontsize=5.0, clip_on=False)
    ax.text(.5, -0.38, "Teal: direction + FDR < 0.05; pale cap: estimable remainder", transform=ax.transAxes, ha="center", va="top", fontsize=5.0, clip_on=False)
    ax.grid(False); ax.spines["bottom"].set_visible(False); label(ax,"c")

    ax = axx[1,1]
    z = cross[(cross.evidence_tier=="discovery_robust") & (cross.dual_genomewide_directional)].copy()
    z["abs_effect"] = z["discovery_log2FC"].abs(); z=z.sort_values("abs_effect", ascending=False).head(10).sort_values("discovery_log2FC")
    y=np.arange(len(z))
    for name,col,lo,hi,c,off in [("Discovery","discovery_log2FC",None,None,COLORS["dark"],-.20),("TCGA","tcga_log2FC",None,None,COLORS["teal"],0),("ICGC","icgc_log2FC",None,None,COLORS["orange"],.20)]:
        ax.scatter(z[col], y+off, s=11, color=c, label=name, zorder=3)
    ax.axvline(0,color=COLORS["dark"],lw=.7); ax.set_yticks(y,z["gene_symbol"]); ax.set_xlabel("log₂ fold change (platform-specific scale)")
    ax.set_title("Examples passing the dual genome-wide sensitivity rule", loc="left"); ax.grid(False); ax.legend(loc="upper left"); label(ax,"d")
    export(fig, "Fig2_dual_fdr_validation_CBC", panels)
    return [
        {"figure":"Fig2","panel":"a","claim":"TCGA effects were concordant with discovery effects among 167 jointly estimable genes","source":"tcga_icgc_identifier_harmonized_sensitivity.tsv"},
        {"figure":"Fig2","panel":"b","claim":"Identifier-harmonized ICGC effects were concordant with discovery effects among the same 167 genes","source":"tcga_icgc_identifier_harmonized_sensitivity.tsv"},
        {"figure":"Fig2","panel":"c","claim":"Candidate-family and genome-wide FDR counts were similar after identifier harmonization","source":"tcga_icgc_identifier_harmonized_sensitivity.tsv"},
        {"figure":"Fig2","panel":"d","claim":"159/167 jointly estimable discovery-robust genes passed the dual genome-wide sensitivity rule","source":"tcga_icgc_identifier_harmonized_sensitivity.tsv"},
    ]


def fig_supp_sc() -> list[dict[str,str]]:
    cells=tsv("results/validation/gse202642/cell_metadata.tsv.gz")
    dom=tsv("results/validation/gse202642/candidate_dominant_celltype.tsv")
    loc=tsv("results/validation/gse202642/candidate_celltype_localization.tsv")
    fig=plt.figure(figsize=(183*MM,120*MM),constrained_layout=False)
    gs=fig.add_gridspec(2,2,width_ratios=(1.18,1),height_ratios=(1,1),left=.07,right=.98,top=.92,bottom=.12,wspace=.34,hspace=.38)
    ax_a=fig.add_subplot(gs[:,0]); ax_b=fig.add_subplot(gs[0,1]); ax_c=fig.add_subplot(gs[1,1])
    panels=[ax_a,ax_b,ax_c]
    celltypes=sorted(cells.cell_type_broad.dropna().unique())
    palette=["#D55E00","#E69F00","#A7A500","#56B4E9","#009E73","#0072B2","#CC79A7","#7A6BA8","#5A9A68","#999999"]
    cmap=dict(zip(celltypes,palette))
    ax=ax_a
    for ct in celltypes:
        d=cells[cells.cell_type_broad.eq(ct)]
        ax.scatter(d.UMAP1,d.UMAP2,s=.22,color=cmap[ct],alpha=.52,linewidth=0,rasterized=True)
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2"); ax.set_title("GSE202642 broad cell-type map (97,255 QC cells)",loc="left"); label(ax,"a",y=1.0315)
    # The broad cell-class names are repeated as colored bars in panel b; the UMAP is
    # intentionally kept uncluttered so the supplementary panel remains readable.
    ax=ax_b
    d=dom[dom.evidence_tier.eq("high_confidence")].dominant_cell_type.value_counts().sort_values()
    yy=np.arange(len(d))*1.35
    short_class={"Cholangiocyte":"Cholangiocyte","Hepatocyte":"Hepatocyte","Kupffer_macrophage":"Kupffer/macrophage","Fibroblast_CAF":"Fibroblast/CAF","B_cell":"B cell","Endothelial":"Endothelial","Neutrophil":"Neutrophil","T_cell":"T cell","Plasma":"Plasma","NK_cell":"NK cell"}
    ax.barh(yy,d.values,color=[cmap.get(x,COLORS["grey"]) for x in d.index])
    # Draw horizontal-bar labels explicitly.  Matplotlib's rotated tick-label
    # text can be exported with a shared bounding box by the PDF collision
    # auditor even when the rendered labels are visibly separated.
    ax.set_yticks(yy)
    ax.set_yticklabels([])
    for y_value, class_name in zip(yy, d.index):
        ax.text(-0.025, y_value, short_class.get(class_name, class_name),
                transform=ax.get_yaxis_transform(), ha="right", va="center",
                fontsize=5.2, clip_on=False)
    ax.set_ylim(-.7,yy[-1]+.7); ax.set_xlabel("Genes assigned to dominant class"); ax.set_title("Dominant localization of discovery robust genes",loc="left"); ax.grid(False); label(ax,"b")
    ax=ax_c
    sel=dom[dom.evidence_tier.eq("high_confidence")].sort_values("gene_symbol").head(14).gene_symbol.tolist()
    heat=loc[loc.gene_symbol.isin(sel)].pivot(index="gene_symbol",columns="cell_type_broad",values="mean_lognorm").reindex(sel)
    z=heat.sub(heat.mean(axis=1),axis=0).div(heat.std(axis=1).replace(0,np.nan),axis=0).clip(-2.5,2.5)
    im=ax.imshow(z.to_numpy(float),aspect="auto",cmap="RdBu_r",vmin=-2.5,vmax=2.5,interpolation="nearest")
    display_class={"B_cell":"B","Cholangiocyte":"Chol","Endothelial":"Endo","Fibroblast_CAF":"CAF","Hepatocyte":"Hep","Kupffer_macrophage":"Kup","NK_cell":"NK","Neutrophil":"Neut","Plasma":"Plas","T_cell":"T"}
    ax.set_xticks(np.arange(len(z.columns)),[display_class.get(x,x.replace("_"," ")) for x in z.columns],ha="center",fontsize=5.2); ax.set_yticks(np.arange(len(z.index)),z.index); ax.tick_params(axis="x",length=0); ax.spines["bottom"].set_visible(False); ax.set_xlabel("Broad cell class (abbreviated; see caption)"); ax.set_ylabel("Gene"); ax.set_title("Expression localization is descriptive",loc="left"); label(ax,"c")
    cbar=fig.colorbar(im,ax=ax,fraction=.035,pad=.02); cbar.set_label("Within-gene z score",fontsize=5.8); cbar.ax.tick_params(labelsize=5.2)
    export(fig,"FigS1_single_cell_localization_CBC",panels)
    return [
        {"figure":"FigS1","panel":"a","claim":"Broad cell classes in 97,255 QC cells","source":"gse202642/cell_metadata.tsv.gz"},
        {"figure":"FigS1","panel":"b","claim":"Discovery robust genes localize across broad classes","source":"candidate_dominant_celltype.tsv"},
        {"figure":"FigS1","panel":"c","claim":"Localization patterns are descriptive","source":"candidate_celltype_localization.tsv"},
    ]


def graphical_abstract() -> None:
    # A data-flow schematic, not a mechanism diagram. Dimensions meet CBC's
    # 1328 x 531 pixel minimum when exported at 300 dpi.
    fig, ax = plt.subplots(figsize=(13.28,5.31),dpi=300)
    ax.set_xlim(0,13.28); ax.set_ylim(0,5.31); ax.axis("off")
    boxes=[(0.35,1.65,2.25,2.0,"Six GEO cohorts\n439 intent-to-analyse\n167 complete pairs",COLORS["navy"]),(3.0,1.65,2.25,2.0,"Paired effects\n17,900 genes\nREML meta-analysis",COLORS["teal"]),(5.65,1.65,2.25,2.0,"Robustness\nLOCO + fixed effect\n438-sample blocked model",COLORS["orange"]),(8.3,1.65,2.25,2.0,"External RNA-seq\nTCGA 50 pairs\nICGC 199 pairs",COLORS["purple"]),(10.95,1.65,2.0,2.0,"159/167 genes\npass dual\ngenome-wide rule",COLORS["green"])]
    for x,y,w,h,text,c in boxes:
        ax.add_patch(plt.Rectangle((x,y),w,h,facecolor=c,alpha=.95,edgecolor="none",joinstyle="round"))
        ax.text(x+w/2,y+h/2,text,ha="center",va="center",fontsize=10,color="white",fontweight="bold",linespacing=1.25)
    for i in range(len(boxes)-1):
        x=boxes[i][0]+boxes[i][2]; x2=boxes[i+1][0]
        ax.annotate("",xy=(x2-.08,2.65),xytext=(x+.08,2.65),arrowprops=dict(arrowstyle="->",lw=2,color=COLORS["dark"]))
    ax.text(6.64,.55,"A recurrent tumour-versus-non-tumour expression core, with multiplicity and missingness reported explicitly",ha="center",va="center",fontsize=12,color=COLORS["dark"])
    # Keep the full canvas: the CBC graphical-abstract minimum is 1328 x 531 px.
    fig.savefig(FIG/"graphical_abstract_CBC.png",dpi=300,pad_inches=.12)
    fig.savefig(FIG/"graphical_abstract_CBC.pdf",pad_inches=.12)
    fig.savefig(FIG/"graphical_abstract_CBC.svg",pad_inches=.12)
    plt.close(fig)


def main() -> None:
    FIG.mkdir(parents=True,exist_ok=True); QA.mkdir(parents=True,exist_ok=True); SRC.mkdir(parents=True,exist_ok=True)
    rows=fig1()+fig2()+fig_supp_sc(); graphical_abstract()
    source_specs = [
        "results/analysis/differential/cohort_differential_summary.tsv",
        "results/analysis/meta/meta_all_genes.tsv.gz",
        "results/analysis/sensitivity/candidate_downgrade_audit.tsv",
        "results/validation/cbc_external_validation_v1/tcga_locked_candidate_multiplicity.tsv",
        "results/validation/cbc_external_validation_v1/icgc_locked_candidate_identifier_harmonized_sensitivity.tsv",
        "results/validation/cbc_external_validation_v1/tcga_icgc_identifier_harmonized_sensitivity.tsv",
        "results/validation/cbc_external_validation_v1/external_concordance_direction.tsv",
        "results/validation/cbc_external_validation_v1/external_effect_size_concordance.tsv",
        "results/validation/cbc_external_validation_v1/tcga_filtered_discovery_robust_candidates_audit.tsv",
        "results/validation/cbc_external_validation_v1/identifier_harmonization_counts.tsv",
        "results/validation/cbc_external_validation_v1/multiplicity_counts_with_denominators.tsv",
        "results/validation/cbc_targeted_revision_audits_v1/allain2016_overlap_summary.tsv",
        "results/validation/cbc_targeted_revision_audits_v1/allain2016_overlap_gene_level.tsv",
        "results/validation/cbc_targeted_revision_audits_v1/discovery_external_sample_identifier_overlap_audit.tsv",
        "00_protocol/cbc_icgc_historical_symbol_map_v1.tsv",
        "results/validation/gse202642/candidate_dominant_celltype.tsv",
        "results/validation/gse202642/candidate_celltype_localization.tsv",
        "results/validation/gse202642/cell_metadata.tsv.gz",
    ]
    for rel in source_specs:
        copy_source(rel)
    use_by_name = {
        "cohort_differential_summary.tsv": "Fig1a",
        "meta_all_genes.tsv.gz": "Fig1b",
        "candidate_downgrade_audit.tsv": "Fig1c-d",
        "tcga_locked_candidate_multiplicity.tsv": "Fig2a; Table S1-S2",
        "icgc_locked_candidate_identifier_harmonized_sensitivity.tsv": "Fig2b; Table S1-S2",
        "tcga_icgc_identifier_harmonized_sensitivity.tsv": "Fig2a-d; Table S1",
        "external_concordance_direction.tsv": "Table S3",
        "external_effect_size_concordance.tsv": "Fig2a-b; Table S3",
        "tcga_filtered_discovery_robust_candidates_audit.tsv": "Table S2",
        "identifier_harmonization_counts.tsv": "Table S1-S2",
        "multiplicity_counts_with_denominators.tsv": "Original-symbol sensitivity",
        "allain2016_overlap_summary.tsv": "Table S4",
        "allain2016_overlap_gene_level.tsv": "Table S4 source data",
        "discovery_external_sample_identifier_overlap_audit.tsv": "Methods sample-overlap audit",
        "cbc_icgc_historical_symbol_map_v1.tsv": "Table S2 identifier mapping",
        "candidate_dominant_celltype.tsv": "FigS1b",
        "candidate_celltype_localization.tsv": "FigS1c",
        "cell_metadata.tsv.gz": "FigS1a",
    }
    manifest_rows = []
    for rel in source_specs:
        path = SRC / Path(rel).name
        manifest_rows.append({
            "file": f"source_data/{path.name}",
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "figure_or_use": use_by_name[path.name],
            "upstream_source": rel,
        })
    pd.DataFrame(manifest_rows).to_csv(OUT/"source_data_manifest_CBC.tsv",sep="\t",index=False)
    pd.DataFrame(rows).to_csv(OUT/"figure_panel_audit.tsv",sep="\t",index=False)
    (OUT/"figure_contracts.md").write_text("""# CBC figure contracts

## Fig. 1
Question: do the six paired discovery cohorts yield effects that remain stable under the frozen robustness checks? The panels show cohort structure, the complete discovery universe, tier definitions, and effect correlations. The claim is tissue-level recurrence, not causality.

## Fig. 2
Question: how much of the locked discovery set is supported in independent paired RNA-seq cohorts after transparent identifier harmonization and under explicitly declared multiplicity families? Panels a-b show effect concordance in the same jointly estimable genes, c compares candidate-family and complete-transcriptome BH, and d gives examples passing the stricter dual genome-wide sensitivity rule.

## Supplementary Fig. S1
Question: where are discovery-robust transcripts detected in one single-cell cohort? The figure is localization context only. The seven tumour and four adjacent tissue samples were not used as an independent differential-expression validation.
""",encoding="utf-8")
    qa_note = (
        "Multi-panel alignment was checked at 1.5 pt tolerance before export."
        if HAS_ALIGNMENT_AUDIT
        else "The optional external panel-alignment audit was skipped because NATURE_FIGURE_SCRIPTS was not configured; inspect exported figures manually."
    )
    (OUT/"figure_qa_notes.md").write_text(f"""# Figure QA notes

Figures were drawn with Python/matplotlib. {qa_note} Collision warnings must be reviewed at final size. The graphical abstract is a data-flow schematic and contains no mechanistic or therapeutic claim.
""",encoding="utf-8")
    (OUT/"qa"/"nature_figure_backend.json").write_text(json.dumps({"backend":"python"},indent=2)+"\n",encoding="utf-8")
    print(f"CBC figures built: {len(rows)} panel records")


if __name__ == "__main__":
    main()
