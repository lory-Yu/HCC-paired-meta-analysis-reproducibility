#!/usr/bin/env python3
"""Compute descriptive external-concordance metrics for the CBC manuscript.

This is a read-only secondary audit.  It consumes the locked post-rejection
external-validation tables and writes new derived tables; it does not change
candidate locks, discovery results, or manuscript files.

The direction tests are exact binomial tests against 0.5 with
Clopper--Pearson intervals.  Because genes are correlated and were selected
before this audit, the p values are descriptive and are not treated as a
gene-level replication probability.  Effect-size associations are Pearson
and Spearman correlations and ordinary least-squares regression of external
log2FC on discovery log2FC.
"""

from __future__ import annotations

import argparse
import hashlib
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta, binomtest, linregress, pearsonr, spearmanr, t


TIERS = ["discovery_robust"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def bool_col(series: pd.Series) -> pd.Series:
    """Parse the TSV boolean convention without treating NA as False."""
    vals = series.astype("string").str.upper()
    return vals.map({"TRUE": True, "FALSE": False, "1": True, "0": False})


def cp_interval(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    if n <= 0:
        return (math.nan, math.nan)
    lo = 0.0 if k == 0 else float(beta.ppf(alpha / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - alpha / 2, k + 1, n - k))
    return lo, hi


def direction_row(
    tier: str,
    scope: str,
    cohort: str,
    metric: str,
    success: int,
    total: int,
    denominator_definition: str,
) -> dict:
    if total:
        bt = binomtest(success, total, 0.5, alternative="two-sided")
        one = binomtest(success, total, 0.5, alternative="greater")
        lo, hi = cp_interval(success, total)
        prop = success / total
    else:
        bt = one = None
        lo = hi = prop = math.nan
    return {
        "evidence_tier": tier,
        "scope": scope,
        "cohort": cohort,
        "metric": metric,
        "success": success,
        "total": total,
        "proportion": prop,
        "clopper_pearson_95_low": lo,
        "clopper_pearson_95_high": hi,
        "exact_binom_p_two_sided": math.nan if bt is None else float(bt.pvalue),
        "exact_binom_p_greater": math.nan if one is None else float(one.pvalue),
        "null_probability": 0.5,
        "denominator_definition": denominator_definition,
    }


def effect_stats(x: pd.Series, y: pd.Series) -> dict:
    frame = pd.DataFrame({"x": x, "y": y}).replace([np.inf, -np.inf], np.nan).dropna()
    n = len(frame)
    if n < 3:
        return {"n": n}
    xv = frame["x"].to_numpy(dtype=float)
    yv = frame["y"].to_numpy(dtype=float)
    pr = pearsonr(xv, yv)
    sr = spearmanr(xv, yv)
    lr = linregress(xv, yv)
    crit = float(t.ppf(0.975, n - 2))
    slope_lo = float(lr.slope - crit * lr.stderr)
    slope_hi = float(lr.slope + crit * lr.stderr)
    intercept_se = float(lr.intercept_stderr)
    intercept_lo = float(lr.intercept - crit * intercept_se)
    intercept_hi = float(lr.intercept + crit * intercept_se)
    return {
        "n": n,
        "pearson_r": float(pr.statistic),
        "pearson_p": float(pr.pvalue),
        "spearman_rho": float(sr.statistic),
        "spearman_p": float(sr.pvalue),
        "ols_slope_validation_on_discovery": float(lr.slope),
        "ols_slope_se": float(lr.stderr),
        "ols_slope_95_low": slope_lo,
        "ols_slope_95_high": slope_hi,
        "ols_slope_p": float(lr.pvalue),
        "ols_intercept": float(lr.intercept),
        "ols_intercept_se": intercept_se,
        "ols_intercept_95_low": intercept_lo,
        "ols_intercept_95_high": intercept_hi,
        "r_squared": float(lr.rvalue**2),
    }


def load_inputs(root: Path):
    out = root / "results" / "validation" / "cbc_external_validation_v1"
    tcga_path = out / "tcga_locked_candidate_multiplicity.tsv"
    icgc_path = out / "icgc_locked_candidate_multiplicity.tsv"
    cross_path = out / "tcga_icgc_locked_candidate_cross_validation.tsv"
    icgc_harmonized_path = out / "icgc_locked_candidate_identifier_harmonized_sensitivity.tsv"
    cross_harmonized_path = out / "tcga_icgc_identifier_harmonized_sensitivity.tsv"
    for path in (tcga_path, icgc_path, cross_path, icgc_harmonized_path, cross_harmonized_path):
        if not path.exists():
            raise FileNotFoundError(path)
    tcga = pd.read_csv(tcga_path, sep="\t")
    icgc = pd.read_csv(icgc_path, sep="\t")
    cross = pd.read_csv(cross_path, sep="\t")
    icgc_harmonized = pd.read_csv(icgc_harmonized_path, sep="\t")
    cross_harmonized = pd.read_csv(cross_harmonized_path, sep="\t")
    for frame in (tcga, icgc, cross, icgc_harmonized, cross_harmonized):
        for col in (
            "estimable",
            "same_direction",
            "jointly_estimable",
            "dual_same_direction",
            "dual_candidate_family_directional",
            "dual_genomewide_directional",
        ):
            if col in frame:
                frame[col] = bool_col(frame[col])
    return (
        out,
        tcga,
        icgc,
        cross,
        icgc_harmonized,
        cross_harmonized,
        (tcga_path, icgc_path, cross_path, icgc_harmonized_path, cross_harmonized_path),
    )


def build_direction_table(tcga: pd.DataFrame, icgc: pd.DataFrame, cross: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for tier in TIERS:
        tc = tcga[tcga["evidence_tier"] == tier]
        ic = icgc[icgc["evidence_tier"] == tier]
        cv = cross[cross["evidence_tier"] == tier]
        tc_est = tc[tc["estimable"] == True]
        ic_est = ic[ic["estimable"] == True]
        joint = cv[cv["jointly_estimable"] == True]
        rows.extend(
            [
                direction_row(
                    tier,
                    "all_estimable",
                    "TCGA-LIHC",
                    "same_direction",
                    int(tc_est["same_direction"].sum()),
                    len(tc_est),
                    "all discovery-robust genes estimable in TCGA",
                ),
                direction_row(
                    tier,
                    "all_estimable",
                    "ICGC-LIRI-JP",
                    "same_direction",
                    int(ic_est["same_direction"].sum()),
                    len(ic_est),
                    "all discovery-robust genes estimable in ICGC",
                ),
                direction_row(
                    tier,
                    "joint_133",
                    "TCGA-LIHC",
                    "same_direction",
                    int(joint["tcga_same_direction"].sum()),
                    len(joint),
                    "genes jointly estimable in both external cohorts",
                ),
                direction_row(
                    tier,
                    "joint_133",
                    "ICGC-LIRI-JP",
                    "same_direction",
                    int(joint["icgc_same_direction"].sum()),
                    len(joint),
                    "genes jointly estimable in both external cohorts",
                ),
                direction_row(
                    tier,
                    "joint_133",
                    "TCGA-and-ICGC",
                    "dual_same_direction",
                    int(joint["dual_same_direction"].sum()),
                    len(joint),
                    "genes jointly estimable and same direction in both cohorts",
                ),
                direction_row(
                    tier,
                    "joint_133",
                    "TCGA-and-ICGC",
                    "dual_candidate_family_FDR_and_same_direction",
                    int(joint["dual_candidate_family_directional"].sum()),
                    len(joint),
                    "genes jointly estimable; candidate-family FDR<0.05 and direction in both cohorts",
                ),
                direction_row(
                    tier,
                    "joint_133",
                    "TCGA-and-ICGC",
                    "dual_genomewide_FDR_and_same_direction",
                    int(joint["dual_genomewide_directional"].sum()),
                    len(joint),
                    "genes jointly estimable; genome-wide FDR<0.05 and direction in both cohorts",
                ),
            ]
        )
    return pd.DataFrame(rows)


def build_effect_table(tcga: pd.DataFrame, icgc: pd.DataFrame, cross: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for tier in TIERS:
        tc = tcga[tcga["evidence_tier"] == tier]
        ic = icgc[icgc["evidence_tier"] == tier]
        cv = cross[cross["evidence_tier"] == tier]
        # Use the cross-validation table for the gene intersection so that
        # each pair uses precisely the same locked gene universe.
        for scope, genes in (
            ("all_estimable", None),
            ("joint_133", cv[cv["jointly_estimable"] == True]),
        ):
            if genes is None:
                tc_use = tc[tc["estimable"] == True]
                ic_use = ic[ic["estimable"] == True]
            else:
                tc_use = genes[["gene_symbol", "discovery_log2FC", "tcga_log2FC"]].rename(
                    columns={"tcga_log2FC": "validation_log2FC"}
                )
                ic_use = genes[["gene_symbol", "discovery_log2FC", "icgc_log2FC"]].rename(
                    columns={"icgc_log2FC": "validation_log2FC"}
                )
            for cohort, frame in (("TCGA-LIHC", tc_use), ("ICGC-LIRI-JP", ic_use)):
                s = effect_stats(frame["discovery_log2FC"], frame["validation_log2FC"])
                rows.append(
                    {
                        "evidence_tier": tier,
                        "scope": scope,
                        "cohort": cohort,
                        **s,
                    }
                )
    return pd.DataFrame(rows)


def append_identifier_harmonized_metrics(
    direction: pd.DataFrame,
    effects: pd.DataFrame,
    icgc: pd.DataFrame,
    cross: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    tier = "discovery_robust"
    ic = icgc[icgc["evidence_tier"] == tier]
    cv = cross[cross["evidence_tier"] == tier]
    ic_est = ic[ic["estimable"] == True]
    joint = cv[cv["jointly_estimable"] == True]
    extra_direction = pd.DataFrame(
        [
            direction_row(
                tier,
                "identifier_harmonized_all_estimable",
                "ICGC-LIRI-JP",
                "same_direction",
                int(ic_est["same_direction"].sum()),
                len(ic_est),
                "all discovery-robust genes estimable after post hoc identifier harmonization",
            ),
            direction_row(
                tier,
                "identifier_harmonized_joint_167",
                "TCGA-LIHC",
                "same_direction",
                int(joint["tcga_same_direction"].sum()),
                len(joint),
                "genes jointly estimable after post hoc ICGC identifier harmonization",
            ),
            direction_row(
                tier,
                "identifier_harmonized_joint_167",
                "ICGC-LIRI-JP",
                "same_direction",
                int(joint["icgc_same_direction"].sum()),
                len(joint),
                "genes jointly estimable after post hoc ICGC identifier harmonization",
            ),
            direction_row(
                tier,
                "identifier_harmonized_joint_167",
                "TCGA-and-ICGC",
                "dual_same_direction",
                int(joint["dual_same_direction"].sum()),
                len(joint),
                "genes jointly estimable and same direction in both cohorts after identifier harmonization",
            ),
            direction_row(
                tier,
                "identifier_harmonized_joint_167",
                "TCGA-and-ICGC",
                "dual_candidate_family_FDR_and_same_direction",
                int(joint["dual_candidate_family_directional"].sum()),
                len(joint),
                "jointly estimable genes; candidate-family FDR<0.05 and direction in both cohorts after harmonization",
            ),
            direction_row(
                tier,
                "identifier_harmonized_joint_167",
                "TCGA-and-ICGC",
                "dual_genomewide_FDR_and_same_direction",
                int(joint["dual_genomewide_directional"].sum()),
                len(joint),
                "jointly estimable genes; genome-wide FDR<0.05 and direction in both cohorts after harmonization",
            ),
        ]
    )
    extra_effects = []
    extra_effects.append(
        {
            "evidence_tier": tier,
            "scope": "identifier_harmonized_all_estimable",
            "cohort": "ICGC-LIRI-JP",
            **effect_stats(ic_est["discovery_log2FC"], ic_est["validation_log2FC"]),
        }
    )
    for cohort, y_col in (("TCGA-LIHC", "tcga_log2FC"), ("ICGC-LIRI-JP", "icgc_log2FC")):
        extra_effects.append(
            {
                "evidence_tier": tier,
                "scope": "identifier_harmonized_joint_167",
                "cohort": cohort,
                **effect_stats(joint["discovery_log2FC"], joint[y_col]),
            }
        )
    return (
        pd.concat([direction, extra_direction], ignore_index=True),
        pd.concat([effects, pd.DataFrame(extra_effects)], ignore_index=True),
    )


def write_review(
    path: Path,
    direction: pd.DataFrame,
    effects: pd.DataFrame,
    input_paths: tuple[Path, ...],
    script_path: Path,
) -> None:
    def f(x, digits=4):
        return "NA" if pd.isna(x) else f"{float(x):.{digits}f}"

    def fp(x):
        return "NA" if pd.isna(x) else f"{float(x):.3e}"

    lines = [
        "# External concordance audit (CBC v1)",
        "",
        "This is a derived, descriptive audit of the locked post-rejection external-validation tables. It does not alter the discovery analysis, candidate locks, or manuscript files.",
        "",
        "## Statistical definitions",
        "",
        "- Direction concordance is the number of external log2FC estimates having the same sign as the frozen discovery log2FC. Exact binomial tests use a 0.5 null; intervals are two-sided 95% Clopper--Pearson intervals.",
        "- Effect-size concordance is assessed with Pearson correlation, Spearman correlation, and ordinary least-squares regression of external log2FC on discovery log2FC.",
        "- The `all_estimable` scope uses every discovery-robust gene estimable in that external cohort. The `joint_133` scope uses only genes estimable in both TCGA-LIHC and ICGC-LIRI-JP.",
        "- The gene-level binomial p values are descriptive: genes are correlated and the locked candidate set was selected before this audit. They are not evidence that genes are independent biological replicates.",
        "",
        "## Direction results",
        "",
        "| scope | cohort | metric | success/total | proportion | CP 95% CI | exact two-sided p |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for _, r in direction.iterrows():
        lines.append(
            f"| {r.scope} | {r.cohort} | {r.metric} | {int(r.success)}/{int(r.total)} | {f(r.proportion,4)} | {f(r.clopper_pearson_95_low,4)}–{f(r.clopper_pearson_95_high,4)} | {fp(r.exact_binom_p_two_sided)} |"
        )
    lines += [
        "",
        "## Effect-size results",
        "",
        "| scope | cohort | n | Pearson r | Spearman rho | OLS slope (95% CI) | intercept (95% CI) | R² |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for _, r in effects.iterrows():
        lines.append(
            f"| {r.scope} | {r.cohort} | {int(r.n)} | {f(r.pearson_r,4)} | {f(r.spearman_rho,4)} | {f(r.ols_slope_validation_on_discovery,4)} ({f(r.ols_slope_95_low,4)}–{f(r.ols_slope_95_high,4)}) | {f(r.ols_intercept,4)} ({f(r.ols_intercept_95_low,4)}–{f(r.ols_intercept_95_high,4)}) | {f(r.r_squared,4)} |"
        )
    lines += [
        "",
        "The slopes are descriptive scale comparisons, not calibration coefficients or causal effects. Values above one indicate that the external effect estimates are larger in magnitude on average than the discovery estimates under the respective platform and model specifications.",
        "",
        "## Provenance",
        "",
        f"- Script: `{script_path}`",
    ]
    for p in input_paths:
        lines.append(f"- Input SHA-256: `{p.name}` = `{sha256(p)}`")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=".")
    args = parser.parse_args()
    root = Path(args.project_root).resolve()
    out, tcga, icgc, cross, icgc_harmonized, cross_harmonized, input_paths = load_inputs(root)
    direction = build_direction_table(tcga, icgc, cross)
    effects = build_effect_table(tcga, icgc, cross)
    direction, effects = append_identifier_harmonized_metrics(
        direction, effects, icgc_harmonized, cross_harmonized
    )
    direction_path = out / "external_concordance_direction.tsv"
    effect_path = out / "external_effect_size_concordance.tsv"
    review_path = out / "external_concordance_review_zh.md"
    direction.to_csv(direction_path, sep="\t", index=False, float_format="%.12g")
    effects.to_csv(effect_path, sep="\t", index=False, float_format="%.12g")
    write_review(review_path, direction, effects, input_paths, Path(__file__).resolve())
    print(f"Wrote {direction_path}")
    print(f"Wrote {effect_path}")
    print(f"Wrote {review_path}")


if __name__ == "__main__":
    main()
