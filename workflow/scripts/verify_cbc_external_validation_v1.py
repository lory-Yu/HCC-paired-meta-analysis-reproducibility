#!/usr/bin/env python3
"""Verify the non-overwriting CBC external-validation multiplicity reconstruction."""

from __future__ import annotations

import csv
import gzip
import hashlib
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "validation" / "cbc_external_validation_v1"
HIGH_SHA = "55110e282bf9d26619de70d6e8c1e5972d050766b034a6fa06009021ae9fd3cd"
INTER_SHA = "23812a47c8df31301aca8a9061ecd7da56dd41afa5f214a3fc6e4a585494aa59"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if path.suffix == ".gz":
        handle_ctx = gzip.open(path, mode="rt", encoding="utf-8", newline="")
    else:
        handle_ctx = path.open(mode="r", encoding="utf-8", newline="")
    with handle_ctx as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def bh_adjust(values: list[float]) -> list[float]:
    n = len(values)
    order = sorted(range(n), key=values.__getitem__)
    adjusted = [math.nan] * n
    running = 1.0
    for rank_index in range(n - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        running = min(running, values[original_index] * n / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted


def as_bool(value: str) -> bool:
    assert value in {"TRUE", "FALSE"}, value
    return value == "TRUE"


def assert_bh_recalculation(rows: list[dict[str, str]], p_col: str, fdr_col: str) -> None:
    p = [float(row[p_col]) for row in rows]
    observed = [float(row[fdr_col]) for row in rows]
    expected = bh_adjust(p)
    max_delta = max(abs(a - b) for a, b in zip(observed, expected))
    assert max_delta < 1e-12, f"BH mismatch in {fdr_col}: {max_delta}"


def main() -> None:
    required = [
        "tcga_complete_filtered_transcriptome.tsv.gz",
        "icgc_complete_filtered_transcriptome.tsv.gz",
        "tcga_locked_candidate_multiplicity.tsv",
        "icgc_locked_candidate_multiplicity.tsv",
        "tcga_icgc_locked_candidate_cross_validation.tsv",
        "multiplicity_counts_with_denominators.tsv",
        "normalization_and_multiplicity_audit.tsv",
        "locked_input_manifest.tsv",
        "multiplicity_reconstruction_review_zh.md",
    ]
    for name in required:
        path = OUT / name
        assert path.is_file() and path.stat().st_size > 0, f"missing/empty: {path}"

    assert sha256(ROOT / "results/analysis/candidate_lock/candidate_set_v2_high_confidence.tsv") == HIGH_SHA
    assert sha256(ROOT / "results/analysis/candidate_lock/candidate_set_v2_intermediate.tsv") == INTER_SHA

    tcga_full = read_tsv(OUT / "tcga_complete_filtered_transcriptome.tsv.gz")
    icgc_full = read_tsv(OUT / "icgc_complete_filtered_transcriptome.tsv.gz")
    assert len(tcga_full) == 15657, len(tcga_full)
    assert len(icgc_full) == 19353, len(icgc_full)
    assert len({row["gene_id"] for row in tcga_full}) == len(tcga_full)
    assert len({row["gene_symbol"] for row in icgc_full}) == len(icgc_full)
    assert_bh_recalculation(tcga_full, "nominal_P", "genomewide_FDR")
    assert_bh_recalculation(icgc_full, "nominal_P", "genomewide_FDR")

    tcga = read_tsv(OUT / "tcga_locked_candidate_multiplicity.tsv")
    icgc = read_tsv(OUT / "icgc_locked_candidate_multiplicity.tsv")
    cross = read_tsv(OUT / "tcga_icgc_locked_candidate_cross_validation.tsv")
    assert len(tcga) == len(icgc) == len(cross) == 292
    assert sum(row["evidence_tier"] == "discovery_robust" for row in cross) == 173
    assert sum(row["evidence_tier"] == "intermediate" for row in cross) == 119

    for rows in (tcga, icgc):
        for tier, expected_n in (("discovery_robust", 173), ("intermediate", 119)):
            tier_rows = [row for row in rows if row["evidence_tier"] == tier]
            assert len(tier_rows) == expected_n
            estimable = [row for row in tier_rows if as_bool(row["estimable"])]
            assert_bh_recalculation(estimable, "nominal_P", "candidate_family_FDR")
            for row in estimable:
                same = as_bool(row["same_direction"])
                assert as_bool(row["nominal_directional"]) == (same and float(row["nominal_P"]) < 0.05)
                assert as_bool(row["candidate_family_directional"]) == (same and float(row["candidate_family_FDR"]) < 0.05)
                assert as_bool(row["genomewide_directional"]) == (same and float(row["genomewide_FDR"]) < 0.05)

    counts = read_tsv(OUT / "multiplicity_counts_with_denominators.tsv")
    index = {(row["evidence_tier"], row["scope"], row["metric"]): row for row in counts}
    for tier, total in (("discovery_robust", 173), ("intermediate", 119)):
        joint = [row for row in cross if row["evidence_tier"] == tier and as_bool(row["jointly_estimable"])]
        assert int(index[(tier, "TCGA-and-ICGC", "jointly_estimable")]["denominator"]) == total
        assert int(index[(tier, "TCGA-and-ICGC", "jointly_estimable")]["numerator"]) == len(joint)
        for metric, column in (
            ("dual_same_direction", "dual_same_direction"),
            ("dual_nominal_P_lt_0.05_and_same_direction", "dual_nominal_directional"),
            ("dual_candidate_family_FDR_lt_0.05_and_same_direction", "dual_candidate_family_directional"),
            ("dual_genomewide_FDR_lt_0.05_and_same_direction", "dual_genomewide_directional"),
        ):
            row = index[(tier, "TCGA-and-ICGC", metric)]
            assert int(row["denominator"]) == len(joint)
            assert int(row["numerator"]) == sum(as_bool(item[column]) for item in joint)

    # A known, result-independent audit anchor for this frozen input/model.
    high_tcga = [row for row in tcga if row["evidence_tier"] == "discovery_robust"]
    assert sum(as_bool(row["estimable"]) for row in high_tcga) == 168
    assert sum(as_bool(row["genomewide_directional"]) for row in high_tcga) == 161

    print("CBC external validation v1: PASS (full-universe and candidate-family BH verified)")


if __name__ == "__main__":
    main()
