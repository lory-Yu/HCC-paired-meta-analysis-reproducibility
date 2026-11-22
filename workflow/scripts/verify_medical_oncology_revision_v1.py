#!/usr/bin/env python3
"""Verify the locked denominator and secondary-sensitivity outputs in v0.2.0."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results" / "validation" / "medical_oncology_revision_v1"


def rows(name: str) -> list[dict[str, str]]:
    with (BASE / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    primary = rows("external_validation_primary_unambiguous_harmonized_gene_level.tsv")
    assert len(primary) == 173, f"primary candidate rows: {len(primary)}"

    joined = {row["gene_symbol"]: row for row in primary}
    assert len(joined) == 173, "primary gene symbols are not unique"

    summary = rows("external_validation_denominator_decision_summary.tsv")
    observed_pairs = {
        (row["metric"], int(row["numerator"]), int(row["denominator"]))
        for row in summary
    }
    expected_pairs = {
        ("icgc_estimable", 167, 173),
        ("jointly_estimable", 163, 173),
        ("dual_candidate_family", 157, 163),
        ("dual_genomewide", 155, 163),
    }
    assert expected_pairs.issubset(observed_pairs), observed_pairs

    ambiguity = rows("historical_symbol_ambiguous_four_decision.tsv")
    assert len(ambiguity) == 4, f"ambiguous-symbol rows: {len(ambiguity)}"

    donor = rows("singlecell_donor_unpaired_countsum_pseudobulk_results.tsv")
    assert len(donor) == 173, f"donor-level candidate rows: {len(donor)}"
    for row in donor:
        if row.get("donor_detail_json"):
            json.loads(row["donor_detail_json"])

    testable = [row for row in donor if row.get("status") == "ok"]
    assert len(testable) == 100, f"donor-level testable rows: {len(testable)}"
    adjusted_columns = [
        key
        for key in donor[0]
        if key.lower() in {"fdr", "padj", "adj_p", "p_adj", "q_value"}
        or "fdr" in key.lower()
    ]
    assert adjusted_columns, "no FDR/adjusted-p column found"
    passed = 0
    for row in testable:
        values = []
        for key in adjusted_columns:
            value = row.get(key, "")
            if value not in {"", "NA", "NaN", "nan"}:
                values.append(float(value))
        if values and min(values) < 0.05:
            passed += 1
    assert passed == 0, f"unexpected donor-level FDR passes: {passed}"

    print("medical_oncology_revision_v1=PASS")
    print("primary_candidates=173")
    print("icgc_estimable=167")
    print("jointly_estimable=163")
    print("dual_candidate_family_fdr=157")
    print("dual_genomewide_fdr=155")
    print("ambiguous_historical_symbols=4")
    print("donor_level_tests=100")
    print("donor_level_fdr_passes=0")


if __name__ == "__main__":
    main()
