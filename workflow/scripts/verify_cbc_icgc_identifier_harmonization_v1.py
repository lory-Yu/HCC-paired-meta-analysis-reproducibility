#!/usr/bin/env python3
import csv
from pathlib import Path


def read_tsv(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def truth(value: str) -> bool:
    return str(value).upper() in {"TRUE", "1"}


root = Path(__file__).resolve().parents[2]
out = root / "results" / "validation" / "cbc_external_validation_v1"
icgc = read_tsv(out / "icgc_locked_candidate_identifier_harmonized_sensitivity.tsv")
cross = read_tsv(out / "tcga_icgc_identifier_harmonized_sensitivity.tsv")
counts = read_tsv(out / "identifier_harmonization_counts.tsv")

assert len(icgc) == 292
assert len(cross) == 292
assert len({(r["gene_id"], r["gene_symbol"], r["evidence_tier"]) for r in icgc}) == 292
assert len({(r["gene_id"], r["gene_symbol"], r["evidence_tier"]) for r in cross}) == 292
assert sum(r["mapping_status"] == "historical_symbol_harmonized" for r in icgc) == 34
assert {r["evidence_tier"] for r in counts} == {"discovery_robust", "intermediate"}

high = [r for r in cross if r["evidence_tier"] == "discovery_robust"]
assert len(high) == 173
assert sum(truth(r["icgc_estimable"]) for r in high) == 171
assert sum(truth(r["jointly_estimable"]) for r in high) == 167
assert sum(truth(r["dual_same_direction"]) for r in high) == 166
assert sum(truth(r["dual_candidate_family_directional"]) for r in high) == 161
assert sum(truth(r["dual_genomewide_directional"]) for r in high) == 159
assert all(not truth(r["dual_genomewide_directional"]) or truth(r["dual_same_direction"]) for r in high)

print("CBC ICGC identifier harmonization sensitivity: PASS")
