#!/usr/bin/env python3
"""Validate raw-file identity against the frozen GEO sample manifest."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
from collections import Counter, defaultdict
from pathlib import Path


RAW_EXTENSIONS = {".cel", ".gpr", ".txt", ".idat"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True)
    args = parser.parse_args()
    root = Path(args.project_root).resolve()
    frozen_path = root / "00_protocol" / "sample_manifest_v1.0.tsv"
    output_dir = root / "results" / "qc"
    provenance_dir = root / "results" / "provenance"
    output_dir.mkdir(parents=True, exist_ok=True)
    provenance_dir.mkdir(parents=True, exist_ok=True)

    with frozen_path.open(newline="", encoding="utf-8") as handle:
        frozen = list(csv.DictReader(handle, delimiter="\t"))
    included = {row["gsm"]: row for row in frozen if row["status"] == "include"}
    expected_by_cohort: dict[str, set[str]] = defaultdict(set)
    for gsm, row in included.items():
        expected_by_cohort[row["gse"]].add(gsm)

    observed: dict[str, list[Path]] = defaultdict(list)
    observed_by_cohort: dict[str, set[str]] = defaultdict(set)
    file_rows: list[dict[str, str]] = []
    for cohort in expected_by_cohort:
        if cohort in {"GSE57957", "GSE76427"}:
            matrix_path = root / "data" / "downloads" / cohort / f"{cohort}_non-normalized.txt.gz"
            if matrix_path.exists():
                with gzip.open(matrix_path, "rt", encoding="utf-8", errors="replace") as handle:
                    header = handle.readline().rstrip("\r\n").split("\t")
                matrix_labels = {name.removesuffix(".AVG_Signal") for name in header if name.endswith(".AVG_Signal")}
                cohort_rows = [row for row in frozen if row["gse"] == cohort and row["status"] == "include"]
                for row in cohort_rows:
                    if cohort == "GSE57957":
                        match = re.search(r"HCC_(\d+)([TN])", row["source_name"], flags=re.IGNORECASE)
                        label = f"{match.group(2).upper()}{match.group(1)}" if match else ""
                    else:
                        number = row["patient_id"].removeprefix("HCC")
                        # GSE76427 uses PT for tumour and ANTT for adjacent
                        # non-tumour tissue in the deposited non-normalized matrix.
                        label = ("PT" if row["group"] == "tumor" else "ANTT") + number
                    if label in matrix_labels:
                        observed_by_cohort[cohort].add(row["gsm"])
                file_rows.append({
                    "cohort": cohort,
                    "gsm": "MULTI_SAMPLE_MATRIX",
                    "path": str(matrix_path),
                    "size_bytes": str(matrix_path.stat().st_size),
                    "sha256": sha256(matrix_path),
                    "expected": "yes",
                })
            continue
        raw_dir = root / "data" / "work" / cohort / "raw"
        for path in sorted(raw_dir.rglob("*")):
            if not path.is_file() or path.name.startswith(".") or path.name.endswith(".gz"):
                continue
            gsm_match = re.search(r"(GSM\d+)", path.name, flags=re.IGNORECASE)
            if not gsm_match:
                continue
            gsm = gsm_match.group(1).upper()
            observed[gsm].append(path)
            observed_by_cohort[cohort].add(gsm)
            file_rows.append({
                "cohort": cohort,
                "gsm": gsm,
                "path": str(path),
                "size_bytes": str(path.stat().st_size),
                "sha256": sha256(path),
                "expected": "yes" if gsm in expected_by_cohort[cohort] else "no",
            })

    validation_rows: list[dict[str, str]] = []
    summary_rows: list[dict[str, str]] = []
    failures = 0
    for cohort, expected in expected_by_cohort.items():
        observed_ids = observed_by_cohort[cohort]
        missing = sorted(expected - observed_ids)
        unexpected = sorted(observed_ids - expected)
        duplicates = [] if cohort in {"GSE57957", "GSE76427"} else sorted(
            gsm for gsm in expected if len([p for p in observed[gsm] if cohort in str(p)]) > 1
        )
        for gsm in sorted(expected | observed_ids):
            if gsm in missing:
                status, reason = "fail", "missing eligible raw sample file"
                failures += 1
            elif gsm in unexpected:
                status, reason = "review", "raw file is not in eligible frozen set"
            elif gsm in duplicates:
                status, reason = "review", "multiple raw files map to one eligible GSM"
            else:
                status, reason = "pass", "exact eligible GSM match"
            validation_rows.append({"cohort": cohort, "gsm": gsm, "status": status, "reason": reason})
        summary_rows.append({
            "cohort": cohort,
            "expected_eligible_gsm": str(len(expected)),
            "observed_eligible_gsm": str(len(expected & observed_ids)),
            "missing_gsm": str(len(missing)),
            "unexpected_gsm": str(len(unexpected)),
            "duplicate_gsm": str(len(duplicates)),
            "identity_gate": "fail" if missing else "pass",
        })

    def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    write_tsv(provenance_dir / "raw_file_manifest.tsv", file_rows,
              ["cohort", "gsm", "path", "size_bytes", "sha256", "expected"])
    write_tsv(output_dir / "input_validation.tsv", validation_rows,
              ["cohort", "gsm", "status", "reason"])
    write_tsv(output_dir / "input_validation_summary.tsv", summary_rows,
              ["cohort", "expected_eligible_gsm", "observed_eligible_gsm", "missing_gsm",
               "unexpected_gsm", "duplicate_gsm", "identity_gate"])

    exclusion_path = output_dir / "qc_exclusions.tsv"
    write_tsv(exclusion_path, [],
              ["cohort", "gsm", "decision", "rule", "evidence_1", "evidence_2", "decision_time", "notes"])
    if failures:
        raise SystemExit(f"identity validation failed for {failures} eligible samples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
