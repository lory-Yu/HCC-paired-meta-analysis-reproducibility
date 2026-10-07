#!/usr/bin/env python3
"""Assemble ICGC-LIRI-JP candidate expression matrix and sample manifests."""
from __future__ import annotations

import gzip
import hashlib
import argparse
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data/external/icgc_liri_jp/raw"
HEADERS = ROOT / "data/external/icgc_liri_jp/headers"
OUT = ROOT / "results/validation/icgc_liri_jp"
OUT.mkdir(parents=True, exist_ok=True)
CAND_HIGH = ROOT / "results/analysis/candidate_lock/candidate_set_v2_high_confidence.tsv"
CAND_INT = ROOT / "results/analysis/candidate_lock/candidate_set_v2_intermediate.tsv"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_header(name: str) -> list[str]:
    with gzip.open(HEADERS / name, "rt") as f:
        return f.readline().rstrip("\n").split("\t")


def read_parts(datatype: str, header: list[str]) -> list[dict]:
    rows = []
    for path in sorted(RAW.glob(f"*/{datatype}/*.gz")):
        with gzip.open(path, "rt") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < len(header):
                    parts += [""] * (len(header) - len(parts))
                rows.append(dict(zip(header, parts)))
    return rows


def load_candidates() -> tuple[dict[str, dict], set[str]]:
    rows = {}
    symbols = set()
    for path, tier in [(CAND_HIGH, "high_confidence"), (CAND_INT, "intermediate_confidence")]:
        with open(path) as f:
            header = f.readline().rstrip("\n").split("\t")
            for line in f:
                parts = line.rstrip("\n").split("\t")
                d = dict(zip(header, [p.strip('"') for p in parts]))
                gid = d["gene_id"]
                sym = d["gene_symbol"]
                rows[gid] = {
                    "gene_id": gid,
                    "gene_symbol": sym,
                    "discovery_log2FC": float(d["pooled_log2FC"]),
                    "discovery_FDR": float(d["FDR"]),
                    "discovery_I2": float(d["I2"]),
                    "evidence_tier": tier,
                }
                symbols.add(sym)
    return rows, symbols


def is_tumor(specimen_type: str) -> bool:
    return specimen_type == "Primary tumour - solid tissue"


def is_normal(specimen_type: str) -> bool:
    return specimen_type in {
        "Normal - solid tissue",
        "Normal - tissue adjacent to primary",
    }


def main(project_root: Path | None = None) -> None:
    global ROOT, RAW, HEADERS, OUT, CAND_HIGH, CAND_INT
    if project_root is not None:
        ROOT = project_root.resolve()
    RAW = ROOT / "data/external/icgc_liri_jp/raw"
    HEADERS = ROOT / "data/external/icgc_liri_jp/headers"
    OUT = ROOT / "results/validation/icgc_liri_jp"
    OUT.mkdir(parents=True, exist_ok=True)
    CAND_HIGH = ROOT / "results/analysis/candidate_lock/candidate_set_v2_high_confidence.tsv"
    CAND_INT = ROOT / "results/analysis/candidate_lock/candidate_set_v2_intermediate.tsv"
    assert sha256(CAND_HIGH) == "55110e282bf9d26619de70d6e8c1e5972d050766b034a6fa06009021ae9fd3cd"
    assert sha256(CAND_INT) == "23812a47c8df31301aca8a9061ecd7da56dd41afa5f214a3fc6e4a585494aa59"
    candidates, symbols = load_candidates()
    print(f"candidates={len(candidates)} symbols={len(symbols)}", flush=True)

    specimen_h = read_header("specimen.tsv.gz")
    sample_h = read_header("sample.tsv.gz")
    exp_h = read_header("exp_seq.tsv.gz")
    specimens = {r["icgc_specimen_id"]: r for r in read_parts("specimen", specimen_h)}
    samples = {r["icgc_sample_id"]: r for r in read_parts("sample", sample_h)}
    print(f"specimens={len(specimens)} samples={len(samples)}", flush=True)

    # Collect expression for candidate symbols only: sample -> symbol -> raw_read_count
    expr: dict[str, dict[str, float]] = defaultdict(dict)
    sample_meta: dict[str, dict] = {}
    n_lines = 0
    for path in sorted(RAW.glob("*/exp_seq/*.gz")):
        with gzip.open(path, "rt") as f:
            for line in f:
                n_lines += 1
                parts = line.rstrip("\n").split("\t")
                if len(parts) < len(exp_h):
                    parts += [""] * (len(exp_h) - len(parts))
                d = dict(zip(exp_h, parts))
                gene = d["gene_id"]
                if gene not in symbols:
                    continue
                sid = d["icgc_sample_id"]
                try:
                    count = float(d["raw_read_count"])
                except ValueError:
                    continue
                # If duplicate gene rows for same sample, sum raw counts (prespecified).
                expr[sid][gene] = expr[sid].get(gene, 0.0) + count
                if sid not in sample_meta:
                    sample_meta[sid] = d
    print(f"scanned_exp_lines={n_lines} samples_with_candidates={len(expr)}", flush=True)

    # Annotate samples with specimen type
    annotated = []
    for sid, meta in sample_meta.items():
        samp = samples.get(sid, {})
        spid = meta.get("icgc_specimen_id") or samp.get("icgc_specimen_id", "")
        sp = specimens.get(spid, {})
        stype = sp.get("specimen_type", "")
        group = "tumor" if is_tumor(stype) else ("normal" if is_normal(stype) else "other")
        annotated.append(
            {
                "icgc_sample_id": sid,
                "icgc_specimen_id": spid,
                "icgc_donor_id": meta.get("icgc_donor_id") or samp.get("icgc_donor_id", ""),
                "submitted_sample_id": meta.get("submitted_sample_id", ""),
                "analysis_id": meta.get("analysis_id", ""),
                "specimen_type": stype,
                "group": group,
                "project_code": meta.get("project_code", "LIRI-JP"),
            }
        )

    # Deterministic selection: lexicographic first sample by submitted_sample_id then icgc_sample_id
    annotated.sort(key=lambda x: (x["icgc_donor_id"], x["group"], x["submitted_sample_id"], x["icgc_sample_id"]))
    keep = {}
    excluded = []
    for row in annotated:
        if row["group"] not in {"tumor", "normal"}:
            excluded.append({**row, "reason": "not_primary_tumour_or_liver_normal"})
            continue
        key = (row["icgc_donor_id"], row["group"])
        if key in keep:
            excluded.append({**row, "reason": "duplicate_donor_group_lexicographic_first_kept"})
            continue
        keep[key] = row

    kept_rows = list(keep.values())
    tumors = {r["icgc_donor_id"]: r for r in kept_rows if r["group"] == "tumor"}
    normals = {r["icgc_donor_id"]: r for r in kept_rows if r["group"] == "normal"}
    paired_donors = sorted(set(tumors) & set(normals))
    print(f"kept_samples={len(kept_rows)} paired_donors={len(paired_donors)}", flush=True)

    # Prefer Normal - solid tissue over adjacent if both somehow kept: already unique by group.
    # Write manifests
    def write_tsv(path: Path, rows: list[dict], cols: list[str] | None = None) -> None:
        if not rows:
            path.write_text("")
            return
        cols = cols or list(rows[0].keys())
        with open(path, "w") as f:
            f.write("\t".join(cols) + "\n")
            for r in rows:
                f.write("\t".join(str(r.get(c, "NA")) for c in cols) + "\n")

    sample_manifest = []
    for row in annotated:
        key = (row["icgc_donor_id"], row["group"])
        sample_manifest.append(
            {
                **row,
                "used_in_paired_primary": row["group"] in {"tumor", "normal"}
                and keep.get(key) is row
                and row["icgc_donor_id"] in tumors
                and row["icgc_donor_id"] in normals,
                "used_in_unpaired_secondary": row["group"] in {"tumor", "normal"} and keep.get(key) is row,
            }
        )
    write_tsv(OUT / "sample_manifest.tsv", sample_manifest)
    write_tsv(OUT / "excluded_samples.tsv", excluded)
    pairs = [
        {
            "icgc_donor_id": d,
            "tumor_sample": tumors[d]["icgc_sample_id"],
            "normal_sample": normals[d]["icgc_sample_id"],
            "tumor_specimen_type": tumors[d]["specimen_type"],
            "normal_specimen_type": normals[d]["specimen_type"],
            "tumor_submitted_sample_id": tumors[d]["submitted_sample_id"],
            "normal_submitted_sample_id": normals[d]["submitted_sample_id"],
        }
        for d in paired_donors
    ]
    write_tsv(OUT / "paired_donors.tsv", pairs)

    # Candidate matrix: genes x samples (log2(raw+1) computed later in R; store raw counts)
    all_used_samples = sorted(
        {r["icgc_sample_id"] for r in kept_rows if r["group"] in {"tumor", "normal"}}
    )
    symbols_sorted = sorted(symbols)
    with gzip.open(OUT / "candidate_expression_matrix.tsv.gz", "wt") as f:
        f.write("gene_symbol\t" + "\t".join(all_used_samples) + "\n")
        for sym in symbols_sorted:
            vals = []
            for sid in all_used_samples:
                v = expr.get(sid, {}).get(sym)
                vals.append("NA" if v is None else f"{v:.6g}")
            f.write(sym + "\t" + "\t".join(vals) + "\n")

    # Provenance summary for this assembly
    with open(OUT / "assembly_summary.tsv", "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"project\tLIRI-JP\n")
        f.write(f"release\tICGC Data Portal Release 28 (open bucket icgc25k-open)\n")
        f.write(f"exp_files\t236\n")
        f.write(f"exp_lines_scanned\t{n_lines}\n")
        f.write(f"donors_with_expression\t{len({r['icgc_donor_id'] for r in annotated})}\n")
        f.write(f"paired_donors\t{len(paired_donors)}\n")
        f.write(f"candidate_symbols\t{len(symbols)}\n")
        f.write(f"expression_unit\traw_read_count\n")
        f.write(f"gene_model\tRefSeq gene_id as symbol\n")

    print("assembly complete", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=None)
    args = parser.parse_args()
    main(args.project_root)
