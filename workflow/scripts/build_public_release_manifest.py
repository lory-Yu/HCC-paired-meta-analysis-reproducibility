#!/usr/bin/env python3
"""Build a relative file manifest and package-level SHA-256 list."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "FILE_MANIFEST.tsv"
CHECKSUMS = ROOT / "SHA256SUMS.sha256"
EXCLUDE = {MANIFEST, CHECKSUMS}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            value.update(block)
    return value.hexdigest()


def release_files() -> list[Path]:
    return sorted(
        (
            path
            for path in ROOT.rglob("*")
            if path.is_file()
            and path not in EXCLUDE
            and ".git" not in path.relative_to(ROOT).parts
        ),
        key=lambda path: path.relative_to(ROOT).as_posix(),
    )


def main() -> None:
    files = release_files()
    with MANIFEST.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["path", "bytes", "sha256", "section"])
        for path in files:
            rel = path.relative_to(ROOT).as_posix()
            writer.writerow([rel, path.stat().st_size, digest(path), rel.split("/", 1)[0]])

    checksum_files = sorted(files + [MANIFEST], key=lambda path: path.relative_to(ROOT).as_posix())
    CHECKSUMS.write_text(
        "".join(f"{digest(path)}  {path.relative_to(ROOT).as_posix()}\n" for path in checksum_files),
        encoding="utf-8",
    )
    print(f"manifest_rows={len(files)}")
    print(f"checksum_rows={len(checksum_files)}")


if __name__ == "__main__":
    main()
