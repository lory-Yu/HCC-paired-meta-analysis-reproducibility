#!/usr/bin/env python3
"""Write deterministic SHA-256 checksums for the complete CBC delivery folder."""

from __future__ import annotations

import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
# Standalone public release keeps the manuscript files under manuscript/; the workspace under results/.
OUT = ROOT / "manuscript" if (ROOT / "manuscript" / "source_data").is_dir() else ROOT / "results" / "manuscript_computational_biology_chemistry"
DEST = OUT / "delivery_checksums_CBC.sha256"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    files = []
    for path in OUT.rglob("*"):
        if not path.is_file() or path == DEST:
            continue
        if path.name.startswith((".~", "~$")):
            # Ignore transient office lock files; they are not delivery artifacts.
            continue
        rel = path.relative_to(OUT)
        # Renders, release archives, release plans and git checkouts are not submission files.
        if rel.parts and rel.parts[0] in {"docx_render", "public_release", "repository_release_plan"}:
            continue
        if ".git" in rel.parts:
            continue
        # Superseded collision reports in the figure directory are retained for
        # history but are not part of the delivery; current reports live in figure_qa/.
        if path.name.endswith(".pdf.collision-audit.json"):
            continue
        files.append(path)
    files.sort(key=lambda p: p.relative_to(OUT).as_posix())
    lines = [f"{sha256(path)}  {path.relative_to(OUT).as_posix()}" for path in files]
    DEST.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {DEST} ({len(lines)} files)")


if __name__ == "__main__":
    main()
