#!/usr/bin/env python3
"""Fail when v0.2.0 additions contain private paths or submission-process material."""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SELF = Path(__file__).resolve()

FORBIDDEN_FILENAME = re.compile(
    r"(?:cover.?letter|decision.?letter|reviewer|editorial|rejection|submission.?form|private.?access)",
    re.IGNORECASE,
)
FORBIDDEN_BYTES = [
    ("author workstation path", re.compile(rb"/(?:home|Users)/[^/\s]+/")),
    ("temporary directory path", re.compile(rb"/tmp/")),
    ("submission-process label", re.compile(rb"post[_ -]?rejection|editorial[_ -]?response", re.IGNORECASE)),
    ("private-token label", re.compile(rb"(?:api[_ -]?key|access[_ -]?token|client[_ -]?secret)\s*[:=]", re.IGNORECASE)),
    ("private telephone number", re.compile(rb"85326195")),
    ("workspace-assistant label", re.compile(rb"ChatGPT|Codex|Cursor", re.IGNORECASE)),
]


def files() -> list[Path]:
    return sorted(
        path
        for path in ROOT.rglob("*")
        if path.is_file()
        and path.resolve() != SELF
        and ".git" not in path.relative_to(ROOT).parts
    )


def main() -> None:
    failures: list[str] = []
    checked = files()
    for path in checked:
        rel = path.relative_to(ROOT).as_posix()
        if FORBIDDEN_FILENAME.search(path.name):
            failures.append(f"forbidden filename: {rel}")
        data = path.read_bytes()
        for label, pattern in FORBIDDEN_BYTES:
            if pattern.search(data):
                failures.append(f"{label}: {rel}")

    if failures:
        raise SystemExit("privacy audit failed:\n" + "\n".join(failures))
    print("release_privacy_audit_v020=PASS")
    print(f"files_checked={len(checked)}")


if __name__ == "__main__":
    main()
