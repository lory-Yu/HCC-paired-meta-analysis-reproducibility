#!/usr/bin/env python3
"""Audit GEO sample-characteristics fields of the six discovery cohorts.

Checks, per cohort: (1) whether each field's content matches its label
(sex-like values under an age label and vice versa); (2) whether sex and age
agree between the tumour and non-tumour sample of the same patient; (3) whether
the GEO tissue field agrees with the frozen sample-manifest group; (4) for
GSE57555, whether date of birth is compatible with recorded age.
Writes a TSV of findings and prints a summary. Read-only with respect to inputs.
"""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
# Full GEO SOFT exports (about 540 MB, with expression tables) are not redistributed in
# the public release; point GEO_SOFT_DIR at re-downloaded exports to rerun the audit.
META = Path(os.environ.get("GEO_SOFT_DIR", ROOT / "00_protocol" / "source_metadata"))
MANIFEST = ROOT / "00_protocol" / "sample_manifest_v1.0.tsv"
OUT = ROOT / "results" / "audit" / "geo_metadata_field_audit_v1.tsv"
COHORTS = ["GSE121248", "GSE45114", "GSE57555", "GSE57957", "GSE76427", "GSE84402"]

SEX_VALUES = {"m", "f", "male", "female", "1", "2"}
SEX_LABEL = re.compile(r"^(sex|gender)", re.I)
AGE_LABEL = re.compile(r"^age", re.I)
TUMOUR = re.compile(r"tumou?r|cancer|carcinoma|hcc", re.I)
NONTUMOUR = re.compile(r"adjacent|non-?tumou?r|non-?cancerous|normal|pericancerous", re.I)


def parse(gse: str) -> pd.DataFrame:
    rows = []
    for block in (META / f"{gse}_full.soft").read_text().split("^SAMPLE = ")[1:]:
        gsm = block.split("\n", 1)[0].strip()
        fields: dict[str, str] = {}
        for line in block.split("\n"):
            m = re.match(r"!Sample_characteristics_ch1 = ([^:]+):\s*(.*)", line)
            if m:
                fields[m.group(1).strip()] = m.group(2).strip()
        for key, tag in (("_title", "!Sample_title = "), ("_source_name", "!Sample_source_name_ch1 = ")):
            hit = re.search(re.escape(tag) + r"(.*)", block)
            fields[key] = hit.group(1).strip() if hit else ""
        rows.append({"gsm": gsm, **fields})
    return pd.DataFrame(rows)


def main() -> None:
    missing = [g for g in COHORTS if not (META / f"{g}_full.soft").is_file()]
    if missing:
        sys.exit(f"GEO SOFT exports not found in {META} for: {', '.join(missing)}. "
                 "Download the full SOFT export of each series from NCBI GEO and set GEO_SOFT_DIR; "
                 "the audit result is archived in results/audit/geo_metadata_field_audit_v1.tsv.")
    manifest = pd.read_csv(MANIFEST, sep="\t", dtype=str)
    findings = []

    def add(gse, check, field, result, n, detail):
        findings.append(dict(cohort=gse, check=check, field=field, result=result, n=n, detail=detail))

    for gse in COHORTS:
        df = parse(gse)
        man = manifest[manifest.gse == gse].set_index("gsm")
        df = df[df.gsm.isin(man.index)].copy()
        df["patient_id"] = df.gsm.map(man.patient_id)
        df["group"] = df.gsm.map(man.group)
        df["analysed"] = df.gsm.map((man.complete_pair == "yes") & (man.status == "include"))
        fields = [c for c in df.columns if c not in {"gsm", "patient_id", "group", "analysed", "_title", "_source_name"}]

        # (1) label/content plausibility for sex- and age-labelled fields.
        sex_field = age_field = None
        mismatched = set()
        for f in fields:
            vals = df[f].dropna().str.lower()
            vals = vals[~vals.isin({"na", ""})]
            if vals.empty:
                continue
            numeric = pd.to_numeric(vals, errors="coerce")
            sexlike = vals.isin(SEX_VALUES).mean()
            if SEX_LABEL.match(f):
                ok = sexlike == 1.0
                if not ok:
                    mismatched.add(f)
                add(gse, "label_content", f, "PASS" if ok else "MISMATCH", len(vals),
                    "sex values" if ok else f"non-sex values, e.g. {', '.join(sorted(vals.unique())[:5])}")
                if ok:
                    sex_field = f
            if AGE_LABEL.match(f):
                ok = numeric.notna().all() and numeric.between(0, 110).all()
                if not ok:
                    mismatched.add(f)
                add(gse, "label_content", f, "PASS" if ok else "MISMATCH", len(vals),
                    f"numeric ages {int(numeric.min())}-{int(numeric.max())}" if ok
                    else f"non-age values, e.g. {', '.join(sorted(vals.unique())[:5])}")
                if ok:
                    age_field = f
        # Content-based sex/age only for fields whose label failed (GSE121248).
        for f in mismatched:
            vals = df[f].dropna().str.lower()
            if AGE_LABEL.match(f) and len(vals) and vals.isin(SEX_VALUES).all():
                sex_field = f
            if SEX_LABEL.match(f) and len(vals) and pd.to_numeric(vals, errors="coerce").between(0, 110).all():
                age_field = f
        if sex_field is None:
            add(gse, "label_content", "(sex)", "NOT_RECORDED", 0, "no sex/gender field in GEO characteristics")

        # (2) within-patient agreement of sex and age.
        for kind, f in (("sex", sex_field), ("age", age_field)):
            if f is None:
                continue
            per = df.dropna(subset=[f]).groupby("patient_id")[f].nunique()
            multi = (df.groupby("patient_id").size() > 1)
            checked = per[multi.reindex(per.index, fill_value=False)]
            bad = checked[checked > 1]
            add(gse, f"within_patient_{kind}", f, "PASS" if bad.empty else "CONFLICT", int(len(checked)),
                "consistent across samples of each patient" if bad.empty else f"conflicting patients: {', '.join(bad.index[:10])}")

        # (2b) sex written in the sample title, where present, vs the sex field.
        if sex_field is not None:
            title_sex = df["_title"].str.lower().str.extract(r"\b(female|male)\b")[0]
            has = title_sex.notna()
            if has.any():
                field_sex = df[sex_field].str.lower().map({"m": "male", "f": "female", "male": "male", "female": "female", "1": "male", "2": "female"})
                bad = (title_sex[has] != field_sex[has]).sum()
                add(gse, "title_vs_sex_field", f"title vs {sex_field}", "PASS" if bad == 0 else "CONFLICT", int(has.sum()),
                    "sex in sample title agrees with sex field" if bad == 0 else f"{bad} disagreements")

        # (3) tissue field vs manifest group.
        tissue_fields = [f for f in fields if re.search(r"tissue|disease state", f, re.I)] + ["_source_name", "_title"]
        best = None
        for f in tissue_fields:
            v = df[f].fillna("")
            t = v.str.contains(TUMOUR) & ~v.str.contains(NONTUMOUR)
            nt = v.str.contains(NONTUMOUR)
            if t.any() and nt.any() and (t | nt).mean() > 0.9:
                best = (f, t, nt)
                break
        if best is None:
            add(gse, "tissue_vs_manifest", "(none)", "NOT_CHECKABLE", len(df),
                "no characteristics field separates tumour from non-tumour; manifest grouping relies on title/source fields")
        else:
            f, t, nt = best
            pred = pd.Series(pd.NA, index=df.index, dtype="object")
            pred[t] = "tumor"; pred[nt] = "control"
            grp = df.group.str.lower().map(lambda g: "tumor" if "tum" in str(g) else "control")
            sub = df[df.group.isin(["tumor", "control"])]
            mism = (pred[sub.index] != grp[sub.index]).sum()
            add(gse, "tissue_vs_manifest", f, "PASS" if mism == 0 else "MISMATCH", len(sub),
                "all manifest groups agree with GEO tissue field" if mism == 0 else f"{mism} disagreements")

        # (4) GSE57555 date of birth vs age (analysed patients; others reported separately).
        if gse == "GSE57555" and {"date of birth", "age"} <= set(fields):
            d = df.dropna(subset=["date of birth", "age"]).drop_duplicates("patient_id").copy()
            d["implied"] = pd.to_datetime(d["date of birth"], format="%Y/%m/%d", errors="coerce").dt.year + pd.to_numeric(d["age"], errors="coerce")
            a = d[d.analysed == True]  # noqa: E712
            add(gse, "dob_vs_age", "date of birth + age", "PASS" if a.implied.max() - a.implied.min() <= 2 else "REVIEW",
                len(a), f"analysed patients: birth year + age = {int(a.implied.min())}-{int(a.implied.max())}")
            rest = d[d.analysed != True]  # noqa: E712
            odd = rest[(rest.implied < a.implied.min() - 5)]
            if len(odd):
                add(gse, "dob_vs_age_not_analysed", "date of birth + age", "NOTE", len(rest),
                    "not analysed (excluded CCC/unpaired): implausible implied year for "
                    + "; ".join(f"patient {pid} ({int(y)})" for pid, y in zip(odd.patient_id, odd.implied)))

    out = pd.DataFrame(findings)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
