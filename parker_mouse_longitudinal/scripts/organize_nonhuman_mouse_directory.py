#!/usr/bin/env python3
"""Organize Mouse_MRI_Directory raw DICOM folders into a subject-first layout.

Creates a symlink-based structure:
  <out_root>/by_subject/<subject>/<timepoint>/<modality>/<reference_name> -> <series_dir>

Writes:
  - organized_manifest.tsv
  - expected_structure_check.tsv
  - cohort_subjects.tsv
  - summary.txt
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Iterable


MODALITY_BY_PREFIX = {
    "3000": "T1gre",
    "4000": "T2tse",
    "6000": "fMRI_6000",
    "7000": "fMRI_7000",
    "8000": "TOF",
    "9000": "T2_RARE",
}


@dataclass
class SeriesRow:
    timepoint: str
    subject_raw: str
    subject_id: str
    exam_id: str
    series_id: str
    series_prefix: str
    modality: str
    source_dir: Path
    link_path: Path
    reference_name: str


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Organize nonhuman mouse MRI raw directory")
    p.add_argument(
        "--raw-root",
        default="${SNHEART_ROOT:-/path/to/SNHEART}/MRI_data/Raw/nonhuman_mouse/Mouse_MRI_Directory/Raw Mouse Data",
        help="Path to extracted 'Raw Mouse Data' directory",
    )
    p.add_argument(
        "--out-root",
        default="${SNHEART_ROOT:-/path/to/SNHEART}/MRI_data/Raw/nonhuman_mouse/Mouse_MRI_Directory_organized",
        help="Output root for organized symlink tree + manifests",
    )
    p.add_argument(
        "--expected-timepoints",
        default="6M,9M,16M",
        help="Comma-separated expected timepoint folder names",
    )
    p.add_argument(
        "--expected-subjects",
        default="100M1,100M3,100M4,98F1,98F2",
        help="Comma-separated expected subject IDs (without carets)",
    )
    return p.parse_args()


def normalize_subject(name: str) -> str:
    trimmed = name.rstrip("^")
    if trimmed:
        return trimmed
    return name


def series_prefix(series_id: str) -> str:
    match = re.match(r"^(\d)", series_id)
    if not match:
        return "unknown"
    return f"{match.group(1)}000"


def sanitize(text: str) -> str:
    out = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")
    return out or "na"


def iter_series(raw_root: Path) -> Iterable[SeriesRow]:
    for tp_dir in sorted(p for p in raw_root.iterdir() if p.is_dir()):
        tp = tp_dir.name
        for subj_dir in sorted(p for p in tp_dir.iterdir() if p.is_dir()):
            subj_raw = subj_dir.name
            subj_id = normalize_subject(subj_raw)
            for exam_dir in sorted(p for p in subj_dir.iterdir() if p.is_dir()):
                exam_id = exam_dir.name
                exam_token = sanitize(exam_id)
                for sdir in sorted(p for p in exam_dir.iterdir() if p.is_dir()):
                    sid = sdir.name
                    spfx = series_prefix(sid)
                    modality = MODALITY_BY_PREFIX.get(spfx, f"OTHER_{spfx}")
                    ref = "_".join(
                        [
                            sanitize(subj_id),
                            sanitize(tp),
                            sanitize(modality),
                            sanitize(sid),
                            exam_token,
                        ]
                    )
                    yield SeriesRow(
                        timepoint=tp,
                        subject_raw=subj_raw,
                        subject_id=subj_id,
                        exam_id=exam_id,
                        series_id=sid,
                        series_prefix=spfx,
                        modality=modality,
                        source_dir=sdir.resolve(),
                        link_path=Path(),
                        reference_name=ref,
                    )


def ensure_symlink(target: Path, link_path: Path) -> None:
    if link_path.exists() or link_path.is_symlink():
        # Re-point stale links to keep reruns deterministic.
        if link_path.is_symlink() and link_path.resolve() == target:
            return
        link_path.unlink()
    link_path.symlink_to(target)


def main() -> int:
    args = parse_args()
    raw_root = Path(args.raw_root).resolve()
    out_root = Path(args.out_root).resolve()
    by_subject = out_root / "by_subject"
    manifests = out_root / "manifests"
    by_subject.mkdir(parents=True, exist_ok=True)
    manifests.mkdir(parents=True, exist_ok=True)

    if not raw_root.is_dir():
        raise SystemExit(f"ERROR: raw root not found: {raw_root}")

    rows: list[SeriesRow] = []
    for row in iter_series(raw_root):
        link_dir = by_subject / row.subject_id / row.timepoint / row.modality
        link_dir.mkdir(parents=True, exist_ok=True)
        link_path = link_dir / row.reference_name
        ensure_symlink(row.source_dir, link_path)
        row.link_path = link_path
        rows.append(row)

    manifest_tsv = manifests / "organized_manifest.tsv"
    with manifest_tsv.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "subject_id",
                "subject_raw",
                "timepoint",
                "exam_id",
                "series_id",
                "series_prefix",
                "modality",
                "reference_name",
                "source_dir",
                "organized_link",
            ]
        )
        for r in rows:
            w.writerow(
                [
                    r.subject_id,
                    r.subject_raw,
                    r.timepoint,
                    r.exam_id,
                    r.series_id,
                    r.series_prefix,
                    r.modality,
                    r.reference_name,
                    str(r.source_dir),
                    str(r.link_path),
                ]
            )

    seen_tp = sorted({r.timepoint for r in rows})
    seen_subj = sorted({r.subject_id for r in rows})
    expected_tp = [x.strip() for x in args.expected_timepoints.split(",") if x.strip()]
    expected_subj = [x.strip() for x in args.expected_subjects.split(",") if x.strip()]
    missing_tp = sorted(set(expected_tp) - set(seen_tp))
    missing_subj = sorted(set(expected_subj) - set(seen_subj))
    extra_subj = sorted(set(seen_subj) - set(expected_subj))

    check_tsv = manifests / "expected_structure_check.tsv"
    with check_tsv.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["check_type", "value", "status"])
        for tp in expected_tp:
            w.writerow(["timepoint", tp, "present" if tp in seen_tp else "missing"])
        for subj in expected_subj:
            w.writerow(["subject", subj, "present" if subj in seen_subj else "missing"])
        for subj in extra_subj:
            w.writerow(["subject_extra", subj, "extra"])

    cohort_tsv = manifests / "cohort_subjects.tsv"
    with cohort_tsv.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["mouse_id", "mouse_root", "timepoint_glob"])
        for subj in seen_subj:
            w.writerow([subj, str(by_subject / subj), "*"])

    summary_txt = manifests / "summary.txt"
    with summary_txt.open("w") as f:
        f.write(f"raw_root={raw_root}\n")
        f.write(f"out_root={out_root}\n")
        f.write(f"series_rows={len(rows)}\n")
        f.write(f"timepoints={','.join(seen_tp)}\n")
        f.write(f"subjects={','.join(seen_subj)}\n")
        f.write(f"missing_timepoints={','.join(missing_tp) if missing_tp else 'none'}\n")
        f.write(f"missing_subjects={','.join(missing_subj) if missing_subj else 'none'}\n")
        f.write(f"extra_subjects={','.join(extra_subj) if extra_subj else 'none'}\n")
        f.write(f"manifest={manifest_tsv}\n")
        f.write(f"structure_check={check_tsv}\n")
        f.write(f"cohort_subjects={cohort_tsv}\n")

    print(f"Organized {len(rows)} series rows")
    print(f"Manifest: {manifest_tsv}")
    print(f"Check: {check_tsv}")
    print(f"Cohort TSV: {cohort_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
