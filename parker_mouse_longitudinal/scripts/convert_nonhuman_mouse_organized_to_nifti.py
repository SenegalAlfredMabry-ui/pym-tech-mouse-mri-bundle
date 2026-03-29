#!/usr/bin/env python3
"""Convert organized nonhuman mouse DICOM series to NIfTI using dcm2niix.

Input is the organizer manifest:
  Mouse_MRI_Directory_organized/manifests/organized_manifest.tsv

Output layout:
  <out_root>/<subject_id>/<timepoint>/<reference_name>.nii.gz
  <out_root>/<subject_id>/<timepoint>/<reference_name>.json
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
import subprocess
import sys


@dataclass
class WorkRow:
    subject_id: str
    timepoint: str
    modality: str
    reference_name: str
    source_dir: Path
    out_dir: Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Convert organized nonhuman mouse DICOM series to NIfTI")
    p.add_argument(
        "--manifest-tsv",
        default="${SNHEART_ROOT:-/path/to/SNHEART}/MRI_data/Raw/nonhuman_mouse/Mouse_MRI_Directory_organized/manifests/organized_manifest.tsv",
        help="Organizer manifest TSV path",
    )
    p.add_argument(
        "--out-root",
        default="${SNHEART_ROOT:-/path/to/SNHEART}/MRI_data/Raw/nonhuman_mouse/NIFTI_longitudinal",
        help="Output root for converted NIfTI files",
    )
    p.add_argument(
        "--dcm2niix-path",
        default="${DCM2NIIX_PATH:-dcm2niix}",
        help="Path to dcm2niix binary",
    )
    p.add_argument(
        "--include-modalities",
        default="T1gre,T2_RARE,T2tse,fMRI_6000,fMRI_7000,TOF",
        help="Comma-separated modalities to convert (default keeps core scan families)",
    )
    p.add_argument("--overwrite", action="store_true", help="Re-convert even when NIfTI output exists")
    p.add_argument("--dry-run", action="store_true", help="Print planned conversions without running dcm2niix")
    return p.parse_args()


def load_rows(manifest_tsv: Path, out_root: Path, allowed_modalities: set[str]) -> list[WorkRow]:
    rows: list[WorkRow] = []
    with manifest_tsv.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for rec in reader:
            modality = rec["modality"]
            if modality not in allowed_modalities:
                continue
            subject_id = rec["subject_id"]
            timepoint = rec["timepoint"]
            reference_name = rec["reference_name"]
            source_dir = Path(rec["source_dir"])
            out_dir = out_root / subject_id / timepoint
            rows.append(
                WorkRow(
                    subject_id=subject_id,
                    timepoint=timepoint,
                    modality=modality,
                    reference_name=reference_name,
                    source_dir=source_dir,
                    out_dir=out_dir,
                )
            )
    return rows


def find_nifti_outputs(out_dir: Path, reference_name: str) -> list[Path]:
    patterns = [
        f"{reference_name}.nii.gz",
        f"{reference_name}.nii",
        f"{reference_name}_*.nii.gz",
        f"{reference_name}_*.nii",
    ]
    seen: set[Path] = set()
    out: list[Path] = []
    for pat in patterns:
        for p in sorted(out_dir.glob(pat)):
            rp = p.resolve()
            if rp not in seen:
                seen.add(rp)
                out.append(rp)
    return out


def preferred_output(paths: list[Path]) -> Path:
    def keyfunc(p: Path) -> tuple[int, int, str]:
        name = p.name.lower()
        phase_penalty = 1 if "_ph" in name else 0
        echo_penalty = 1
        if "_e1" in name:
            echo_penalty = 0
        return (phase_penalty, echo_penalty, name)

    return sorted(paths, key=keyfunc)[0]


def run_conversion(row: WorkRow, dcm2niix_path: Path, overwrite: bool, dry_run: bool) -> tuple[str, str]:
    row.out_dir.mkdir(parents=True, exist_ok=True)
    existing = find_nifti_outputs(row.out_dir, row.reference_name)
    if existing and not overwrite:
        return ("SKIP_EXISTS", str(preferred_output(existing)))

    cmd = [
        str(dcm2niix_path),
        "-f",
        row.reference_name,
        "-m",
        "n",
        "-p",
        "n",
        "-z",
        "y",
        "-ba",
        "n",
        "-o",
        str(row.out_dir),
        str(row.source_dir),
    ]
    if dry_run:
        return ("DRY_RUN", " ".join(cmd))

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        msg = (proc.stderr or proc.stdout).strip().splitlines()
        tail = msg[-1] if msg else "dcm2niix failed"
        return ("FAIL", tail)

    outputs = find_nifti_outputs(row.out_dir, row.reference_name)
    if outputs:
        return ("OK", str(preferred_output(outputs)))
    return ("FAIL", "dcm2niix returned success but output NIfTI not found")


def main() -> int:
    args = parse_args()
    manifest_tsv = Path(args.manifest_tsv).resolve()
    out_root = Path(args.out_root).resolve()
    dcm2niix_path = Path(args.dcm2niix_path).resolve()
    allowed_modalities = {m.strip() for m in args.include_modalities.split(",") if m.strip()}

    if not manifest_tsv.is_file():
        print(f"ERROR: manifest TSV missing: {manifest_tsv}", file=sys.stderr)
        return 1
    if not dcm2niix_path.is_file():
        print(f"ERROR: dcm2niix missing: {dcm2niix_path}", file=sys.stderr)
        return 1

    rows = load_rows(manifest_tsv, out_root, allowed_modalities)
    if not rows:
        print("ERROR: no conversion rows selected from manifest", file=sys.stderr)
        return 1

    logs_dir = out_root / "manifests"
    logs_dir.mkdir(parents=True, exist_ok=True)
    conv_tsv = logs_dir / "conversion_manifest.tsv"
    with conv_tsv.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "subject_id",
                "timepoint",
                "modality",
                "reference_name",
                "source_dir",
                "status",
                "detail",
            ]
        )
        ok = 0
        fail = 0
        skip = 0
        dry = 0
        for row in rows:
            status, detail = run_conversion(row, dcm2niix_path, args.overwrite, args.dry_run)
            if status == "OK":
                ok += 1
            elif status == "FAIL":
                fail += 1
            elif status == "SKIP_EXISTS":
                skip += 1
            elif status == "DRY_RUN":
                dry += 1
            w.writerow(
                [
                    row.subject_id,
                    row.timepoint,
                    row.modality,
                    row.reference_name,
                    str(row.source_dir),
                    status,
                    detail,
                ]
            )

    print(f"Rows selected: {len(rows)}")
    print(f"OK={ok} SKIP_EXISTS={skip} DRY_RUN={dry} FAIL={fail}")
    print(f"Conversion manifest: {conv_tsv}")

    if fail > 0:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
