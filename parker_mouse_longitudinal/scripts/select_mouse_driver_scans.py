#!/usr/bin/env python3
"""Select one driver scan per mouse timepoint directory.

Selection rule (descending priority):
1) Exclude obvious phase images (_ph).
2) Prefer named high-res 3D sequences (ISO100/Gre3D).
3) Prefer near-isotropic voxels (lower anisotropy ratio).
4) Prefer smaller voxel volume.
5) Prefer larger matrix (more voxels).
6) Stable filename tie-break.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Iterable

import nibabel as nib


@dataclass
class Candidate:
    mouse_id: str
    timepoint: str
    path: Path
    vx: float
    vy: float
    vz: float
    nx: int
    ny: int
    nz: int
    anisotropy: float
    voxel_volume: float
    voxels: int
    sequence_rank: int
    series_id: str
    series_prefix: str
    modality_hint: str


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Pick high-resolution driver scans per mouse timepoint")
    p.add_argument("--mouse-root", required=True, help="Mouse directory containing timepoint folders")
    p.add_argument("--output-tsv", required=True, help="Output TSV path")
    p.add_argument(
        "--timepoint-glob",
        default="*",
        help="Glob for timepoint directories under mouse root (default: *)",
    )
    p.add_argument(
        "--series-prefix-priority",
        default="3000,9000,4000",
        help="Preferred series-prefix order for driver scans (default: 3000,9000,4000)",
    )
    p.add_argument(
        "--exclude-prefixes",
        default="6000,7000,8000",
        help="Series prefixes to avoid for driver selection when alternatives exist (default: 6000,7000,8000)",
    )
    return p.parse_args()


def sequence_rank(name: str) -> int:
    lower = name.lower()
    if "_ph" in lower:
        return 999
    if "iso100" in lower or "gre3dinvivo" in lower:
        return 0
    if "rarevfl" in lower:
        return 1
    if "turborare" in lower or "tse" in lower:
        return 2
    return 3


def parse_series_id(name: str) -> str:
    match = re.search(r"_(\d{4,6})_", name)
    if not match:
        return ""
    return match.group(1)


def series_prefix(series_id: str) -> str:
    if not series_id:
        return ""
    return f"{series_id[0]}000"


def modality_hint(prefix: str, name: str) -> str:
    mapped = {
        "3000": "T1gre",
        "4000": "T2tse",
        "6000": "fMRI",
        "7000": "fMRI",
        "8000": "TOF",
        "9000": "T2_RARE",
    }
    if prefix in mapped:
        return mapped[prefix]

    lower = name.lower()
    if "gre3d" in lower:
        return "T1gre"
    if "rarevfl" in lower or "rare" in lower:
        return "T2_RARE"
    if "turborare" in lower or "tse" in lower:
        return "T2tse"
    return "unknown"


def normalize_token(text: str) -> str:
    out = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")
    return out or "na"


def nifti_candidates(tp_dir: Path, mouse_id: str) -> Iterable[Candidate]:
    for nii in sorted(tp_dir.glob("*.nii.gz")):
        srank = sequence_rank(nii.name)
        if srank >= 999:
            continue
        sid = parse_series_id(nii.name)
        spfx = series_prefix(sid)
        mhint = modality_hint(spfx, nii.name)
        img = nib.load(str(nii))
        zooms = img.header.get_zooms()[:3]
        shape = img.shape[:3]
        vx, vy, vz = (abs(float(z)) for z in zooms)
        nx, ny, nz = (int(s) for s in shape)
        vmin = min(vx, vy, vz)
        vmax = max(vx, vy, vz)
        anis = (vmax / vmin) if vmin > 0 else 1e9
        vvol = vx * vy * vz
        nvox = nx * ny * nz
        yield Candidate(
            mouse_id=mouse_id,
            timepoint=tp_dir.name,
            path=nii.resolve(),
            vx=vx,
            vy=vy,
            vz=vz,
            nx=nx,
            ny=ny,
            nz=nz,
            anisotropy=anis,
            voxel_volume=vvol,
            voxels=nvox,
            sequence_rank=srank,
            series_id=sid,
            series_prefix=spfx,
            modality_hint=mhint,
        )


def main() -> int:
    args = parse_args()
    mouse_root = Path(args.mouse_root).resolve()
    out_tsv = Path(args.output_tsv).resolve()
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    mouse_id = mouse_root.name

    preferred = [p.strip() for p in args.series_prefix_priority.split(",") if p.strip()]
    preferred_order = {prefix: idx for idx, prefix in enumerate(preferred)}
    excluded = {p.strip() for p in args.exclude_prefixes.split(",") if p.strip()}

    tp_dirs = sorted(p for p in mouse_root.glob(args.timepoint_glob) if p.is_dir())
    rows = []

    for tp in tp_dirs:
        cands = list(nifti_candidates(tp, mouse_id))
        if not cands:
            continue
        non_excluded = [c for c in cands if c.series_prefix not in excluded]
        if not non_excluded:
            non_excluded = cands

        preferred_pool = [c for c in non_excluded if c.series_prefix in preferred_order]
        pool = preferred_pool if preferred_pool else non_excluded

        pool.sort(
            key=lambda c: (
                preferred_order.get(c.series_prefix, len(preferred_order) + 1),
                c.sequence_rank,
                c.anisotropy,
                c.voxel_volume,
                -c.voxels,
                c.path.name,
            )
        )
        best = pool[0]
        rows.append(best)

    with out_tsv.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "timepoint",
                "driver_nifti",
                "sequence_rank",
                "vx_mm",
                "vy_mm",
                "vz_mm",
                "anisotropy",
                "voxel_volume_mm3",
                "nx",
                "ny",
                "nz",
                "n_voxels",
                "mouse_id",
                "series_id",
                "series_prefix",
                "modality_hint",
                "driver_reference",
            ]
        )
        for r in rows:
            sid_token = r.series_id if r.series_id else "seriesNA"
            ref = "_".join(
                [
                    normalize_token(r.mouse_id),
                    normalize_token(r.timepoint),
                    normalize_token(r.modality_hint),
                    normalize_token(sid_token),
                ]
            )
            w.writerow(
                [
                    r.timepoint,
                    str(r.path),
                    r.sequence_rank,
                    f"{r.vx:.6f}",
                    f"{r.vy:.6f}",
                    f"{r.vz:.6f}",
                    f"{r.anisotropy:.6f}",
                    f"{r.voxel_volume:.6f}",
                    r.nx,
                    r.ny,
                    r.nz,
                    r.voxels,
                    r.mouse_id,
                    r.series_id,
                    r.series_prefix,
                    r.modality_hint,
                    ref,
                ]
            )

    print(f"Wrote {len(rows)} selected scans to {out_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
