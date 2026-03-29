#!/usr/bin/env python3
"""QC and run-metrics for integrated multimodal mouse longitudinal outputs.

Inputs:
- alignment_manifest.tsv from mouse_longitudinal_fmri_t2tse.sh

Outputs:
- qc_gates.tsv: row-level pass/fail with coverage + centroid checks
- timepoint_qc_gates.tsv: aggregated pass/fail by timepoint
- fmri_run_metrics.tsv: FD/DVARS/tSNR and censor percentages for 4D fMRI runs
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import nibabel as nib
import numpy as np


@dataclass
class GateConfig:
    t2_min_coverage: float = 0.40
    t2_max_centroid_mm: float = 5.0
    fmri_min_coverage: float = 0.02
    fmri_max_centroid_mm: float = 15.0
    fd_threshold_mm: float = 0.20
    dvars_z_threshold: float = 2.0
    tsnr_min_median: float = 10.0


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compute multimodal QC gates and fMRI metrics")
    p.add_argument("--alignment-manifest", required=True, help="alignment_manifest.tsv path")
    p.add_argument("--out-dir", required=True, help="Output directory for QC tables")
    p.add_argument("--fd-threshold-mm", type=float, default=0.20, help="FD proxy threshold (mm)")
    p.add_argument(
        "--dvars-z-threshold",
        type=float,
        default=2.0,
        help="DVARS robust z threshold for censoring",
    )
    p.add_argument("--tsnr-min-median", type=float, default=10.0, help="tSNR median minimum for fmri pass")
    return p.parse_args()


def load_data(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    img = nib.load(str(path))
    data = np.asarray(img.get_fdata(dtype=np.float32))
    return data, img.affine


def mean3d(data: np.ndarray) -> np.ndarray:
    if data.ndim == 3:
        return data
    if data.ndim == 4:
        return np.mean(data, axis=3)
    raise ValueError(f"Unsupported ndim: {data.ndim}")


def mask_nonzero(data3d: np.ndarray) -> np.ndarray:
    return np.isfinite(data3d) & (np.abs(data3d) > 0)


def centroid_mm(mask: np.ndarray, affine: np.ndarray) -> np.ndarray | None:
    idx = np.argwhere(mask)
    if idx.size == 0:
        return None
    center_vox = idx.mean(axis=0)
    center_mm = nib.affines.apply_affine(affine, center_vox)
    return np.asarray(center_mm, dtype=np.float64)


def corrcoef_safe(a: np.ndarray, b: np.ndarray) -> float:
    if a.size < 10 or b.size < 10:
        return float("nan")
    av = a - np.mean(a)
    bv = b - np.mean(b)
    den = np.sqrt(np.sum(av * av) * np.sum(bv * bv))
    if den == 0:
        return float("nan")
    return float(np.sum(av * bv) / den)


def robust_z(values: np.ndarray) -> np.ndarray:
    if values.size == 0:
        return values
    med = np.median(values)
    mad = np.median(np.abs(values - med))
    if mad == 0:
        std = np.std(values)
        if std == 0:
            return np.zeros_like(values)
        return (values - med) / std
    return 0.6745 * (values - med) / mad


def compute_fd_dvars_tsnr(img4d: np.ndarray, mask3d: np.ndarray) -> Dict[str, float]:
    n_vol = img4d.shape[3]
    vox = img4d[mask3d, :]
    if vox.size == 0 or n_vol < 2:
        return {
            "n_volumes": float(n_vol),
            "mean_fd_proxy_mm": float("nan"),
            "p95_fd_proxy_mm": float("nan"),
            "mean_dvars": float("nan"),
            "p95_dvars": float("nan"),
            "mean_tsnr": float("nan"),
            "median_tsnr": float("nan"),
            "pct_censored_fd": float("nan"),
            "pct_censored_dvars": float("nan"),
            "pct_censored_union": float("nan"),
        }

    # FD proxy from intensity-weighted center-of-mass shifts (in voxel units).
    fds = []
    for t in range(n_vol):
        vol = img4d[..., t]
        w = np.maximum(vol, 0)
        s = np.sum(w)
        if s <= 0:
            fds.append(np.array([np.nan, np.nan, np.nan]))
            continue
        coords = np.indices(vol.shape, dtype=np.float32)
        cm = np.array([np.sum(coords[i] * w) / s for i in range(3)], dtype=np.float64)
        fds.append(cm)
    fds = np.asarray(fds)
    diffs = np.linalg.norm(np.diff(fds, axis=0), axis=1)
    diffs = diffs[np.isfinite(diffs)]

    # DVARS from temporal derivative over masked voxels.
    dvox = np.diff(vox, axis=1)
    dvars = np.sqrt(np.mean(dvox * dvox, axis=0))

    # tSNR from mean/std across time.
    m = np.mean(vox, axis=1)
    s = np.std(vox, axis=1)
    tsnr = np.divide(m, s, out=np.full_like(m, np.nan), where=s > 0)

    return {
        "n_volumes": float(n_vol),
        "mean_fd_proxy_mm": float(np.mean(diffs)) if diffs.size else float("nan"),
        "p95_fd_proxy_mm": float(np.percentile(diffs, 95)) if diffs.size else float("nan"),
        "mean_dvars": float(np.mean(dvars)) if dvars.size else float("nan"),
        "p95_dvars": float(np.percentile(dvars, 95)) if dvars.size else float("nan"),
        "mean_tsnr": float(np.nanmean(tsnr)) if tsnr.size else float("nan"),
        "median_tsnr": float(np.nanmedian(tsnr)) if tsnr.size else float("nan"),
        "_fd_series": diffs,
        "_dvars_series": dvars,
    }


def write_tsv(path: Path, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def main() -> int:
    args = parse_args()
    cfg = GateConfig(
        fd_threshold_mm=args.fd_threshold_mm,
        dvars_z_threshold=args.dvars_z_threshold,
        tsnr_min_median=args.tsnr_min_median,
    )

    manifest = Path(args.alignment_manifest).resolve()
    out_dir = Path(args.out_dir).resolve()
    if not manifest.is_file():
        raise SystemExit(f"ERROR: missing manifest: {manifest}")

    rows = list(csv.DictReader(manifest.open(), delimiter="\t"))
    if not rows:
        raise SystemExit("ERROR: alignment manifest has no data rows")

    gate_rows: List[Dict[str, object]] = []
    fmri_rows: List[Dict[str, object]] = []
    by_timepoint: Dict[str, List[bool]] = {}

    for r in rows:
        tp = r["timepoint"]
        role = r["role"]
        label = r["label"]
        src = Path(r["source_image"])
        struct = Path(r["structural_native"])
        aligned_struct = Path(r["aligned_struct"])
        aligned_template = Path(r["aligned_template"])
        aligned_allen = Path(r["aligned_allen"])

        missing = [str(p) for p in [src, struct, aligned_struct, aligned_template, aligned_allen] if not p.is_file()]
        if missing:
            gate_rows.append(
                {
                    "timepoint": tp,
                    "role": role,
                    "label": label,
                    "source_image": str(src),
                    "source_ndim": r.get("source_ndim", ""),
                    "coverage_frac": float("nan"),
                    "centroid_mm": float("nan"),
                    "ncc_struct": float("nan"),
                    "gate_pass": 0,
                    "gate_reasons": "missing_files:" + ",".join(missing),
                }
            )
            by_timepoint.setdefault(tp, []).append(False)
            continue

        src_data, src_aff = load_data(src)
        src3 = mean3d(src_data)
        struct_data, struct_aff = load_data(struct)
        struct3 = mean3d(struct_data)
        al_struct_data, al_struct_aff = load_data(aligned_struct)
        al_struct3 = mean3d(al_struct_data)

        struct_mask = mask_nonzero(struct3)
        mod_mask = mask_nonzero(al_struct3)
        overlap = struct_mask & mod_mask
        coverage_frac = float(np.sum(overlap) / np.sum(struct_mask)) if np.sum(struct_mask) > 0 else float("nan")

        c_struct = centroid_mm(struct_mask, struct_aff)
        c_mod = centroid_mm(mod_mask, al_struct_aff)
        if c_struct is None or c_mod is None:
            centroid_dist = float("nan")
        else:
            centroid_dist = float(np.linalg.norm(c_mod - c_struct))

        ncc_mask = overlap
        ncc_val = float("nan")
        if np.sum(ncc_mask) > 100:
            ncc_val = corrcoef_safe(struct3[ncc_mask], al_struct3[ncc_mask])

        finite_ok = bool(np.isfinite(al_struct3).all())
        reasons: List[str] = []
        if not finite_ok:
            reasons.append("nonfinite_values")

        gate_pass = True
        if role == "t2tse":
            if not np.isfinite(coverage_frac) or coverage_frac < cfg.t2_min_coverage:
                gate_pass = False
                reasons.append("low_coverage_t2")
            if not np.isfinite(centroid_dist) or centroid_dist > cfg.t2_max_centroid_mm:
                gate_pass = False
                reasons.append("high_centroid_dist_t2")
        elif role == "fmri":
            if not np.isfinite(coverage_frac) or coverage_frac < cfg.fmri_min_coverage:
                gate_pass = False
                reasons.append("low_coverage_fmri")
            if not np.isfinite(centroid_dist) or centroid_dist > cfg.fmri_max_centroid_mm:
                gate_pass = False
                reasons.append("high_centroid_dist_fmri")
        else:
            reasons.append("unknown_role")

        if not finite_ok:
            gate_pass = False

        gate_rows.append(
            {
                "timepoint": tp,
                "role": role,
                "label": label,
                "source_image": str(src),
                "source_ndim": r.get("source_ndim", ""),
                "coverage_frac": coverage_frac,
                "centroid_mm": centroid_dist,
                "ncc_struct": ncc_val,
                "gate_pass": int(gate_pass),
                "gate_reasons": ";".join(reasons),
            }
        )
        by_timepoint.setdefault(tp, []).append(gate_pass)

        # fMRI run-level metrics for 4D inputs
        try:
            src_ndim = int(float(r.get("source_ndim", "1")))
        except ValueError:
            src_ndim = 1

        if role == "fmri" and src_ndim > 1:
            al_tpl_data, al_tpl_aff = load_data(aligned_template)
            if al_tpl_data.ndim != 4:
                continue
            mean_vol = np.mean(al_tpl_data, axis=3)
            fmri_mask = mask_nonzero(mean_vol)
            met = compute_fd_dvars_tsnr(al_tpl_data, fmri_mask)
            fd_series = met.pop("_fd_series", np.array([], dtype=np.float32))
            dvars_series = met.pop("_dvars_series", np.array([], dtype=np.float32))

            dvars_z = robust_z(dvars_series) if dvars_series.size else np.array([], dtype=np.float32)
            censor_fd = fd_series > cfg.fd_threshold_mm if fd_series.size else np.array([], dtype=bool)
            censor_dvars = np.abs(dvars_z) > cfg.dvars_z_threshold if dvars_z.size else np.array([], dtype=bool)
            n_diff = max(fd_series.size, dvars_series.size, 1)
            n_use = min(fd_series.size, dvars_series.size) if fd_series.size and dvars_series.size else 0
            if n_use > 0:
                union = censor_fd[:n_use] | censor_dvars[:n_use]
                pct_union = 100.0 * float(np.mean(union))
            else:
                pct_union = float("nan")

            pct_fd = 100.0 * float(np.mean(censor_fd)) if censor_fd.size else float("nan")
            pct_dvars = 100.0 * float(np.mean(censor_dvars)) if censor_dvars.size else float("nan")
            tr_sec = nib.load(str(aligned_template)).header.get_zooms()[3]

            fmri_pass = True
            fmri_reasons: List[str] = []
            if np.isfinite(met["median_tsnr"]) and met["median_tsnr"] < cfg.tsnr_min_median:
                fmri_pass = False
                fmri_reasons.append("low_tsnr")
            if np.isfinite(pct_union) and pct_union > 40.0:
                fmri_pass = False
                fmri_reasons.append("high_censor_pct")

            fmri_rows.append(
                {
                    "timepoint": tp,
                    "label": label,
                    "aligned_template_4d": str(aligned_template),
                    "n_volumes": int(met["n_volumes"]),
                    "tr_sec": tr_sec,
                    "mean_fd_proxy_mm": met["mean_fd_proxy_mm"],
                    "p95_fd_proxy_mm": met["p95_fd_proxy_mm"],
                    "mean_dvars": met["mean_dvars"],
                    "p95_dvars": met["p95_dvars"],
                    "mean_tsnr": met["mean_tsnr"],
                    "median_tsnr": met["median_tsnr"],
                    "pct_censored_fd": pct_fd,
                    "pct_censored_dvars": pct_dvars,
                    "pct_censored_union": pct_union,
                    "fmri_run_pass": int(fmri_pass),
                    "fmri_run_reasons": ";".join(fmri_reasons),
                }
            )

    # Timepoint aggregate gate.
    tp_rows: List[Dict[str, object]] = []
    for tp, vals in sorted(by_timepoint.items()):
        vals_arr = np.asarray(vals, dtype=bool)
        tp_rows.append(
            {
                "timepoint": tp,
                "n_rows": int(vals_arr.size),
                "n_pass": int(np.sum(vals_arr)),
                "n_fail": int(np.sum(~vals_arr)),
                "timepoint_pass": int(np.all(vals_arr)),
            }
        )

    write_tsv(
        out_dir / "qc_gates.tsv",
        gate_rows,
        [
            "timepoint",
            "role",
            "label",
            "source_image",
            "source_ndim",
            "coverage_frac",
            "centroid_mm",
            "ncc_struct",
            "gate_pass",
            "gate_reasons",
        ],
    )
    write_tsv(
        out_dir / "timepoint_qc_gates.tsv",
        tp_rows,
        ["timepoint", "n_rows", "n_pass", "n_fail", "timepoint_pass"],
    )
    write_tsv(
        out_dir / "fmri_run_metrics.tsv",
        fmri_rows,
        [
            "timepoint",
            "label",
            "aligned_template_4d",
            "n_volumes",
            "tr_sec",
            "mean_fd_proxy_mm",
            "p95_fd_proxy_mm",
            "mean_dvars",
            "p95_dvars",
            "mean_tsnr",
            "median_tsnr",
            "pct_censored_fd",
            "pct_censored_dvars",
            "pct_censored_union",
            "fmri_run_pass",
            "fmri_run_reasons",
        ],
    )

    print(f"Wrote: {out_dir / 'qc_gates.tsv'}")
    print(f"Wrote: {out_dir / 'timepoint_qc_gates.tsv'}")
    print(f"Wrote: {out_dir / 'fmri_run_metrics.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

