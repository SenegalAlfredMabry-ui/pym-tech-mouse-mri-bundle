#!/usr/bin/env python3
"""Aggregate parcel stats, Jacobian-by-parcel, and QC metrics for mouse longitudinal runs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Tuple

import nibabel as nib
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build tables + QC report for mouse longitudinal outputs")
    p.add_argument("--manifest", required=True, help="TSV with columns: method,timepoint,image_allen,jacobian,roundtrip_parcels")
    p.add_argument("--allen-parcels", required=True, help="Allen parcel labels in Allen space")
    p.add_argument("--template-parcels", required=True, help="Allen parcel labels warped to template space")
    p.add_argument("--atlas-labels", default="", help="Optional labels TSV/CSV with columns parcel_id,label")
    p.add_argument("--out-dir", required=True, help="Directory for output tables/report")
    return p.parse_args()


def df_to_markdown(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._\n"
    cols = [str(c) for c in df.columns]
    lines = []
    lines.append("| " + " | ".join(cols) + " |")
    lines.append("| " + " | ".join(["---"] * len(cols)) + " |")
    for _, row in df.iterrows():
        vals = []
        for c in df.columns:
            v = row[c]
            if pd.isna(v):
                vals.append("")
            elif isinstance(v, float):
                vals.append(f"{v:.6g}")
            else:
                vals.append(str(v))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines) + "\n"


def load_labels(path: str) -> Dict[int, str]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        return {}
    sep = "," if p.suffix.lower() == ".csv" else "\t"
    df = pd.read_csv(p, sep=sep)
    id_col = None
    for c in df.columns:
        if c.lower() in {"parcel_id", "id", "label_id", "index"}:
            id_col = c
            break
    name_col = None
    for c in df.columns:
        if c.lower() in {"label", "name", "region", "parcel"}:
            name_col = c
            break
    if id_col is None or name_col is None:
        return {}
    out = {}
    for _, row in df.iterrows():
        try:
            out[int(row[id_col])] = str(row[name_col])
        except Exception:
            continue
    return out


def unique_labels(arr: np.ndarray) -> List[int]:
    vals = np.unique(arr.astype(np.int64))
    vals = vals[vals > 0]
    return [int(v) for v in vals.tolist()]


def dice(a: np.ndarray, b: np.ndarray) -> float:
    inter = np.logical_and(a, b).sum()
    s = a.sum() + b.sum()
    if s == 0:
        return float("nan")
    return float(2.0 * inter / s)


def centroid_mm(mask: np.ndarray, affine: np.ndarray) -> np.ndarray:
    ijk = np.argwhere(mask)
    if ijk.size == 0:
        return np.array([np.nan, np.nan, np.nan], dtype=float)
    c = ijk.mean(axis=0)
    xyz = nib.affines.apply_affine(affine, c)
    return np.asarray(xyz, dtype=float)


def parcel_means(data: np.ndarray, labels: np.ndarray, ids: List[int]) -> Dict[int, float]:
    out: Dict[int, float] = {}
    for pid in ids:
        m = labels == pid
        if not np.any(m):
            out[pid] = float("nan")
        else:
            out[pid] = float(np.nanmean(data[m]))
    return out


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(args.manifest, sep="\t")
    req_cols = {"method", "timepoint", "image_allen", "jacobian", "roundtrip_parcels"}
    missing = req_cols - set(manifest.columns)
    if missing:
        raise ValueError(f"Manifest missing columns: {sorted(missing)}")

    allen_img = nib.load(args.allen_parcels)
    allen_labels = np.asarray(np.rint(allen_img.get_fdata()), dtype=np.int64)
    tpl_img = nib.load(args.template_parcels)
    tpl_labels = np.asarray(np.rint(tpl_img.get_fdata()), dtype=np.int64)

    allen_ids = unique_labels(allen_labels)
    tpl_ids = unique_labels(tpl_labels)
    all_ids = sorted(set(allen_ids) | set(tpl_ids))

    label_names = load_labels(args.atlas_labels)

    parcel_rows: List[Dict[str, object]] = []
    jac_rows: List[Dict[str, object]] = []
    qc_rows: List[Dict[str, object]] = []

    for _, row in manifest.iterrows():
        method = str(row["method"])
        tp = str(row["timepoint"])

        image_allen = nib.load(str(row["image_allen"]))
        image_allen_data = np.asarray(image_allen.get_fdata(), dtype=np.float32)
        jac_img = nib.load(str(row["jacobian"]))
        jac_data = np.asarray(jac_img.get_fdata(), dtype=np.float32)
        rt_img = nib.load(str(row["roundtrip_parcels"]))
        rt_labels = np.asarray(np.rint(rt_img.get_fdata()), dtype=np.int64)

        # ROI intensity in Allen space
        pmeans = parcel_means(image_allen_data, allen_labels, all_ids)
        for pid, val in pmeans.items():
            parcel_rows.append(
                {
                    "method": method,
                    "timepoint": tp,
                    "parcel_id": pid,
                    "parcel_label": label_names.get(pid, f"parcel_{pid}"),
                    "mean_intensity": val,
                }
            )

        # Jacobian means in template-space parcels
        jmeans = parcel_means(jac_data, tpl_labels, all_ids)
        for pid, val in jmeans.items():
            jac_rows.append(
                {
                    "method": method,
                    "timepoint": tp,
                    "parcel_id": pid,
                    "parcel_label": label_names.get(pid, f"parcel_{pid}"),
                    "mean_jacobian": val,
                }
            )

        # QC: parcel dice + landmark (centroid) error
        per_dice = []
        per_dist = []
        for pid in all_ids:
            m_ref = tpl_labels == pid
            m_rt = rt_labels == pid
            d = dice(m_ref, m_rt)
            c_ref = centroid_mm(m_ref, tpl_img.affine)
            c_rt = centroid_mm(m_rt, tpl_img.affine)
            if np.any(np.isnan(c_ref)) or np.any(np.isnan(c_rt)):
                dist = float("nan")
            else:
                dist = float(np.linalg.norm(c_ref - c_rt))
            per_dice.append(d)
            per_dist.append(dist)

        qc_rows.append(
            {
                "method": method,
                "timepoint": tp,
                "mean_dice": float(np.nanmean(per_dice)),
                "median_dice": float(np.nanmedian(per_dice)),
                "mean_landmark_error_mm": float(np.nanmean(per_dist)),
                "median_landmark_error_mm": float(np.nanmedian(per_dist)),
            }
        )

    parcel_df = pd.DataFrame(parcel_rows)
    jac_df = pd.DataFrame(jac_rows)
    qc_df = pd.DataFrame(qc_rows)

    parcel_path = out_dir / "parcel_stats.tsv"
    jac_path = out_dir / "jacobian_by_parcel.tsv"
    qc_path = out_dir / "qc_landmarks_dice.tsv"

    parcel_df.to_csv(parcel_path, sep="\t", index=False)
    jac_df.to_csv(jac_path, sep="\t", index=False)
    qc_df.to_csv(qc_path, sep="\t", index=False)

    report_path = out_dir / "QC_report.md"
    with report_path.open("w") as f:
        f.write("# Mouse Longitudinal QC Report\n\n")
        f.write("## Summary by Method\n\n")
        if not qc_df.empty:
            summary = (
                qc_df.groupby("method", dropna=False)
                .agg(
                    n_timepoints=("timepoint", "count"),
                    mean_dice=("mean_dice", "mean"),
                    mean_landmark_error_mm=("mean_landmark_error_mm", "mean"),
                )
                .reset_index()
            )
            f.write(df_to_markdown(summary))
            f.write("\n\n")
        else:
            f.write("No QC rows were generated.\n\n")

        f.write("## Timepoint QC\n\n")
        if qc_df.empty:
            f.write("No timepoint metrics available.\n")
        else:
            f.write(df_to_markdown(qc_df.sort_values(["method", "timepoint"])))
            f.write("\n")

    print(f"Wrote: {parcel_path}")
    print(f"Wrote: {jac_path}")
    print(f"Wrote: {qc_path}")
    print(f"Wrote: {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
