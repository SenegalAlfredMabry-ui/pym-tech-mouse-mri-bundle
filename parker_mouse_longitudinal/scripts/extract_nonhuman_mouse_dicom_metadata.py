#!/usr/bin/env python3
"""Extract key DICOM metadata for nonhuman mouse scans.

Designed for environments without pydicom by parsing early DICOM headers only.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import struct
from collections import defaultdict
from pathlib import Path


TARGET_TAGS = {
    (0x0010, 0x0010): "patient_name",
    (0x0010, 0x0020): "patient_id",
    (0x0008, 0x0020): "study_date",
    (0x0008, 0x0030): "study_time",
    (0x0008, 0x0021): "series_date",
    (0x0008, 0x0031): "series_time",
    (0x0008, 0x0022): "acq_date",
    (0x0008, 0x0032): "acq_time",
    (0x0020, 0x000D): "study_uid",
    (0x0020, 0x000E): "series_uid",
    (0x0020, 0x0011): "series_number",
    (0x0008, 0x103E): "series_description",
    (0x0018, 0x1030): "protocol_name",
    (0x0018, 0x0024): "sequence_name",
    (0x0008, 0x0070): "manufacturer",
    (0x0008, 0x1090): "scanner_model",
    (0x0018, 0x0087): "field_strength_T",
    (0x0008, 0x1010): "station_name",
    (0x0028, 0x0010): "rows",
    (0x0028, 0x0011): "cols",
    (0x0028, 0x0030): "pixel_spacing_mm",
    (0x0018, 0x0088): "spacing_between_slices_mm",
    (0x0018, 0x0050): "slice_thickness_mm",
}

LONG_VR = {
    b"OB",
    b"OW",
    b"OF",
    b"SQ",
    b"UT",
    b"UN",
    b"UR",
    b"OD",
    b"OL",
    b"UC",
    b"SV",
    b"UV",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract nonhuman mouse DICOM metadata and session chronology."
    )
    parser.add_argument(
        "--root",
        default="${SNHEART_ROOT:-/path/to/SNHEART}/MRI_data/Raw/nonhuman_mouse",
        help="Root directory containing extracted nonhuman mouse data.",
    )
    parser.add_argument(
        "--dicom-subdir",
        default="Schmitz^Trip-KI-Char Imaging",
        help="Subdirectory under root that contains DICOM folders.",
    )
    parser.add_argument(
        "--max-header-bytes",
        type=int,
        default=4 * 1024 * 1024,
        help="Max bytes to read from each DICOM file while parsing headers.",
    )
    return parser.parse_args()


def clean_text(raw: bytes) -> str:
    return raw.decode("latin1", errors="ignore").strip("\x00 ").strip()


def parse_date_time(date_s: str, time_s: str) -> str:
    if not date_s:
        return ""
    time_part = (time_s or "").split(".")[0]
    time_part = (time_part + "000000")[:6]
    try:
        return dt.datetime.strptime(date_s + time_part, "%Y%m%d%H%M%S").isoformat()
    except ValueError:
        return ""


def read_tag(data: bytes, off: int, endian: str):
    if off + 4 > len(data):
        return None, None, off
    group, elem = struct.unpack(endian + "HH", data[off : off + 4])
    return group, elem, off + 4


def parse_element(data: bytes, off: int, explicit: bool, endian: str):
    group, elem, off = read_tag(data, off, endian)
    if group is None:
        return None

    if explicit:
        if off + 2 > len(data):
            return None
        vr = data[off : off + 2]
        off += 2
        if vr in LONG_VR:
            if off + 6 > len(data):
                return None
            off += 2  # reserved bytes
            length = struct.unpack(endian + "I", data[off : off + 4])[0]
            off += 4
        else:
            if off + 2 > len(data):
                return None
            length = struct.unpack(endian + "H", data[off : off + 2])[0]
            off += 2
    else:
        vr = None
        if off + 4 > len(data):
            return None
        length = struct.unpack(endian + "I", data[off : off + 4])[0]
        off += 4

    if length == 0xFFFFFFFF:
        return {
            "tag": (group, elem),
            "vr": vr,
            "length": length,
            "value": b"",
            "next": off,
            "undef": True,
        }

    end = off + length
    if end > len(data):
        return None
    return {
        "tag": (group, elem),
        "vr": vr,
        "length": length,
        "value": data[off:end],
        "next": end,
        "undef": False,
    }


def parse_dicom_header(path: Path, max_header_bytes: int) -> dict[str, str]:
    out = {v: "" for v in TARGET_TAGS.values()}
    with path.open("rb") as f:
        raw = f.read(max_header_bytes)

    off = 128
    if len(raw) < 132 or raw[128:132] != b"DICM":
        off = 0
    else:
        off = 132

    # File meta information is explicit little-endian.
    ts_uid = ""
    m_off = off
    while m_off < len(raw):
        tag = parse_element(raw, m_off, explicit=True, endian="<")
        if not tag:
            break
        group, elem = tag["tag"]
        if group != 0x0002:
            break
        if (group, elem) == (0x0002, 0x0010):
            ts_uid = clean_text(tag["value"])
        m_off = tag["next"]

    if ts_uid == "1.2.840.10008.1.2":
        explicit, endian = False, "<"  # implicit little
    elif ts_uid == "1.2.840.10008.1.2.2":
        explicit, endian = True, ">"  # explicit big
    else:
        explicit, endian = True, "<"  # explicit little default

    off = m_off
    seen = set()
    while off < len(raw):
        tag = parse_element(raw, off, explicit=explicit, endian=endian)
        if not tag:
            break
        off = tag["next"]
        if tag["undef"]:
            break

        t = tag["tag"]
        if t == (0x7FE0, 0x0010):
            break
        if t in TARGET_TAGS:
            key = TARGET_TAGS[t]
            val = tag["value"]
            if key in {"rows", "cols"} and len(val) >= 2:
                try:
                    out[key] = str(struct.unpack(endian + "H", val[:2])[0])
                except Exception:
                    out[key] = clean_text(val)
            else:
                out[key] = clean_text(val)
            seen.add(key)
            if len(seen) == len(TARGET_TAGS):
                break
    return out


def int_or_large(value: str) -> int:
    return int(value) if value.isdigit() else 10**9


def main() -> None:
    args = parse_args()
    root = Path(args.root).resolve()
    base = root / args.dicom_subdir
    if not base.exists():
        raise SystemExit(f"DICOM subdir not found: {base}")

    series_csv = root / "mouse_dicom_series_inventory.csv"
    study_csv = root / "mouse_dicom_study_inventory.csv"
    summary_txt = root / "mouse_scan_summary.txt"

    dcms = sorted(base.rglob("*.dcm"))
    if not dcms:
        raise SystemExit(f"No DICOM files found under: {base}")

    series_rows: list[dict[str, str]] = []
    for dcm in dcms:
        md = parse_dicom_header(dcm, args.max_header_bytes)
        rel = dcm.relative_to(root)
        parts = rel.parts
        top_mouse = parts[1] if len(parts) > 1 else ""
        exam_folder = parts[2] if len(parts) > 2 else ""
        series_folder = parts[3] if len(parts) > 3 else ""

        acq_dt = parse_date_time(
            md["acq_date"] or md["series_date"] or md["study_date"],
            md["acq_time"] or md["series_time"] or md["study_time"],
        )
        study_dt = parse_date_time(md["study_date"], md["study_time"])
        series_dt = parse_date_time(
            md["series_date"] or md["study_date"],
            md["series_time"] or md["acq_time"] or md["study_time"],
        )

        series_rows.append(
            {
                "top_mouse_folder": top_mouse,
                "exam_folder": exam_folder,
                "series_folder": series_folder,
                "dicom_path": str(rel),
                "patient_name": md["patient_name"],
                "patient_id": md["patient_id"],
                "study_instance_uid": md["study_uid"],
                "series_instance_uid": md["series_uid"],
                "study_date": md["study_date"],
                "study_time": md["study_time"],
                "series_date": md["series_date"],
                "series_time": md["series_time"],
                "acquisition_date": md["acq_date"],
                "acquisition_time": md["acq_time"],
                "acq_datetime": acq_dt,
                "study_datetime": study_dt,
                "series_datetime": series_dt,
                "series_number": md["series_number"],
                "series_description": md["series_description"],
                "protocol_name": md["protocol_name"],
                "sequence_name": md["sequence_name"],
                "manufacturer": md["manufacturer"],
                "scanner_model": md["scanner_model"],
                "station_name": md["station_name"],
                "field_strength_T": md["field_strength_T"],
                "rows": md["rows"],
                "cols": md["cols"],
                "pixel_spacing_mm": md["pixel_spacing_mm"],
                "spacing_between_slices_mm": md["spacing_between_slices_mm"],
                "slice_thickness_mm": md["slice_thickness_mm"],
            }
        )

    series_rows_sorted = sorted(
        series_rows,
        key=lambda r: (
            r["top_mouse_folder"],
            r["exam_folder"],
            int_or_large(r["series_number"]),
        ),
    )
    with series_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(series_rows_sorted[0].keys()))
        writer.writeheader()
        writer.writerows(series_rows_sorted)

    study_map: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in series_rows:
        study_map[(row["top_mouse_folder"], row["exam_folder"])].append(row)

    study_rows: list[dict[str, str]] = []
    for (top_mouse, exam_folder), rows in study_map.items():
        datetimes = sorted(
            {
                r["acq_datetime"] or r["series_datetime"] or r["study_datetime"]
                for r in rows
                if (r["acq_datetime"] or r["series_datetime"] or r["study_datetime"])
            }
        )
        series_nums = sorted(
            {int(r["series_number"]) for r in rows if r["series_number"].isdigit()}
        )
        descs = sorted({r["series_description"] for r in rows if r["series_description"]})
        protocols = sorted({r["protocol_name"] for r in rows if r["protocol_name"]})

        study_rows.append(
            {
                "top_mouse_folder": top_mouse,
                "exam_folder": exam_folder,
                "patient_name": rows[0]["patient_name"],
                "patient_id": rows[0]["patient_id"],
                "study_instance_uids": ";".join(
                    sorted({r["study_instance_uid"] for r in rows if r["study_instance_uid"]})
                ),
                "n_series": str(len(rows)),
                "first_datetime": datetimes[0] if datetimes else "",
                "last_datetime": datetimes[-1] if datetimes else "",
                "series_numbers": ",".join(str(x) for x in series_nums),
                "series_descriptions": " | ".join(descs),
                "protocol_names": " | ".join(protocols),
                "manufacturer": rows[0]["manufacturer"],
                "scanner_model": rows[0]["scanner_model"],
                "station_name": rows[0]["station_name"],
                "field_strength_T": rows[0]["field_strength_T"],
                "pixel_spacing_mm": rows[0]["pixel_spacing_mm"],
                "slice_thickness_mm": rows[0]["slice_thickness_mm"],
            }
        )

    study_rows_sorted = sorted(study_rows, key=lambda r: r["first_datetime"] or "9999")
    with study_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(study_rows_sorted[0].keys()))
        writer.writeheader()
        writer.writerows(study_rows_sorted)

    with summary_txt.open("w") as f:
        f.write("Nonhuman mouse MRI DICOM inventory\n")
        f.write(f"Root: {base}\n")
        f.write(f"Total DICOM files parsed: {len(series_rows)}\n")
        f.write(f"Total top-level exam folders: {len(study_rows_sorted)}\n\n")
        if study_rows_sorted:
            first = study_rows_sorted[0]
            f.write("Earliest exam by metadata:\n")
            f.write(
                f"  {first['top_mouse_folder']} / {first['exam_folder']} / "
                f"first_datetime={first['first_datetime']}\n\n"
            )
        f.write("Chronological exam order:\n")
        for idx, row in enumerate(study_rows_sorted, start=1):
            f.write(
                f"{idx}. {row['top_mouse_folder']} | {row['exam_folder']} | "
                f"first={row['first_datetime']} | n_series={row['n_series']}\n"
            )
        f.write("\nImportant notes:\n")
        f.write("- Data contains five top-level mouse IDs under the extracted OneDrive package.\n")
        f.write("- Series order within exam uses SeriesNumber and acquisition/series times.\n")
        f.write(
            "- If multiple exams share date and lack distinct times, exact first-vs-second "
            "ordering cannot be proven from current headers alone.\n"
        )

    print(f"WROTE {series_csv}")
    print(f"WROTE {study_csv}")
    print(f"WROTE {summary_txt}")
    if study_rows_sorted:
        first = study_rows_sorted[0]
        print(
            "EARLIEST_EXAM "
            f"{first['top_mouse_folder']} {first['exam_folder']} {first['first_datetime']}"
        )


if __name__ == "__main__":
    main()
