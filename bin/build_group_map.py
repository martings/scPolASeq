#!/usr/bin/env python3

import argparse
import csv
from collections import Counter
from pathlib import Path


CANONICAL_COLUMNS = [
    "sample_id",
    "library_id",
    "barcode_raw",
    "barcode_corrected",
    "cell_id",
    "cluster_id",
    "cell_type",
    "condition",
    "batch",
    "label_source",
]


def canonical_cell_id(sample_id: str, library_id: str, barcode: str) -> str:
    parts = [sample_id]
    if library_id:
        parts.append(library_id)
    parts.append(barcode)
    return ":".join(parts)


def load_rows(paths):
    rows = []
    for path in paths:
        candidate = Path(path)
        if not candidate.exists() or candidate.stat().st_size == 0:
            continue
        with candidate.open("r", newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            rows.extend(reader)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotation-files", nargs="+", required=True)
    parser.add_argument("--grouping-levels", required=True)
    parser.add_argument("--out-cell-annotations", required=True)
    parser.add_argument("--out-group-map", required=True)
    parser.add_argument("--out-group-summary", required=True)
    args = parser.parse_args()

    rows = load_rows(args.annotation_files)
    rows.sort(key=lambda row: ((row.get("sample_id") or ""), (row.get("library_id") or ""), (row.get("barcode_corrected") or "")))

    with Path(args.out_cell_annotations).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CANONICAL_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    grouping_levels = [level.strip() for level in args.grouping_levels.split(",") if level.strip()]
    group_rows = []
    summary = Counter()
    for row in rows:
        sample_id = row.get("sample_id", "")
        library_id = row.get("library_id", "")
        barcode = row.get("barcode_corrected") or row.get("barcode_raw") or ""
        condition = row.get("condition", "") or "NA"
        cluster_id = row.get("cluster_id", "")
        cell_type = row.get("cell_type", "") or "unlabeled"

        for level in grouping_levels:
            group_id = None
            if level == "cell":
                group_id = row.get("cell_id") or canonical_cell_id(sample_id, library_id, barcode)
            elif level == "cluster" and cluster_id:
                group_id = cluster_id
            elif level == "cell_type" and cell_type:
                group_id = cell_type
            elif level == "sample_condition_celltype":
                group_id = f"{sample_id}__{condition}__{cell_type}"

            if group_id:
                group_rows.append(
                    {
                        "sample_id": sample_id,
                        "library_id": library_id,
                        "barcode_corrected": barcode,
                        "group_level": level,
                        "group_id": group_id,
                    }
                )
                summary[(level, group_id)] += 1

    with Path(args.out_group_map).open("w", newline="", encoding="utf-8") as handle:
        fieldnames = ["sample_id", "library_id", "barcode_corrected", "group_level", "group_id"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(group_rows)

    with Path(args.out_group_summary).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["group_level", "group_id", "n_barcodes"])
        for (level, group_id), count in sorted(summary.items()):
            writer.writerow([level, group_id, count])


if __name__ == "__main__":
    main()
