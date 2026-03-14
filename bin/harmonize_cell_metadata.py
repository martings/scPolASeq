#!/usr/bin/env python3

import argparse
import csv
import json
from pathlib import Path


CANONICAL_COLUMNS = [
    "sample_id",
    "barcode_raw",
    "barcode_corrected",
    "cell_id",
    "cluster_id",
    "cell_type",
    "condition",
    "batch",
    "label_source",
]


def load_table(path: Path):
    if not path.exists() or path.name == "NO_FILE" or path.stat().st_size == 0:
        return []
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def barcode_key(row):
    sample_id = row.get("sample_id", "")
    barcode = row.get("barcode_corrected") or row.get("barcode_raw") or row.get("barcode") or ""
    return sample_id, barcode


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cell-metadata", required=True)
    parser.add_argument("--cluster-assignments", required=True)
    parser.add_argument("--cell-type-labels", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args()

    metadata_rows = load_table(Path(args.cell_metadata))
    cluster_rows = load_table(Path(args.cluster_assignments))
    cell_type_rows = load_table(Path(args.cell_type_labels))

    by_key = {}
    for row in metadata_rows:
        key = barcode_key(row)
        by_key[key] = {
            "sample_id": row.get("sample_id", ""),
            "barcode_raw": row.get("barcode_raw") or row.get("barcode") or row.get("barcode_corrected", ""),
            "barcode_corrected": row.get("barcode_corrected") or row.get("barcode") or row.get("barcode_raw", ""),
            "cell_id": row.get("cell_id") or f"{row.get('sample_id', '')}:{row.get('barcode_corrected') or row.get('barcode_raw') or row.get('barcode', '')}",
            "cluster_id": row.get("cluster_id", ""),
            "cell_type": row.get("cell_type", ""),
            "condition": row.get("condition", ""),
            "batch": row.get("batch", ""),
            "label_source": row.get("label_source", "external_metadata"),
        }

    for row in cluster_rows:
        key = barcode_key(row)
        existing = by_key.setdefault(
            key,
            {
                "sample_id": row.get("sample_id", ""),
                "barcode_raw": row.get("barcode_raw") or row.get("barcode") or row.get("barcode_corrected", ""),
                "barcode_corrected": row.get("barcode_corrected") or row.get("barcode") or row.get("barcode_raw", ""),
                "cell_id": f"{row.get('sample_id', '')}:{row.get('barcode_corrected') or row.get('barcode_raw') or row.get('barcode', '')}",
                "cluster_id": "",
                "cell_type": "",
                "condition": row.get("condition", ""),
                "batch": row.get("batch", ""),
                "label_source": "cluster_assignments",
            },
        )
        existing["cluster_id"] = row.get("cluster_id") or row.get("cluster") or existing["cluster_id"]
        existing["label_source"] = "cluster_assignments"

    for row in cell_type_rows:
        key = barcode_key(row)
        existing = by_key.setdefault(
            key,
            {
                "sample_id": row.get("sample_id", ""),
                "barcode_raw": row.get("barcode_raw") or row.get("barcode") or row.get("barcode_corrected", ""),
                "barcode_corrected": row.get("barcode_corrected") or row.get("barcode") or row.get("barcode_raw", ""),
                "cell_id": f"{row.get('sample_id', '')}:{row.get('barcode_corrected') or row.get('barcode_raw') or row.get('barcode', '')}",
                "cluster_id": "",
                "cell_type": "",
                "condition": row.get("condition", ""),
                "batch": row.get("batch", ""),
                "label_source": "cell_type_labels",
            },
        )
        existing["cell_type"] = row.get("cell_type") or row.get("label") or existing["cell_type"]
        if not existing["label_source"]:
            existing["label_source"] = "cell_type_labels"

    out_rows = sorted(by_key.values(), key=lambda row: (row["sample_id"], row["barcode_corrected"]))
    with Path(args.out).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CANONICAL_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    manifest = {
        "cell_metadata_rows": len(metadata_rows),
        "cluster_assignment_rows": len(cluster_rows),
        "cell_type_rows": len(cell_type_rows),
        "harmonized_rows": len(out_rows),
    }
    Path(args.manifest).write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
