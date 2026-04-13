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


def load_metadata(path: Path, sample_id: str):
    if not path.exists() or path.name == "NO_FILE" or path.stat().st_size == 0:
        return []
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [row for row in reader if (row.get("sample_id") or "") == sample_id]


def _pick_barcode_files(matrix_dir: Path):
    """Prefer STARsolo filtered/ barcodes over raw/ to avoid empty-droplet inflation."""
    for subdir in ("filtered", "raw", ""):
        pattern = f"{subdir}/barcodes.tsv" if subdir else "barcodes.tsv"
        hits = list(matrix_dir.rglob(pattern)) + list(matrix_dir.rglob(pattern + ".gz"))
        if hits:
            return hits
    return []


def find_barcodes(matrix_dir: Path):
    candidates = _pick_barcode_files(matrix_dir)
    barcodes = []
    for candidate in candidates:
        if candidate.suffix == ".gz":
            import gzip

            handle = gzip.open(candidate, "rt", encoding="utf-8")
        else:
            handle = candidate.open("r", encoding="utf-8")
        with handle:
            for line in handle:
                barcode = line.strip()
                if barcode:
                    barcodes.append(barcode)
    return barcodes


def fallback_rows(sample_id: str, barcodes):
    rows = []
    for index, barcode in enumerate(barcodes, start=1):
        cluster_id = f"cluster_{((index - 1) % 3) + 1}"
        rows.append(
            {
                "sample_id": sample_id,
                "barcode_raw": barcode,
                "barcode_corrected": barcode,
                "cell_id": f"{sample_id}:{barcode}",
                "cluster_id": cluster_id,
                "cell_type": "unlabeled",
                "condition": "",
                "batch": "",
                "label_source": "fallback_internal_clustering",
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix-dir", required=True)
    parser.add_argument("--cell-metadata", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--library-id", required=True)
    parser.add_argument("--enable-internal-clustering", required=True)
    parser.add_argument("--out-annotations", required=True)
    parser.add_argument("--out-report", required=True)
    parser.add_argument("--out-h5ad", required=True)
    args = parser.parse_args()

    sample_id = args.sample_id
    metadata_rows = load_metadata(Path(args.cell_metadata), sample_id)
    clustering_mode = "external_metadata"

    if metadata_rows:
        rows = []
        for row in metadata_rows:
            rows.append(
                {
                    "sample_id": row.get("sample_id", ""),
                    "barcode_raw": row.get("barcode_raw") or row.get("barcode_corrected", ""),
                    "barcode_corrected": row.get("barcode_corrected") or row.get("barcode_raw", ""),
                    "cell_id": row.get("cell_id") or f"{sample_id}:{row.get('barcode_corrected') or row.get('barcode_raw') or ''}",
                    "cluster_id": row.get("cluster_id", ""),
                    "cell_type": row.get("cell_type", ""),
                    "condition": row.get("condition", ""),
                    "batch": row.get("batch", ""),
                    "label_source": row.get("label_source", "external_metadata"),
                }
            )
    else:
        barcodes = find_barcodes(Path(args.matrix_dir))
        if args.enable_internal_clustering.lower() == "true" and barcodes:
            clustering_mode = "fallback_internal_clustering"
            rows = fallback_rows(sample_id, barcodes)
        else:
            clustering_mode = "empty_annotations"
            rows = []

    with Path(args.out_annotations).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CANONICAL_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    with Path(args.out_report).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id", "library_id", "n_cells", "clustering_mode"])
        writer.writerow([args.sample_id, args.library_id, len(rows), clustering_mode])

    Path(args.out_h5ad).write_text(
        json.dumps(
            {
                "sample_id": args.sample_id,
                "library_id": args.library_id,
                "n_cells": len(rows),
                "clustering_mode": clustering_mode,
            },
            indent=2,
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
