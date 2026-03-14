#!/usr/bin/env python3

import argparse
import csv
import gzip
from pathlib import Path


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def write_rows(rows, out_path: Path):
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = ["sample_id", "library_id", "barcode_raw", "barcode_corrected", "source"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def extract_from_solo(solo_dir: Path, sample_id: str, library_id: str):
    candidates = list(solo_dir.rglob("barcodes.tsv")) + list(solo_dir.rglob("barcodes.tsv.gz"))
    rows = []
    seen = set()
    for candidate in candidates:
        with open_text(candidate) as handle:
            for line in handle:
                barcode = line.strip()
                if not barcode or barcode in seen:
                    continue
                seen.add(barcode)
                rows.append(
                    {
                        "sample_id": sample_id,
                        "library_id": library_id,
                        "barcode_raw": barcode,
                        "barcode_corrected": barcode,
                        "source": "starsolo",
                    }
                )
    return rows


def extract_from_bam(bam_path: Path, sample_id: str, library_id: str):
    try:
        import pysam
    except ImportError:
        return []

    rows = []
    seen = set()
    with pysam.AlignmentFile(str(bam_path), "rb") as bam_handle:
        for read in bam_handle.fetch(until_eof=True):
            tags = dict(read.tags)
            barcode = tags.get("CB") or tags.get("CR")
            if not barcode or barcode in seen:
                continue
            seen.add(barcode)
            rows.append(
                {
                    "sample_id": sample_id,
                    "library_id": library_id,
                    "barcode_raw": tags.get("CR", barcode),
                    "barcode_corrected": barcode,
                    "source": "bam_tags",
                }
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--library", required=True)
    parser.add_argument("--solo-dir")
    parser.add_argument("--bam")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    rows = []
    if args.solo_dir:
        rows = extract_from_solo(Path(args.solo_dir), args.sample, args.library)
    elif args.bam:
        rows = extract_from_bam(Path(args.bam), args.sample, args.library)
    write_rows(rows, Path(args.out))


if __name__ == "__main__":
    main()
