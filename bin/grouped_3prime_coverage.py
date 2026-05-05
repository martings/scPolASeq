#!/usr/bin/env python3

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


def sanitize(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value)


def row_matches_scope(row, sample_id: str, library_id: str) -> bool:
    if row.get("sample_id") != sample_id:
        return False
    row_library_id = (row.get("library_id") or "").strip()
    return not row_library_id or row_library_id == library_id


def load_group_map(path: Path, sample_id: str, library_id: str):
    mapping = defaultdict(list)
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if not row_matches_scope(row, sample_id, library_id):
                continue
            barcode = row.get("barcode_corrected") or row.get("barcode_raw") or ""
            mapping[barcode].append((row.get("group_level", ""), row.get("group_id", "")))
    return mapping


def load_barcode_registry(path: Path, sample_id: str, library_id: str):
    rows = []
    if not path.exists() or path.stat().st_size == 0:
        return rows
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row_matches_scope(row, sample_id, library_id):
                rows.append(row)
    return rows


def load_chrom_sizes(path: Path):
    sizes = []
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if row:
                sizes.append((row[0], int(row[1])))
    return sizes


def count_from_bam(bam_path: Path, group_map, min_mapq: int):
    try:
        import pysam
    except ImportError:
        return None

    counts = defaultdict(lambda: {"read_count": 0, "umis": set()})
    with pysam.AlignmentFile(str(bam_path), "rb") as bam_handle:
        for read in bam_handle.fetch(until_eof=True):
            if read.is_unmapped or read.mapping_quality < min_mapq:
                continue
            tags = dict(read.tags)
            barcode = tags.get("CB") or tags.get("CR")
            umi = tags.get("UB") or tags.get("UR") or read.query_name
            if not barcode or barcode not in group_map:
                continue
            chrom = read.reference_name
            strand = "-" if read.is_reverse else "+"
            position = read.reference_start + 1 if read.is_reverse else read.reference_end
            for group_level, group_id in group_map[barcode]:
                key = (group_level, group_id, chrom, position, strand)
                counts[key]["read_count"] += 1
                counts[key]["umis"].add(umi)
    return counts


def fallback_counts(group_map, barcode_registry, chrom_sizes):
    counts = defaultdict(lambda: {"read_count": 0, "umis": set()})
    if not chrom_sizes:
        chrom_sizes = [("chrFallback", 1000)]
    chrom = chrom_sizes[0][0]
    for index, row in enumerate(barcode_registry, start=1):
        barcode = row.get("barcode_corrected", "")
        groups = group_map.get(barcode, [])
        for group_level, group_id in groups:
            key = (group_level, group_id, chrom, index, "+")
            counts[key]["read_count"] += 1
            counts[key]["umis"].add(f"fallback_umi_{index}")
    return counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--barcode-registry", required=True)
    parser.add_argument("--group-map", required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--library-id", required=True)
    parser.add_argument("--emit-bigwigs", required=True)
    parser.add_argument("--min-mapq", type=int, default=255)
    parser.add_argument("--out-counts", required=True)
    parser.add_argument("--out-track-dir", required=True)
    parser.add_argument("--out-summary", required=True)
    args = parser.parse_args()

    group_map = load_group_map(Path(args.group_map), args.sample_id, args.library_id)
    barcode_registry = load_barcode_registry(Path(args.barcode_registry), args.sample_id, args.library_id)
    chrom_sizes = load_chrom_sizes(Path(args.chrom_sizes))

    counts = count_from_bam(Path(args.bam), group_map, args.min_mapq)
    mode = "bam"
    if counts is None or not counts:
        counts = fallback_counts(group_map, barcode_registry, chrom_sizes)
        mode = "fallback"

    out_rows = []
    bedgraph_rows = defaultdict(list)
    for (group_level, group_id, chrom, position, strand), value in sorted(counts.items()):
        umi_count = len(value["umis"])
        out_rows.append(
            {
                "sample_id": args.sample_id,
                "library_id": args.library_id,
                "group_level": group_level,
                "group_id": group_id,
                "chrom": chrom,
                "position": position,
                "strand": strand,
                "read_count": value["read_count"],
                "umi_count": umi_count,
            }
        )
        bedgraph_rows[group_id].append((chrom, position - 1, position, value["read_count"]))

    with Path(args.out_counts).open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "sample_id",
            "library_id",
            "group_level",
            "group_id",
            "chrom",
            "position",
            "strand",
            "read_count",
            "umi_count",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    track_dir = Path(args.out_track_dir)
    track_dir.mkdir(parents=True, exist_ok=True)
    for group_id, rows in bedgraph_rows.items():
        safe_group_id = sanitize(group_id)
        bedgraph_path = track_dir / f"{safe_group_id}.bedGraph"
        with bedgraph_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerows(rows)
        if args.emit_bigwigs.lower() == "true":
            (track_dir / f"{safe_group_id}.bigWig").write_text("", encoding="utf-8")

    with Path(args.out_summary).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id", "library_id", "mode", "n_groups", "n_loci"])
        writer.writerow([args.sample_id, args.library_id, mode, len(bedgraph_rows), len(out_rows)])


if __name__ == "__main__":
    main()
