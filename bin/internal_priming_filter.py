#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def load_blacklist(path: Path):
    intervals = []
    if not path.exists() or path.name == "NO_FILE" or path.stat().st_size == 0:
        return intervals
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            intervals.append((row[0], int(row[1]) + 1, int(row[2])))
    return intervals


def overlaps(intervals, chrom, position):
    for blacklist_chrom, start, end in intervals:
        if chrom == blacklist_chrom and start <= position <= end:
            return True
    return False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--site-catalog", required=True)
    parser.add_argument("--genome-fasta", required=True)
    parser.add_argument("--priming-blacklist", required=True)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--out-metrics", required=True)
    args = parser.parse_args()

    blacklist = load_blacklist(Path(args.priming_blacklist))
    flagged = 0
    total = 0
    out_rows = []

    with Path(args.site_catalog).open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            total += 1
            position = int(row.get("end") or row.get("start") or 0)
            if overlaps(blacklist, row.get("chrom", ""), position):
                row["priming_flag"] = "blacklist_overlap"
                flagged += 1
            else:
                row["priming_flag"] = row.get("priming_flag") or "pass"
            out_rows.append(row)

    with Path(args.out_tsv).open("w", newline="", encoding="utf-8") as handle:
        fieldnames = ["site_id", "gene_id", "chrom", "start", "end", "strand", "site_class", "site_source", "priming_flag"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    with Path(args.out_metrics).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["total_sites", "flagged_sites"])
        writer.writerow([total, flagged])


if __name__ == "__main__":
    main()
