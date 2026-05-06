#!/usr/bin/env python3
"""
Convert SCAPTURE evaluated BED15 peak files to the normalized site_catalog TSV
used by the scPolASeq APA pipeline.

SCAPTURE BED15 columns:
  0-5:   chrom, start, end, name, score, strand  (standard BED6)
  6-11:  thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
  12:    number of reference poly(A) sites supporting the PAS
  13:    DeepPASS prediction label (positive/negative) in current SCAPTURE
  14:    DeepPASS prediction score (0–1) in current SCAPTURE

Output site_catalog columns:
  site_id, gene_id, chrom, start, end, strand, site_class, site_source, priming_flag
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

FIELDNAMES = [
    "site_id", "gene_id", "chrom", "start", "end", "strand",
    "site_class", "site_source", "priming_flag",
]


def _gene_from_name(name: str) -> str:
    """Extract gene symbol from the SCAPTURE BED name field.

    SCAPTURE encodes peak names in several observed formats:
      - "GeneSymbol:chrom:start:end:strand"
      - "GeneSymbol_peak_N"
      - "GeneSymbol|..."
    Fall back to the full name if no delimiter is recognised.
    """
    for sep in (":", "_peak", "|"):
        if sep in name:
            return name.split(sep)[0]
    return name


def _deeppass_score(cols: list[str]) -> float | None:
    """Return a DeepPASS score from known SCAPTURE BED column layouts.

    Current SCAPTURE (1.0): col 13 = "positive"/"negative" label,
    col 14 = 200 bp genomic sequence. Map the label to 1.0 / 0.0 so the
    downstream threshold filter works uniformly.
    Legacy SCAPTURE: col 13 or 14 may be a numeric probability (0–1).
    """
    if len(cols) > 13:
        label = cols[13].strip().lower()
        if label == "positive":
            return 1.0
        if label == "negative":
            return 0.0
    # Legacy: numeric probability in col 14 or 13
    for idx in (14, 13):
        if len(cols) <= idx:
            continue
        try:
            score = float(cols[idx])
        except ValueError:
            continue
        if 0.0 <= score <= 1.0:
            return score
    return None


def parse_bed15(path: Path, site_class: str, threshold: float):
    """Yield site_catalog row dicts from one SCAPTURE evaluated BED15 file."""
    if not path.exists() or path.stat().st_size == 0:
        return
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            cols = line.split("\t")
            if len(cols) < 14:
                continue
            try:
                start = int(cols[1])
                end = int(cols[2])
            except ValueError:
                continue
            chrom  = cols[0]
            name   = cols[3]
            strand = cols[5] if cols[5] in ("+", "-") else "."
            deeppass = _deeppass_score(cols)
            if deeppass is None:
                continue
            if deeppass < threshold:
                continue
            gene_id = _gene_from_name(name)
            site_id = f"{gene_id}:{chrom}:{start}-{end}:{strand}"
            yield {
                "site_id":      site_id,
                "gene_id":      gene_id,
                "chrom":        chrom,
                "start":        start,
                "end":          end,
                "strand":       strand,
                "site_class":   site_class,
                "site_source":  "scapture",
                "priming_flag": "unchecked",
            }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exonic",    required=True, help="Exonic evaluated BED15")
    parser.add_argument("--intronic",  required=True, help="Intronic evaluated BED15")
    parser.add_argument("--extended",  required=True, help="3′-extended evaluated BED15")
    parser.add_argument("--threshold", type=float, default=0.5,
                        help="Minimum DeepPASS score to retain a site (default: 0.5)")
    parser.add_argument("--out-tsv",   required=True, help="Output site_catalog TSV path")
    args = parser.parse_args()

    inputs = [
        (Path(args.exonic),   "scapture_exonic"),
        (Path(args.intronic), "scapture_intronic"),
        (Path(args.extended), "scapture_3prime_extended"),
    ]

    rows: list[dict] = []
    seen: set[tuple] = set()

    for bed_path, site_class in inputs:
        for row in parse_bed15(bed_path, site_class, args.threshold):
            key = (row["chrom"], row["start"], row["end"], row["strand"])
            if key in seen:
                continue
            seen.add(key)
            rows.append(row)

    rows.sort(key=lambda r: (r["chrom"], r["start"], r["strand"]))

    with Path(args.out_tsv).open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"[run_scapture_filter] wrote {len(rows)} sites to {args.out_tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
