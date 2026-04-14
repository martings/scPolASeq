#!/usr/bin/env python3
"""
bam_to_bedgraph.py — Strand-specific coverage in bedGraph format using pysam.

Replaces bedtools genomecov for the COVERAGE_TRACKS module.
Outputs 0-based half-open intervals [start, end) with integer coverage.

Usage:
    python bam_to_bedgraph.py input.bam --strand + --output out.fwd.bedGraph
    python bam_to_bedgraph.py input.bam --strand - --output out.rev.bedGraph
"""
import argparse
import sys
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger()


def bam_to_bedgraph(bam_path: str, strand: str, output_path: str, min_mapq: int = 0):
    try:
        import pysam
    except ImportError:
        log.error("pysam is required but not installed.")
        sys.exit(1)

    want_reverse = strand == "-"
    # chrom_order preserves BAM header order for sorting
    chrom_order = {}
    chrom_len = {}

    log.info(f"Reading {bam_path} (strand={strand}, min_mapq={min_mapq})")

    # First pass: collect per-position coverage using pileup
    # pysam.pileup streams through the file column by column (memory efficient)
    coverage = {}  # chrom -> {pos: depth}

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for i, ref in enumerate(bam.header.references):
            chrom_order[ref] = i
            chrom_len[ref] = bam.header.lengths[i]

        for col in bam.pileup(stepper="nofilter", min_base_quality=0,
                               ignore_overlaps=False, ignore_orphans=False):
            chrom = col.reference_name
            pos   = col.reference_pos

            depth = 0
            for pileup_read in col.pileups:
                if pileup_read.is_del or pileup_read.is_refskip:
                    continue
                read = pileup_read.alignment
                if read.is_unmapped:
                    continue
                if read.mapping_quality < min_mapq:
                    continue
                if read.is_reverse != want_reverse:
                    continue
                depth += 1

            if depth > 0:
                if chrom not in coverage:
                    coverage[chrom] = {}
                coverage[chrom][pos] = depth

    n_total = sum(len(v) for v in coverage.values())
    log.info(f"Covered positions: {n_total:,} across {len(coverage)} chromosomes")

    if n_total == 0:
        log.warning("No coverage found — writing empty bedGraph")
        Path(output_path).touch()
        return

    # Write bedGraph: merge adjacent positions with identical coverage
    sorted_chroms = sorted(coverage.keys(), key=lambda c: chrom_order.get(c, 999999))
    written = 0
    with open(output_path, "w") as out:
        for chrom in sorted_chroms:
            positions = sorted(coverage[chrom])
            if not positions:
                continue

            start = positions[0]
            end   = start + 1
            val   = coverage[chrom][start]

            for pos in positions[1:]:
                depth = coverage[chrom][pos]
                if pos == end and depth == val:
                    end += 1
                else:
                    out.write(f"{chrom}\t{start}\t{end}\t{val}\n")
                    written += 1
                    start = pos
                    end   = pos + 1
                    val   = depth

            out.write(f"{chrom}\t{start}\t{end}\t{val}\n")
            written += 1

    log.info(f"Written {written:,} bedGraph intervals to {output_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bam",           help="Input BAM (must be sorted+indexed)")
    parser.add_argument("--strand",      required=True, choices=["+", "-"],
                        help="Strand to compute coverage for")
    parser.add_argument("--output",      required=True, help="Output bedGraph path")
    parser.add_argument("--min-mapq",    type=int, default=0,
                        help="Minimum mapping quality (default: 0)")
    args = parser.parse_args()

    if not Path(args.bam).exists():
        log.error(f"BAM file not found: {args.bam}")
        sys.exit(1)

    bam_to_bedgraph(args.bam, args.strand, args.output, args.min_mapq)


if __name__ == "__main__":
    main()
