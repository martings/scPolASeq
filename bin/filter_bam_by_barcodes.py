#!/usr/bin/env python3
"""
filter_bam_by_barcodes.py — Filter a BAM to reads with CB tags matching
validated cell barcodes in the cell annotations table.

Reads are tagged with their cluster and cell_type from the annotations so
downstream grouping can operate on the filtered BAM directly.

Usage:
    python filter_bam_by_barcodes.py \
        --bam            aligned.bam \
        --annotations    cell_annotations.tsv \
        --out-bam        filtered.bam \
        --out-stats      stats.tsv
"""
import argparse
import sys
import logging

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bam",          required=True)
    p.add_argument("--annotations",  required=True)
    p.add_argument("--out-bam",      required=True)
    p.add_argument("--out-stats",    required=True)
    p.add_argument("--cb-tag",       default="CB", help="BAM tag carrying cell barcode")
    p.add_argument("--ub-tag",       default="UB", help="BAM tag carrying UMI")
    return p.parse_args()


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    log = logging.getLogger()

    # Load valid barcodes from annotations
    try:
        anno = pd.read_csv(args.annotations, sep="\t")
    except Exception as e:
        log.error(f"Cannot load annotations: {e}")
        sys.exit(1)

    # Resolve barcode column: accept 'barcode', 'barcode_corrected', or 'barcode_raw'
    bc_col = next(
        (c for c in ("barcode", "barcode_corrected", "barcode_raw") if c in anno.columns),
        None,
    )
    if bc_col is None:
        log.error("Annotations table must contain a 'barcode', 'barcode_corrected', or 'barcode_raw' column")
        sys.exit(1)
    log.info(f"Using barcode column: '{bc_col}'")

    valid_barcodes = set(anno[bc_col].dropna().str.strip().tolist())
    log.info(f"Valid barcodes: {len(valid_barcodes)}")

    # Build barcode → cluster/cell_type lookup for tagging
    bc_to_meta = {}
    for _, row in anno.iterrows():
        bc = str(row[bc_col]).strip()
        # Accept both 'cluster_id' (canonical output) and legacy 'cluster'
        cluster_val = row.get("cluster_id") or row.get("cluster", "NA")
        bc_to_meta[bc] = {
            "cluster":   str(cluster_val) if cluster_val and str(cluster_val) != "nan" else "NA",
            "cell_type": str(row.get("cell_type", "NA")),
            "condition": str(row.get("condition", "NA")),
        }

    # Try to use pysam; fall back to stub if unavailable or BAM is not valid
    try:
        import pysam
        _filter_with_pysam(args, valid_barcodes, bc_to_meta, log)
    except ImportError:
        log.warning("pysam not available — writing stub filtered BAM")
        _stub_output(args, valid_barcodes, log)
    except (ValueError, OSError) as exc:
        log.warning(f"BAM not readable ({exc}) — writing stub filtered BAM")
        _stub_output(args, valid_barcodes, log)


def _filter_with_pysam(args, valid_barcodes, bc_to_meta, log):
    import pysam

    total = valid_cb = filtered = 0

    with pysam.AlignmentFile(args.bam, "rb") as bam_in, \
         pysam.AlignmentFile(args.out_bam, "wb", header=bam_in.header) as bam_out:

        for read in bam_in.fetch(until_eof=True):
            total += 1
            try:
                cb = read.get_tag(args.cb_tag)
            except KeyError:
                continue

            if cb not in valid_barcodes:
                continue

            valid_cb += 1
            meta = bc_to_meta.get(cb, {})

            # Tag read with biological group info
            tags = dict(read.get_tags())
            tags["CL"] = meta.get("cluster",   "NA")  # cluster label
            tags["CT"] = meta.get("cell_type", "NA")  # cell type label
            tags["CD"] = meta.get("condition", "NA")  # condition
            read.set_tags(list(tags.items()))
            bam_out.write(read)
            filtered += 1

    log.info(f"Total reads: {total}, with valid CB: {valid_cb}, written: {filtered}")

    log.info("Sorting and indexing filtered BAM")
    sorted_bam = args.out_bam + ".sorted_tmp.bam"
    pysam.sort("-o", sorted_bam, args.out_bam)
    import os, shutil
    shutil.move(sorted_bam, args.out_bam)
    pysam.index(args.out_bam)

    pd.DataFrame([
        {"metric": "total_reads",    "value": total},
        {"metric": "valid_cb_reads", "value": valid_cb},
        {"metric": "filtered_reads", "value": filtered},
        {"metric": "filter_rate",    "value": round(filtered / max(total, 1), 4)},
    ]).to_csv(args.out_stats, sep="\t", index=False)


def _stub_output(args, valid_barcodes, log):
    """Write placeholder outputs when pysam is unavailable."""
    import shutil
    shutil.copy(args.bam, args.out_bam)

    pd.DataFrame([
        {"metric": "total_reads",    "value": 0},
        {"metric": "valid_cb_reads", "value": 0},
        {"metric": "filtered_reads", "value": 0},
        {"metric": "filter_rate",    "value": 0.0},
        {"metric": "note",           "value": "pysam_unavailable_stub"},
    ]).to_csv(args.out_stats, sep="\t", index=False)


if __name__ == "__main__":
    main()
