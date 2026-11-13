#!/usr/bin/env python3
"""Aggregate staged SCAPTURE metrics into the monolith's manifest schema."""

import argparse
import csv
from pathlib import Path


def read_pairs(path):
    with Path(path).open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows or rows[0] not in (["metric", "value"], ["field", "value"]):
        raise ValueError(f"Invalid key/value header: {path}")
    if any(len(row) != 2 for row in rows[1:]):
        raise ValueError(f"Invalid key/value row: {path}")
    values = dict(rows[1:])
    if len(values) != len(rows) - 1:
        raise ValueError(f"Duplicate key: {path}")
    return values


def read_reference(path):
    with Path(path).open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows or rows[0] != ["category", "integrated_rows"]:
        raise ValueError("Invalid PAS reference metrics header")
    values = dict(rows[1:])
    categories = ("exonic", "intronic", "3primeExtended")
    return sum(int(values[category]) > 0 for category in categories)


def main():
    parser = argparse.ArgumentParser()
    for name in ("library_id", "bam", "gtf", "fasta", "fai", "prepared_gtf",
                 "gene_map", "duplicate_report", "pascall_metrics",
                 "reference_metrics", "catalog_metrics", "quant_manifest",
                 "species", "read_length", "peak_width", "normality_seed",
                 "threshold", "output"):
        parser.add_argument("--" + name.replace("_", "-"), required=True)
    args = parser.parse_args()
    pascall = read_pairs(args.pascall_metrics)
    catalog = read_pairs(args.catalog_metrics)
    quant = read_pairs(args.quant_manifest)
    if quant.get("status") != "complete":
        raise ValueError("PASquant manifest does not report completion")
    duplicate_count = max(0, sum(1 for _ in Path(args.duplicate_report).open()) - 1)
    fields = [
        ("library_id", args.library_id), ("bam", args.bam), ("gtf", args.gtf),
        ("scapture_gtf", args.prepared_gtf), ("gene_name_map", args.gene_map),
        ("genome_fasta", args.fasta), ("genome_fai", args.fai),
        ("rescue_pascall_dir", "none"), ("species", args.species),
        ("read_length", args.read_length), ("peak_width", args.peak_width),
        ("normality_seed", args.normality_seed),
        ("deeppass_threshold", args.threshold),
    ]
    categories = (("exonic", "exonic"), ("intronic", "intronic"),
                  ("3primeExtended", "extended"))
    fields.extend((f"{label}_valid", pascall[f"{category}_valid"])
                  for category, label in categories)
    fields.extend((
        ("valid_peaks", pascall["valid_peaks"]),
    ))
    fields.extend((f"{label}_scored", pascall[f"{category}_scored"])
                  for category, label in categories)
    fields.extend((
        ("scored_peaks", pascall["scored_peaks"]),
        ("site_catalog_rows", catalog["site_catalog_rows"]),
        ("merged_categories", str(read_reference(args.reference_metrics))),
        ("duplicate_original_gene_names", str(duplicate_count)),
        ("source", "SCAPTURE"),
        ("status", "complete"),
        ("quant_matrix", quant["quant_matrix"]),
        ("matrix_site_rows", quant["matrix_site_rows"]),
        ("matrix_cell_columns", quant["matrix_cell_columns"]),
    ))
    with Path(args.output).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("field", "value"))
        writer.writerows(fields)


if __name__ == "__main__":
    main()
