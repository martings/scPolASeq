#!/usr/bin/env python3

"""
Prepare a SCAPTURE-specific GTF with collision-free gene tokens.

SCAPTURE uses the `gene_name` token from the annotation-derived intermediate
files as both a logical gene identifier and a per-gene output filename stem.
Many real annotations contain non-unique gene_name values (for example Y_RNA,
Metazoa_SRP, U6, SNORA70 families), which causes concurrent PAScall jobs to
overwrite one another.

To make SCAPTURE safe without perturbing the rest of the pipeline, this helper
rewrites `gene_name` to the stable `gene_id` for every feature that has a
gene_id. The original gene_name is retained in a sidecar mapping table and in
an extra GTF attribute for traceability.
"""

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Set, Tuple


def parse_attributes(raw: str) -> List[Tuple[str, str]]:
    attrs: List[Tuple[str, str]] = []
    for field in raw.strip().split(";"):
        field = field.strip()
        if not field:
            continue
        if " " not in field:
            attrs.append((field, ""))
            continue
        key, value = field.split(" ", 1)
        attrs.append((key, value.strip().strip('"')))
    return attrs


def format_attributes(attrs: List[Tuple[str, str]]) -> str:
    parts = []
    for key, value in attrs:
        if value == "":
            parts.append(key)
        else:
            parts.append(f'{key} "{value}"')
    return "; ".join(parts) + ";"


def replace_gene_name(attrs: List[Tuple[str, str]]) -> Tuple[List[Tuple[str, str]], str, str]:
    gene_id = ""
    original_gene_name = ""
    rewritten: List[Tuple[str, str]] = []
    saw_gene_name = False
    saw_original_tag = False

    for key, value in attrs:
        if key == "gene_id" and value:
            gene_id = value
        if key == "gene_name" and value:
            original_gene_name = value
        if key == "scapture_original_gene_name":
            saw_original_tag = True

    if not gene_id:
        return attrs, "", original_gene_name

    scapture_gene_name = gene_id
    original_gene_name = original_gene_name or gene_id

    for key, value in attrs:
        if key == "gene_name":
            rewritten.append((key, scapture_gene_name))
            saw_gene_name = True
        else:
            rewritten.append((key, value))

    if not saw_gene_name:
        rewritten.append(("gene_name", scapture_gene_name))

    if not saw_original_tag:
        rewritten.append(("scapture_original_gene_name", original_gene_name))

    return rewritten, gene_id, original_gene_name


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", required=True, help="Input annotation GTF")
    parser.add_argument("--out-gtf", required=True, help="Rewritten GTF for SCAPTURE")
    parser.add_argument("--out-map", required=True, help="Gene token mapping TSV")
    parser.add_argument("--out-duplicate-report", required=True, help="Duplicate original gene_name report TSV")
    args = parser.parse_args()

    in_gtf = Path(args.gtf)
    out_gtf = Path(args.out_gtf)
    out_map = Path(args.out_map)
    out_dup = Path(args.out_duplicate_report)

    gene_rows: Dict[str, Tuple[str, str]] = {}

    with in_gtf.open("r", encoding="utf-8") as src, out_gtf.open("w", encoding="utf-8") as dst:
        for line in src:
            if not line or line.startswith("#"):
                dst.write(line)
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                dst.write(line)
                continue

            attrs = parse_attributes(cols[8])
            rewritten, gene_id, original_gene_name = replace_gene_name(attrs)
            if gene_id:
                gene_rows.setdefault(gene_id, (gene_id, original_gene_name or gene_id))
                cols[8] = format_attributes(rewritten)
                dst.write("\t".join(cols) + "\n")
            else:
                dst.write(line)

    with out_map.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene_id", "scapture_gene_name", "original_gene_name"])
        for gene_id, original_gene_name in sorted(gene_rows.values()):
            writer.writerow([gene_id, gene_id, original_gene_name])

    counts: Dict[str, Set[str]] = {}
    for gene_id, original_gene_name in gene_rows.values():
        counts.setdefault(original_gene_name, set()).add(gene_id)

    duplicate_rows = sorted(
        ((len(gene_ids), original_name) for original_name, gene_ids in counts.items() if len(gene_ids) > 1),
        key=lambda row: (-row[0], row[1]),
    )

    with out_dup.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene_count", "original_gene_name"])
        writer.writerows(duplicate_rows)


if __name__ == "__main__":
    main()
