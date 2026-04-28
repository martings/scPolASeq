#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def load_terminal_exons(path: Path):
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_known_polya(path: Path):
    if not path.exists() or path.name == "NO_FILE" or path.stat().st_size == 0:
        return []
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--terminal-exons", required=True)
    parser.add_argument("--known-polya", required=True)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--out-bed", required=True)
    args = parser.parse_args()

    terminal_exons = load_terminal_exons(Path(args.terminal_exons))
    known_polya = load_known_polya(Path(args.known_polya))

    site_rows = []
    seen = {}

    for exon in terminal_exons:
        strand = exon["strand"]
        site_position = int(exon["end"]) if strand == "+" else int(exon["start"])
        key = (exon["gene_id"], exon["chrom"], site_position, strand)
        row = {
            "site_id": f"{exon['gene_id']}:{exon['chrom']}:{site_position}:{strand}",
            "gene_id": exon["gene_id"],
            "chrom": exon["chrom"],
            "start": site_position,
            "end": site_position,
            "strand": strand,
            "site_class": "annotated_terminal_end",
            "site_source": "annotation",
            "priming_flag": "unchecked",
        }
        site_rows.append(row)
        seen[key] = row

    for atlas in known_polya:
        # gene_id: support PolyASite (gene_id), PolyA_DB (Gene/Alias), gene_name fallback
        gene_id = (
            atlas.get("gene_id")
            or atlas.get("Gene")
            or atlas.get("gene")
            or atlas.get("gene_name")
            or "unknown_gene"
        )
        # chrom: BED 'chrom', PolyA_DB 'Chr'/'chr'
        chrom = atlas.get("chrom") or atlas.get("Chr") or atlas.get("chr") or ""
        # coordinates: canonical → PolyASite BED chromStart/chromEnd → PolyA_DB Position
        start = int(
            atlas.get("start")
            or atlas.get("chromStart")
            or atlas.get("Position")
            or atlas.get("position")
            or atlas.get("end")
            or atlas.get("chromEnd")
            or 0
        )
        end = int(
            atlas.get("end")
            or atlas.get("chromEnd")
            or atlas.get("Position")
            or atlas.get("position")
            or start
        )
        strand = atlas.get("strand") or atlas.get("Strand") or "+"
        atlas_source = atlas.get("source") or atlas.get("reference_source") or "atlas"
        key = (gene_id, chrom, end if strand == "+" else start, strand)
        if key in seen:
            seen[key]["site_source"] = f"annotation+atlas:{atlas_source}"
            continue
        row = {
            "site_id": f"{gene_id}:{chrom}:{start}-{end}:{strand}",
            "gene_id": gene_id,
            "chrom": chrom,
            "start": start,
            "end": end,
            "strand": strand,
            "site_class": "known_pas",
            "site_source": f"atlas:{atlas_source}",
            "priming_flag": "unchecked",
        }
        site_rows.append(row)
        seen[key] = row

    site_rows.sort(key=lambda row: (row["gene_id"], row["chrom"], row["start"], row["strand"]))

    with Path(args.out_tsv).open("w", newline="", encoding="utf-8") as handle:
        fieldnames = ["site_id", "gene_id", "chrom", "start", "end", "strand", "site_class", "site_source", "priming_flag"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(site_rows)

    with Path(args.out_bed).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for row in site_rows:
            writer.writerow([row["chrom"], int(row["start"]) - 1, row["end"], row["site_id"], 0, row["strand"]])


if __name__ == "__main__":
    main()
