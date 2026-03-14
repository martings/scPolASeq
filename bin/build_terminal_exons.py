#!/usr/bin/env python3

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


ATTRIBUTE_RE = re.compile(r'(\S+)\s+"([^"]+)"')


def parse_attributes(field: str):
    return {key: value for key, value in ATTRIBUTE_RE.findall(field)}


def load_exons(gtf_path: Path):
    transcripts = defaultdict(list)
    with gtf_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            chrom, source, feature, start, end, score, strand, frame, attributes = fields
            attrs = parse_attributes(attributes)
            gene_id = attrs.get("gene_id", "")
            transcript_id = attrs.get("transcript_id", "")
            exon_number = attrs.get("exon_number", "")
            transcripts[transcript_id].append(
                {
                    "gene_id": gene_id,
                    "transcript_id": transcript_id,
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "exon_number": exon_number,
                }
            )
    return transcripts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--out-bed", required=True)
    args = parser.parse_args()

    transcripts = load_exons(Path(args.gtf))
    terminal_rows = []
    for transcript_id, exons in transcripts.items():
        strand = exons[0]["strand"]
        if strand == "+":
            terminal = sorted(exons, key=lambda exon: (exon["end"], exon["start"]))[-1]
        else:
            terminal = sorted(exons, key=lambda exon: (exon["start"], exon["end"]))[0]
        terminal_rows.append(
            {
                "gene_id": terminal["gene_id"],
                "transcript_id": transcript_id,
                "chrom": terminal["chrom"],
                "start": terminal["start"],
                "end": terminal["end"],
                "strand": terminal["strand"],
                "exon_number": terminal["exon_number"],
                "terminal_exon_rank": 1,
            }
        )

    terminal_rows.sort(key=lambda row: (row["gene_id"], row["transcript_id"]))

    with Path(args.out_tsv).open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "gene_id",
            "transcript_id",
            "chrom",
            "start",
            "end",
            "strand",
            "exon_number",
            "terminal_exon_rank",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(terminal_rows)

    with Path(args.out_bed).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for row in terminal_rows:
            writer.writerow(
                [
                    row["chrom"],
                    row["start"] - 1,
                    row["end"],
                    f"{row['gene_id']}|{row['transcript_id']}",
                    0,
                    row["strand"],
                ]
            )


if __name__ == "__main__":
    main()
