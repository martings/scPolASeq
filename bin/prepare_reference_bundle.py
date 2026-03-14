#!/usr/bin/env python3

import argparse
import csv
import json
import shutil
from pathlib import Path


def read_fasta_lengths(path: Path):
    lengths = []
    name = None
    seq_len = 0
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths.append((name, seq_len))
                name = line[1:].split()[0]
                seq_len = 0
            else:
                seq_len += len(line)
    if name is not None:
        lengths.append((name, seq_len))
    return lengths


def copy_or_touch(source: Path, target: Path):
    if source.exists() and source.name != "NO_FILE" and source.stat().st_size > 0:
        shutil.copyfile(source, target)
    else:
        target.write_text("", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome-fasta", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--known-polya", required=True)
    parser.add_argument("--priming-blacklist", required=True)
    parser.add_argument("--out-fasta", required=True)
    parser.add_argument("--out-gtf", required=True)
    parser.add_argument("--out-fai", required=True)
    parser.add_argument("--out-dict", required=True)
    parser.add_argument("--out-sizes", required=True)
    parser.add_argument("--out-known-polya", required=True)
    parser.add_argument("--out-blacklist", required=True)
    parser.add_argument("--out-manifest", required=True)
    args = parser.parse_args()

    genome_fasta = Path(args.genome_fasta)
    out_fasta = Path(args.out_fasta)
    out_gtf = Path(args.out_gtf)
    shutil.copyfile(genome_fasta, out_fasta)
    shutil.copyfile(Path(args.gtf), out_gtf)

    lengths = read_fasta_lengths(out_fasta)
    with Path(args.out_sizes).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerows(lengths)

    with Path(args.out_fai).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for contig, length in lengths:
            writer.writerow([contig, length, 0, 0, 0])

    with Path(args.out_dict).open("w", encoding="utf-8") as handle:
        handle.write("@HD\tVN:1.6\tSO:unsorted\n")
        for contig, length in lengths:
            handle.write(f"@SQ\tSN:{contig}\tLN:{length}\n")

    copy_or_touch(Path(args.known_polya), Path(args.out_known_polya))
    copy_or_touch(Path(args.priming_blacklist), Path(args.out_blacklist))

    manifest = {
        "genome_fasta": str(out_fasta.resolve()),
        "gtf": str(out_gtf.resolve()),
        "contigs": len(lengths),
    }
    Path(args.out_manifest).write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
