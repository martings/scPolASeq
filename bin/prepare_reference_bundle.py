#!/usr/bin/env python3

import argparse
import csv
import json
import shutil
from pathlib import Path


def read_fasta_index(path: Path):
    """Return faidx-compatible rows for a FASTA file."""
    entries = []
    name = None
    seq_len = 0
    seq_offset = 0
    line_bases = 0
    line_width = 0

    def flush_entry():
        if name is not None:
            entries.append((name, seq_len, seq_offset, line_bases, line_width))

    with path.open("rb") as handle:
        for raw_line in handle:
            line = raw_line.rstrip(b"\r\n")
            if not line:
                continue
            if line.startswith(b">"):
                flush_entry()
                name = line[1:].split()[0].decode("utf-8")
                seq_len = 0
                seq_offset = handle.tell()
                line_bases = 0
                line_width = 0
            else:
                if line_bases == 0:
                    line_bases = len(line)
                    line_width = len(raw_line)
                seq_len += len(line)
    flush_entry()
    return entries


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

    fasta_index = read_fasta_index(out_fasta)
    with Path(args.out_sizes).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerows((contig, length) for contig, length, *_ in fasta_index)

    with Path(args.out_fai).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerows(fasta_index)

    with Path(args.out_dict).open("w", encoding="utf-8") as handle:
        handle.write("@HD\tVN:1.6\tSO:unsorted\n")
        for contig, length, *_ in fasta_index:
            handle.write(f"@SQ\tSN:{contig}\tLN:{length}\n")

    copy_or_touch(Path(args.known_polya), Path(args.out_known_polya))
    copy_or_touch(Path(args.priming_blacklist), Path(args.out_blacklist))

    manifest = {
        "genome_fasta": str(out_fasta.resolve()),
        "gtf": str(out_gtf.resolve()),
        "genome_fai": str(Path(args.out_fai).resolve()),
        "contigs": len(fasta_index),
    }
    Path(args.out_manifest).write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
