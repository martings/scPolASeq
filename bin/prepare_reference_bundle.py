#!/usr/bin/env python3
"""Stage, validate, and hash reference inputs and STAR index metadata."""

import argparse
import csv
import hashlib
import json
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# Binary genome index format written and consumed by STAR 2.7.4a through 2.7.11b.
SUPPORTED_STAR_GENOME_VERSIONS = {"2.7.4a"}
REQUIRED_STAR_INDEX_FILES = ("Genome", "SA", "SAindex", "genomeParameters.txt")


def is_sentinel(path: Path) -> bool:
    return path.name.startswith("NO_FILE")


def sha256_file(path: Path) -> str:
    with path.open("rb") as handle:
        if hasattr(hashlib, "file_digest"):
            return hashlib.file_digest(handle, "sha256").hexdigest()
        digest = hashlib.sha256()
        for chunk in iter(lambda: handle.read(1048576), b""):
            digest.update(chunk)
        return digest.hexdigest()


def require_non_empty_file(path: Path, label: str) -> None:
    if not path.exists() or not path.is_file() or path.stat().st_size <= 0:
        raise SystemExit(f"ERROR: Required {label} file is missing or empty: {path}")


def read_fasta_index(path: Path) -> List[Tuple[str, int, int, int, int]]:
    """Return faidx-compatible rows for a FASTA file."""
    entries = []
    name = None
    seq_len = 0
    seq_offset = 0
    line_bases = 0
    line_width = 0

    def flush_entry() -> None:
        if name is not None:
            entries.append((name, seq_len, seq_offset, line_bases, line_width))

    with path.open("rb") as handle:
        for raw_line in handle:
            line = raw_line.rstrip(b"\r\n")
            if not line:
                continue
            if line.startswith(b">"):
                flush_entry()
                header = line[1:].split()
                if not header:
                    raise SystemExit(f"ERROR: Empty FASTA header in {path}")
                name = header[0].decode("utf-8")
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
    if not entries:
        raise SystemExit(f"ERROR: FASTA file contains no sequences: {path}")
    return entries


def parse_chrom_sizes(path: Path, label: str) -> Dict[str, int]:
    require_non_empty_file(path, label)
    sizes: Dict[str, int] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, start=1):
            stripped = raw.strip()
            if not stripped:
                continue
            parts = stripped.split("\t")
            if len(parts) < 2:
                raise SystemExit(
                    f"ERROR: Malformed row {line_no} in {label} ({path}): expected '<contig>\\t<length>'"
                )
            contig = parts[0].strip()
            try:
                length = int(parts[1].strip())
            except ValueError as exc:
                raise SystemExit(
                    f"ERROR: Invalid contig length on row {line_no} in {label} ({path}): {parts[1]!r}"
                ) from exc
            if not contig or length <= 0:
                raise SystemExit(
                    f"ERROR: Non-positive or empty contig entry on row {line_no} in {label} ({path})"
                )
            sizes[contig] = length
    if not sizes:
        raise SystemExit(f"ERROR: {label} contains no contig entries: {path}")
    return sizes


def parse_genome_parameters(path: Path) -> Dict[str, str]:
    require_non_empty_file(path, "STAR genomeParameters.txt")
    params: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            if "\t" in line:
                key, val = line.split("\t", 1)
                params[key.strip()] = val.strip()
    return params


def validate_star_index(
    index_dir: Path,
    fasta_contigs: Optional[Dict[str, int]] = None,
    star_binary_version: Optional[str] = None,
) -> Dict[str, str]:
    if not index_dir.exists() or not index_dir.is_dir():
        raise SystemExit(f"ERROR: STAR index directory does not exist: {index_dir}")

    for filename in REQUIRED_STAR_INDEX_FILES:
        require_non_empty_file(index_dir / filename, f"STAR index ({filename})")

    params = parse_genome_parameters(index_dir / "genomeParameters.txt")
    version_genome = params.get("versionGenome")
    if not version_genome:
        raise SystemExit(
            f"ERROR: STAR index genomeParameters.txt missing 'versionGenome': {index_dir / 'genomeParameters.txt'}"
        )
    if version_genome not in SUPPORTED_STAR_GENOME_VERSIONS:
        supported = ", ".join(sorted(SUPPORTED_STAR_GENOME_VERSIONS))
        star_label = star_binary_version or "STAR 2.7.11b"
        raise SystemExit(
            f"ERROR: Incompatible STAR index versionGenome={version_genome!r} in {index_dir}; "
            f"{star_label} requires format [{supported}]. Rebuild without --star_index."
        )

    chr_name_length = index_dir / "chrNameLength.txt"
    if chr_name_length.exists() and fasta_contigs is not None:
        index_contigs = parse_chrom_sizes(chr_name_length, "STAR index chrNameLength.txt")
        if index_contigs != fasta_contigs:
            raise SystemExit(
                f"ERROR: STAR index contigs/lengths ({index_dir / 'chrNameLength.txt'}) "
                f"do not match genome_fasta ({len(index_contigs)} vs {len(fasta_contigs)} contigs)."
            )

    return params


def summarize_star_index(
    index_dir: Path,
    mode: str,
    fasta_contigs: Optional[Dict[str, int]],
    star_binary_version: Optional[str],
    source_path: Optional[str],
) -> dict:
    params = validate_star_index(index_dir, fasta_contigs=fasta_contigs, star_binary_version=star_binary_version)
    file_digests: Dict[str, dict] = {}
    agg = hashlib.sha256()
    for name in sorted(set(REQUIRED_STAR_INDEX_FILES) | {"chrNameLength.txt"}):
        fpath = index_dir / name
        if fpath.exists() and fpath.is_file():
            digest = sha256_file(fpath)
            size = fpath.stat().st_size
            file_digests[name] = {"bytes": size, "sha256": digest}
            agg.update(f"{name}:{size}:{digest}\n".encode("utf-8"))

    return {
        "mode": mode,
        "source": source_path,
        "compatible": True,
        "star_binary_version": star_binary_version or "unknown",
        "version_genome": params.get("versionGenome"),
        "expected_star_version_genome": sorted(SUPPORTED_STAR_GENOME_VERSIONS),
        "sjdb_overhang": params.get("sjdbOverhang"),
        "genome_file_sizes": params.get("genomeFileSizes"),
        "parameters_sha256": file_digests["genomeParameters.txt"]["sha256"],
        "index_sha256": agg.hexdigest(),
        "files": file_digests,
    }


def copy_optional_input(source: Path, target: Path, label: str) -> dict:
    if is_sentinel(source):
        target.write_text("", encoding="utf-8")
        return {"present": False, "source": None, "path": str(target.resolve()), "bytes": 0, "sha256": None}
    require_non_empty_file(source, label)
    shutil.copyfile(source, target)
    return {
        "present": True,
        "source": str(source),
        "path": str(target.resolve()),
        "bytes": target.stat().st_size,
        "sha256": sha256_file(target),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome-fasta")
    parser.add_argument("--gtf")
    parser.add_argument("--known-polya")
    parser.add_argument("--priming-blacklist")
    parser.add_argument("--chrom-sizes", default=None)
    parser.add_argument("--prebuilt-star-index", default=None)
    parser.add_argument("--reference-root", default=None)
    parser.add_argument("--reference-build", default=None)
    parser.add_argument("--star-binary-version", default=None)
    parser.add_argument("--out-fasta")
    parser.add_argument("--out-gtf")
    parser.add_argument("--out-fai")
    parser.add_argument("--out-dict")
    parser.add_argument("--out-sizes")
    parser.add_argument("--out-known-polya")
    parser.add_argument("--out-blacklist")
    parser.add_argument("--out-manifest", required=True)
    parser.add_argument("--finalize-star-index", default=None)
    parser.add_argument("--star-index-mode", choices=["prebuilt", "generated"], default="generated")
    parser.add_argument("--star-index-source", default=None)
    args = parser.parse_args()

    manifest_path = Path(args.out_manifest)

    if args.finalize_star_index:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        fasta_contigs = manifest.get("_fasta_contigs")
        if fasta_contigs is not None:
            fasta_contigs = {k: int(v) for k, v in fasta_contigs.items()}
        manifest.pop("_fasta_contigs", None)
        star_ver = args.star_binary_version or manifest.get("star_binary_version")
        manifest["star_binary_version"] = star_ver or "unknown"
        manifest["star_index"] = summarize_star_index(
            Path(args.finalize_star_index),
            mode=args.star_index_mode,
            fasta_contigs=fasta_contigs,
            star_binary_version=star_ver,
            source_path=args.star_index_source,
        )
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        return

    for req in ("genome_fasta", "gtf", "known_polya", "priming_blacklist",
                "out_fasta", "out_gtf", "out_fai", "out_dict", "out_sizes",
                "out_known_polya", "out_blacklist"):
        if not getattr(args, req):
            raise SystemExit(f"ERROR: Missing required argument --{req.replace('_', '-')}")

    genome_fasta = Path(args.genome_fasta)
    gtf = Path(args.gtf)
    require_non_empty_file(genome_fasta, "genome_fasta")
    require_non_empty_file(gtf, "gtf")

    out_fasta = Path(args.out_fasta)
    out_gtf = Path(args.out_gtf)
    shutil.copyfile(genome_fasta, out_fasta)
    shutil.copyfile(gtf, out_gtf)

    fasta_index = read_fasta_index(out_fasta)
    fasta_contigs = {contig: length for contig, length, *_ in fasta_index}

    chrom_sizes_source = None
    if args.chrom_sizes:
        cs_path = Path(args.chrom_sizes)
        if not is_sentinel(cs_path):
            supplied_sizes = parse_chrom_sizes(cs_path, "chrom_sizes")
            if supplied_sizes != fasta_contigs:
                raise SystemExit(
                    f"ERROR: Supplied chrom_sizes ({cs_path}) contigs/lengths do not match genome_fasta "
                    f"({len(supplied_sizes)} vs {len(fasta_contigs)} contigs)."
                )
            chrom_sizes_source = str(cs_path)

    out_sizes = Path(args.out_sizes)
    with out_sizes.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows((contig, length) for contig, length, *_ in fasta_index)

    out_fai = Path(args.out_fai)
    with out_fai.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(fasta_index)

    with Path(args.out_dict).open("w", encoding="utf-8") as handle:
        handle.write("@HD\tVN:1.6\tSO:unsorted\n")
        for contig, length, *_ in fasta_index:
            handle.write(f"@SQ\tSN:{contig}\tLN:{length}\n")

    polya_meta = copy_optional_input(Path(args.known_polya), Path(args.out_known_polya), "known_polya")
    blacklist_meta = copy_optional_input(
        Path(args.priming_blacklist), Path(args.out_blacklist), "priming_blacklist"
    )

    if args.prebuilt_star_index:
        prebuilt_path = Path(args.prebuilt_star_index)
        if not is_sentinel(prebuilt_path):
            validate_star_index(
                prebuilt_path,
                fasta_contigs=fasta_contigs,
                star_binary_version=args.star_binary_version,
            )

    manifest = {
        "manifest_version": "1.1",
        "reference_root": args.reference_root or None,
        "reference_build": args.reference_build or None,
        "star_binary_version": args.star_binary_version or "pending",
        "genome_fasta": str(out_fasta.resolve()),
        "gtf": str(out_gtf.resolve()),
        "genome_fai": str(out_fai.resolve()),
        "contigs": len(fasta_index),
        "inputs": {
            "genome_fasta": {
                "source": str(genome_fasta),
                "path": str(out_fasta.resolve()),
                "bytes": out_fasta.stat().st_size,
                "sha256": sha256_file(out_fasta),
                "contigs": len(fasta_index),
            },
            "gtf": {
                "source": str(gtf),
                "path": str(out_gtf.resolve()),
                "bytes": out_gtf.stat().st_size,
                "sha256": sha256_file(out_gtf),
            },
            "chrom_sizes": {
                "source": chrom_sizes_source or "generated_from_fasta",
                "path": str(out_sizes.resolve()),
                "bytes": out_sizes.stat().st_size,
                "sha256": sha256_file(out_sizes),
            },
            "known_polya": polya_meta,
            "priming_blacklist": blacklist_meta,
        },
        "_fasta_contigs": fasta_contigs,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
