#!/usr/bin/env python3
"""
Prepare the effective polyA reference table for the pipeline.

The script supports three modes:

1. legacy: no PolyA_DB/PolyASite sources are requested, so --known-polya is
   normalized/copied to the internal TSV contract.
2. single source: only --polya-db or --polyasite is requested, so that source is
   normalized and written directly without collapsing nearby sites.
3. dual source: both sources are requested, so normalized rows are merged within
   --merge-distance bp on chrom/strand while preserving source provenance.

Output columns:
    gene_id  gene_name  chrom  start  end  strand  score  source
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import logging
import sys
import urllib.request
import zipfile
from pathlib import Path
from urllib.parse import urlparse


POLYADB_URL = "https://exon.apps.wistar.org/PolyA_DB/v3/download/3.2/human_pas.zip"
POLYASITE_URL = (
    "https://polyasite.unibas.ch/download/atlas/3.0/GRCh38.GENCODE_42/"
    "atlas.clusters.3.0.GRCh38.GENCODE_42.bed.gz"
)

POLYADB_CACHE_NAME = "human_pas.polyadb3.2.zip"
POLYASITE_CACHE_NAME = "atlas.clusters.3.0.GRCh38.GENCODE_42.bed.gz"

FIELDNAMES_OUT = [
    "gene_id",
    "gene_name",
    "chrom",
    "start",
    "end",
    "strand",
    "score",
    "source",
]

BED_COLS = [
    "chrom",
    "chromStart",
    "chromEnd",
    "name",
    "score",
    "strand",
    "rep",
    "frac_samples",
    "nr_prots",
    "annotation",
    "gene_name",
    "gene_id",
    "repSite_signals",
]

SKIP_VALUES = {"", "skip", "none", "null", "false", "no", "0", "no_file"}
DOWNLOAD_VALUES = {"download", "auto", "true", "yes", "1"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--known-polya", required=True)
    parser.add_argument("--polya-db", default="skip")
    parser.add_argument("--polya-db-file", default="")
    parser.add_argument("--polyasite", default="skip")
    parser.add_argument("--polyasite-file", default="")
    parser.add_argument("--polyadb-url", default=POLYADB_URL)
    parser.add_argument("--polyasite-url", default=POLYASITE_URL)
    parser.add_argument("--cache-dir", required=True)
    parser.add_argument("--merge-distance", type=int, default=25)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--out-manifest", required=True)
    return parser.parse_args()


def clean(value, default: str = "") -> str:
    text = str(value).strip() if value is not None else ""
    return default if text.lower() in {"", "na", "nan", "none", "null"} else text


def normalize_chrom(value) -> str:
    chrom = clean(value)
    if not chrom:
        return chrom
    if chrom in {"M", "MT"}:
        return "chrM"
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def safe_int(value, default: int = 0) -> int:
    try:
        return int(float(clean(value, str(default))))
    except (TypeError, ValueError):
        return default


def normalize_strand(value) -> str:
    strand = clean(value, "+")
    return strand if strand in {"+", "-", "."} else "+"


def is_missing_file(path: Path) -> bool:
    return (not path.exists()) or path.name.startswith("NO_FILE") or path.stat().st_size == 0


def is_url(spec: str) -> bool:
    parsed = urlparse(spec)
    return parsed.scheme in {"http", "https", "ftp"}


def cache_filename(url: str, default_name: str) -> str:
    name = Path(urlparse(url).path).name
    return name or default_name


def download_if_missing(url: str, cache_dir: Path, filename: str, label: str) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    target = cache_dir / filename
    if target.exists() and target.stat().st_size > 0:
        logging.info("%s: using cached download %s", label, target)
        return target

    logging.info("%s: downloading %s", label, url)
    tmp = target.with_suffix(target.suffix + ".tmp")
    with urllib.request.urlopen(url, timeout=120) as response, tmp.open("wb") as handle:
        handle.write(response.read())
    tmp.replace(target)
    logging.info("%s: cached %s (%d bytes)", label, target, target.stat().st_size)
    return target


def resolve_source(
    *,
    spec: str,
    staged_file: str,
    default_url: str,
    default_cache_name: str,
    cache_dir: Path,
    label: str,
) -> Path | None:
    if staged_file:
        staged_path = Path(staged_file)
        if not is_missing_file(staged_path):
            logging.info("%s: using staged file %s", label, staged_path)
            return staged_path

    spec = clean(spec)
    spec_lower = spec.lower()
    if spec_lower in SKIP_VALUES:
        logging.info("%s: skipped", label)
        return None

    if spec_lower in DOWNLOAD_VALUES:
        return download_if_missing(default_url, cache_dir, default_cache_name, label)

    if is_url(spec):
        return download_if_missing(
            spec,
            cache_dir,
            cache_filename(spec, default_cache_name),
            label,
        )

    path = Path(spec)
    if is_missing_file(path):
        raise FileNotFoundError(
            f"{label}: file not found: {path}. Use 'download' to fetch the default URL."
        )
    logging.info("%s: using local file %s", label, path)
    return path


def read_text(path: Path) -> str:
    data = path.read_bytes()
    if data[:2] == b"\x1f\x8b":
        return gzip.decompress(data).decode("utf-8", errors="replace")
    return data.decode("utf-8", errors="replace")


def iter_zip_text_members(path: Path):
    with zipfile.ZipFile(path) as archive:
        candidates = [
            name
            for name in archive.namelist()
            if not name.endswith("/")
            and name.lower().endswith((".tsv", ".txt", ".bed"))
            and "readme" not in name.lower()
        ]
        if not candidates:
            raise ValueError(f"{path} contains no TSV/TXT/BED reference file")
        selected = sorted(candidates)[0]
        logging.info("PolyA_DB: reading %s from %s", selected, path)
        with archive.open(selected) as handle:
            return handle.read().decode("utf-8", errors="replace")


def read_rows_with_optional_header(text: str, bed_cols: list[str] | None = None):
    lines = [line for line in text.splitlines() if line.strip() and not line.startswith("track")]
    if not lines:
        return []

    first = lines[0].lstrip("#")
    if "\t" in first and any(token in first.split("\t") for token in ("chrom", "Chr", "PAS_ID")):
        return list(csv.DictReader(io.StringIO("\n".join([first] + lines[1:])), delimiter="\t"))

    if bed_cols is None:
        return list(csv.DictReader(io.StringIO("\n".join(lines)), delimiter="\t"))

    return [dict(zip(bed_cols, line.split("\t"))) for line in lines if not line.startswith("#")]


def normalize_known_polya(path: Path, source_label: str = "known_polya") -> list[dict[str, str]]:
    if is_missing_file(path):
        logging.info("known_polya: skipped empty/sentinel input")
        return []

    rows = read_rows_with_optional_header(read_text(path), BED_COLS)
    normalized = []
    for row in rows:
        chrom = normalize_chrom(row.get("chrom") or row.get("Chr") or row.get("chr"))
        if not chrom:
            continue
        start = safe_int(
            row.get("start")
            or row.get("chromStart")
            or row.get("Position")
            or row.get("position")
            or row.get("end")
            or row.get("chromEnd")
        )
        end = safe_int(
            row.get("end")
            or row.get("chromEnd")
            or row.get("Position")
            or row.get("position")
            or row.get("start")
            or row.get("chromStart"),
            start,
        )
        normalized.append(
            {
                "gene_id": clean(
                    row.get("gene_id")
                    or row.get("Gene")
                    or row.get("gene")
                    or row.get("gene_name"),
                    "unknown_gene",
                ),
                "gene_name": clean(row.get("gene_name") or row.get("Alias") or row.get("GeneSymbol")),
                "chrom": chrom,
                "start": str(start),
                "end": str(end),
                "strand": normalize_strand(row.get("strand") or row.get("Strand")),
                "score": clean(row.get("score") or row.get("Score"), "0"),
                "source": clean(row.get("source") or row.get("reference_source"), source_label),
            }
        )
    logging.info("known_polya: normalized %d rows", len(normalized))
    return normalized


def normalize_polyadb(path: Path) -> list[dict[str, str]]:
    text = iter_zip_text_members(path) if zipfile.is_zipfile(path) else read_text(path)
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    cols = reader.fieldnames or []

    # v3.2: PAS_ID="chr:pos:strand", direct cols Chromosome/Position/Strand, no GeneSymbol
    # v4.x: PAS_ID="chr:strand:pos", col GeneSymbol, no Chromosome
    is_v32 = "PAS_ID" in cols and "Chromosome" in cols
    is_v4  = "PAS_ID" in cols and "GeneSymbol" in cols and not is_v32

    logging.info("PolyA_DB: detected format %s", "v3.2" if is_v32 else "v4" if is_v4 else "legacy")
    rows_out = []

    for row in reader:
        if is_v32:
            # Use direct columns — more reliable than parsing PAS_ID
            chrom  = normalize_chrom(row.get("Chromosome") or row.get("Chr"))
            pos    = safe_int(row.get("Position") or row.get("position"))
            strand = normalize_strand(row.get("Strand") or row.get("strand"))
            # Prefer Ensembl ID; fall back to symbol; "na" rows become unknown_gene
            gene_id = clean(row.get("Ensemble ID") or row.get("Ensembl ID"))
            if not gene_id or gene_id.lower() == "na":
                gene_id = clean(row.get("Gene Symbol"), "unknown_gene")
            gene_name = clean(row.get("Gene Symbol") or row.get("Gene Name"))
            score = clean(row.get("Mean RPM") or row.get("Score"), "0")
        elif is_v4:
            pas_id = clean(row.get("PAS_ID"))
            parts = pas_id.split(":")
            if len(parts) < 3:
                continue
            chrom  = normalize_chrom(parts[0])
            strand = normalize_strand(parts[1])
            pos    = safe_int(parts[2])
            gene_name = clean(row.get("GeneSymbol"))
            gene_id   = gene_name or "unknown_gene"
            score = clean(row.get("PolyaStrength_percentile"), "0")
        else:
            chrom = normalize_chrom(row.get("Chr") or row.get("chrom"))
            strand = normalize_strand(row.get("Strand") or row.get("strand"))
            pos = safe_int(
                row.get("Position")
                or row.get("position")
                or row.get("start")
                or row.get("chromStart")
            )
            gene_id = clean(row.get("Gene") or row.get("gene_id") or row.get("gene"), "unknown_gene")
            gene_name = clean(row.get("Alias") or row.get("gene_name") or row.get("GeneSymbol"))
            score = clean(row.get("Score") or row.get("score"), "0")

        if not chrom:
            continue
        rows_out.append(
            {
                "gene_id": gene_id,
                "gene_name": gene_name,
                "chrom": chrom,
                "start": str(pos),
                "end": str(pos),
                "strand": strand,
                "score": score,
                "source": "polyadb3.2",
            }
        )

    logging.info("PolyA_DB: normalized %d rows", len(rows_out))
    return rows_out


def normalize_polyasite(path: Path) -> list[dict[str, str]]:
    rows = read_rows_with_optional_header(read_text(path), BED_COLS)
    rows_out = []

    for row in rows:
        chrom = normalize_chrom(row.get("chrom") or row.get("Chr") or row.get("chr"))
        if not chrom:
            continue
        start = safe_int(row.get("start") or row.get("chromStart"))
        end = safe_int(row.get("end") or row.get("chromEnd"), start)
        rows_out.append(
            {
                "gene_id": clean(row.get("gene_id") or row.get("Gene") or row.get("gene_name"), "unknown_gene"),
                "gene_name": clean(row.get("gene_name") or row.get("Alias") or row.get("name")),
                "chrom": chrom,
                "start": str(start),
                "end": str(end),
                "strand": normalize_strand(row.get("strand") or row.get("Strand")),
                "score": clean(row.get("score") or row.get("Score"), "0"),
                "source": "polyasite3.0",
            }
        )

    logging.info("PolyASite: normalized %d rows", len(rows_out))
    return rows_out


def anchor(row: dict[str, str]) -> int:
    start = safe_int(row["start"])
    end = safe_int(row["end"], start)
    return end if row["strand"] == "+" else start


def sort_key(row: dict[str, str]):
    return (row["chrom"], safe_int(row["start"]), safe_int(row["end"]), row["strand"], row["gene_id"])


def source_priority(row: dict[str, str]) -> tuple[int, int, str]:
    source = row.get("source", "")
    if source.startswith("polyadb"):
        priority = 0
    elif source.startswith("polyasite"):
        priority = 1
    else:
        priority = 2
    return (priority, anchor(row), row.get("gene_id", ""))


def collapse_cluster(cluster: list[dict[str, str]]) -> dict[str, str]:
    representative = sorted(cluster, key=source_priority)[0].copy()
    sources = sorted({clean(row.get("source"), "unknown") for row in cluster})
    representative["source"] = "+".join(sources)

    gene_ids = [clean(representative.get("gene_id"))] if clean(representative.get("gene_id")) else []
    gene_names = [clean(representative.get("gene_name"))] if clean(representative.get("gene_name")) else []
    for row in cluster:
        gene_id = clean(row.get("gene_id"))
        gene_name = clean(row.get("gene_name"))
        if gene_id and gene_id not in gene_ids:
            gene_ids.append(gene_id)
        if gene_name and gene_name not in gene_names:
            gene_names.append(gene_name)

    representative["gene_id"] = gene_ids[0] if gene_ids else "unknown_gene"
    representative["gene_name"] = gene_names[0] if gene_names else ""
    return representative


def merge_rows(rows: list[dict[str, str]], merge_distance: int) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault((row["chrom"], row["strand"]), []).append(row)

    merged = []
    for key in sorted(grouped):
        cluster: list[dict[str, str]] = []
        last_anchor = None
        for row in sorted(grouped[key], key=lambda item: (anchor(item), safe_int(item["start"]), item["gene_id"])):
            row_anchor = anchor(row)
            if not cluster or (last_anchor is not None and row_anchor - last_anchor <= merge_distance):
                cluster.append(row)
            else:
                merged.append(collapse_cluster(cluster))
                cluster = [row]
            last_anchor = row_anchor
        if cluster:
            merged.append(collapse_cluster(cluster))

    logging.info("Dual reference: merged %d rows into %d rows", len(rows), len(merged))
    return merged


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=FIELDNAMES_OUT,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(sorted(rows, key=sort_key))


def write_manifest(path: Path, manifest_rows: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["field", "value"])
        writer.writerows(manifest_rows)


def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
    )
    args = parse_args()
    cache_dir = Path(args.cache_dir)

    polyadb_path = resolve_source(
        spec=args.polya_db,
        staged_file=args.polya_db_file,
        default_url=args.polyadb_url,
        default_cache_name=POLYADB_CACHE_NAME,
        cache_dir=cache_dir,
        label="PolyA_DB",
    )
    polyasite_path = resolve_source(
        spec=args.polyasite,
        staged_file=args.polyasite_file,
        default_url=args.polyasite_url,
        default_cache_name=POLYASITE_CACHE_NAME,
        cache_dir=cache_dir,
        label="PolyASite",
    )

    requested_sources = [path for path in (polyadb_path, polyasite_path) if path is not None]
    source_rows: list[dict[str, str]] = []
    mode = "legacy"

    if requested_sources:
        mode = "single_source" if len(requested_sources) == 1 else "dual_reference"
        if polyadb_path is not None:
            source_rows.extend(normalize_polyadb(polyadb_path))
        if polyasite_path is not None:
            source_rows.extend(normalize_polyasite(polyasite_path))
        output_rows = merge_rows(source_rows, args.merge_distance) if len(requested_sources) == 2 else source_rows
    else:
        output_rows = normalize_known_polya(Path(args.known_polya))

    write_rows(Path(args.out_tsv), output_rows)
    write_manifest(
        Path(args.out_manifest),
        [
            ("mode", mode),
            ("known_polya", Path(args.known_polya).name),
            ("polya_db", str(polyadb_path or "skip")),
            ("polyasite", str(polyasite_path or "skip")),
            ("merge_distance", str(args.merge_distance)),
            ("input_rows", str(len(source_rows) if requested_sources else len(output_rows))),
            ("output_rows", str(len(output_rows))),
        ],
    )
    logging.info("Wrote %d rows to %s", len(output_rows), args.out_tsv)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        logging.exception("polyA reference preparation failed")
        sys.exit(1)
