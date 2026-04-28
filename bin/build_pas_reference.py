#!/usr/bin/env python3

import argparse
import csv
import logging
import re
from collections import defaultdict
from pathlib import Path


REQUIRED_SITE_COLUMNS = {
    "site_id",
    "gene_id",
    "chrom",
    "start",
    "end",
    "strand",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a contract-stable PAS reference from repo-native inputs."
    )
    parser.add_argument("--site-catalog", required=True)
    parser.add_argument("--terminal-exons", required=True)
    parser.add_argument("--known-polya", required=True)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--out-manifest", required=True)
    parser.add_argument("--log", required=True)
    parser.add_argument("--out-bed")
    parser.add_argument("--merge-distance", type=int, default=25)
    return parser.parse_args()


def setup_logging(log_path: Path):
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.FileHandler(log_path, mode="w", encoding="utf-8"),
            logging.StreamHandler(),
        ],
    )


def is_missing_file(path: Path) -> bool:
    return (not path.exists()) or path.name == "NO_FILE" or path.stat().st_size == 0


def load_tsv_rows(path: Path):
    if is_missing_file(path):
        return []
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def normalize_chrom(value: str) -> str:
    chrom = str(value or "").strip()
    if not chrom:
        return chrom
    if chrom in {"M", "MT"}:
        return "chrM"
    if chrom.startswith("chr"):
        return chrom
    return f"chr{chrom}"


def safe_int(value, default=0):
    if value is None:
        return default
    text = str(value).strip()
    if not text:
        return default
    return int(float(text))


def anchor_for_interval(start: int, end: int, strand: str) -> int:
    return end if strand == "+" else start


def normalize_site_catalog(path: Path):
    rows = load_tsv_rows(path)
    normalized = []
    for row in rows:
        if not REQUIRED_SITE_COLUMNS.issubset(row.keys()):
            missing = sorted(REQUIRED_SITE_COLUMNS - set(row.keys()))
            raise ValueError(f"site_catalog missing required columns: {missing}")
        start = safe_int(row.get("start"))
        end = safe_int(row.get("end"), start)
        strand = (row.get("strand") or "+").strip() or "+"
        normalized.append(
            {
                "gene_id": (row.get("gene_id") or "unknown_gene").strip(),
                "chrom": normalize_chrom(row.get("chrom")),
                "start": start,
                "end": end,
                "strand": strand,
                "anchor": anchor_for_interval(start, end, strand),
                "site_id": (row.get("site_id") or "").strip(),
                "source_kind": "site_catalog",
                "site_source": (row.get("site_source") or "site_catalog").strip(),
            }
        )
    return normalized


def normalize_terminal_exons(path: Path):
    rows = load_tsv_rows(path)
    normalized = []
    for row in rows:
        start = safe_int(row.get("start"))
        end = safe_int(row.get("end"), start)
        strand = (row.get("strand") or "+").strip() or "+"
        anchor = anchor_for_interval(start, end, strand)
        gene_id = (row.get("gene_id") or "unknown_gene").strip()
        chrom = normalize_chrom(row.get("chrom"))
        normalized.append(
            {
                "gene_id": gene_id,
                "chrom": chrom,
                "start": anchor,
                "end": anchor,
                "strand": strand,
                "anchor": anchor,
                "site_id": f"{gene_id}:{chrom}:{anchor}:{strand}",
                "source_kind": "terminal_exon",
                "site_source": "terminal_exon",
            }
        )
    return normalized


def normalize_known_polya(path: Path):
    rows = load_tsv_rows(path)
    normalized = []
    for row in rows:
        # gene_id: canonical column first, then PolyA_DB 'Gene'/'Alias', then gene_name
        gene_id = (
            row.get("gene_id")
            or row.get("Gene")
            or row.get("gene")
            or row.get("gene_name")
            or "unknown_gene"
        ).strip()
        # chrom: BED uses 'chrom', PolyA_DB uses 'Chr'/'chr'
        chrom = normalize_chrom(
            row.get("chrom") or row.get("Chr") or row.get("chr")
        )
        # start/end: canonical names first, then PolyASite BED chromStart/chromEnd,
        # then PolyA_DB single-base 'Position'/'position'
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
        strand = (row.get("strand") or row.get("Strand") or "+").strip() or "+"
        anchor = anchor_for_interval(start, end, strand)
        normalized.append(
            {
                "gene_id": gene_id,
                "chrom": chrom,
                "start": anchor,
                "end": anchor,
                "strand": strand,
                "anchor": anchor,
                "site_id": f"{gene_id}:{chrom}:{anchor}:{strand}",
                "source_kind": "known_polya",
                "site_source": (row.get("source") or row.get("reference_source") or "known_polya").strip(),
            }
        )
    return normalized


def representative_priority(row):
    if row["source_kind"] == "site_catalog" and "atlas" in row["site_source"]:
        return 0
    if row["source_kind"] == "known_polya":
        return 1
    if row["source_kind"] == "site_catalog":
        return 2
    return 3


def choose_representative(rows):
    strand = rows[0]["strand"]
    if strand == "+":
        coord_key = lambda row: (row["anchor"], row["start"], row["end"])
    else:
        coord_key = lambda row: (-row["anchor"], -row["start"], -row["end"])
    return sorted(rows, key=lambda row: (representative_priority(row), coord_key(row), row["site_id"]))[0]


def evidence_labels(row):
    if row["source_kind"] == "known_polya":
        labels = [
            label.strip()
            for label in (row.get("site_source") or "known_polya").split("+")
            if label.strip()
        ]
        return labels or ["known_polya"]
    return [row["source_kind"]]


def cluster_candidates(candidates, merge_distance):
    grouped = defaultdict(list)
    for row in candidates:
        key = (row["gene_id"], row["chrom"], row["strand"])
        grouped[key].append(row)

    clusters = []
    for key in sorted(grouped):
        rows = sorted(grouped[key], key=lambda row: row["anchor"])
        current = []
        for row in rows:
            if not current:
                current = [row]
                continue
            if abs(row["anchor"] - current[-1]["anchor"]) <= merge_distance:
                current.append(row)
            else:
                clusters.append(current)
                current = [row]
        if current:
            clusters.append(current)
    return clusters


def sanitize_token(text: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", text or "unknown")
    return cleaned.strip("_") or "unknown"


def build_reference_rows(site_catalog_rows, terminal_exon_rows, known_polya_rows, merge_distance):
    candidates = site_catalog_rows + known_polya_rows + terminal_exon_rows
    if not candidates:
        return []

    output_rows = []
    for cluster in cluster_candidates(candidates, merge_distance):
        rep = choose_representative(cluster)
        anchor = rep["anchor"]
        sources = sorted({label for row in cluster for label in evidence_labels(row)})
        source_label = "+".join(sources)
        gene_token = sanitize_token(rep["gene_id"])
        chrom_token = sanitize_token(rep["chrom"])
        output_rows.append(
            {
                "pas_reference_id": f"pasref_{gene_token}_{chrom_token}_{anchor}_{rep['strand']}",
                "site_id": rep["site_id"] or f"{rep['gene_id']}:{rep['chrom']}:{anchor}:{rep['strand']}",
                "gene_id": rep["gene_id"],
                "chrom": rep["chrom"],
                "start": anchor,
                "end": anchor,
                "strand": rep["strand"],
                "reference_source": source_label,
            }
        )

    output_rows.sort(key=lambda row: (row["gene_id"], row["chrom"], row["start"], row["strand"]))
    return output_rows


def write_reference_rows(path: Path, rows):
    fieldnames = [
        "pas_reference_id",
        "site_id",
        "gene_id",
        "chrom",
        "start",
        "end",
        "strand",
        "reference_source",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_bed(path: Path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for row in rows:
            writer.writerow(
                [
                    row["chrom"],
                    max(int(row["start"]) - 1, 0),
                    row["end"],
                    row["pas_reference_id"],
                    0,
                    row["strand"],
                ]
            )


def write_manifest(path: Path, args, site_catalog_rows, terminal_exon_rows, known_polya_rows, output_rows):
    manifest_rows = [
        ("site_catalog", Path(args.site_catalog).name),
        ("terminal_exons", Path(args.terminal_exons).name),
        ("known_polya", Path(args.known_polya).name),
        ("merge_distance", str(args.merge_distance)),
        ("site_catalog_rows", str(len(site_catalog_rows))),
        ("terminal_exon_rows", str(len(terminal_exon_rows))),
        ("known_polya_rows", str(len(known_polya_rows))),
        ("pas_reference_rows", str(len(output_rows))),
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["field", "value"])
        writer.writerows(manifest_rows)


def main():
    args = parse_args()
    log_path = Path(args.log)
    setup_logging(log_path)

    try:
        site_catalog_path = Path(args.site_catalog)
        terminal_exons_path = Path(args.terminal_exons)
        known_polya_path = Path(args.known_polya)

        logging.info("Loading PAS reference inputs")
        site_catalog_rows = normalize_site_catalog(site_catalog_path)
        terminal_exon_rows = normalize_terminal_exons(terminal_exons_path)
        known_polya_rows = normalize_known_polya(known_polya_path)

        logging.info(
            "Loaded %s site_catalog rows, %s terminal_exon rows, %s known_polya rows",
            len(site_catalog_rows),
            len(terminal_exon_rows),
            len(known_polya_rows),
        )

        output_rows = build_reference_rows(
            site_catalog_rows,
            terminal_exon_rows,
            known_polya_rows,
            args.merge_distance,
        )

        out_tsv = Path(args.out_tsv)
        write_reference_rows(out_tsv, output_rows)
        write_manifest(
            Path(args.out_manifest),
            args,
            site_catalog_rows,
            terminal_exon_rows,
            known_polya_rows,
            output_rows,
        )

        if args.out_bed:
            write_bed(Path(args.out_bed), output_rows)

        logging.info("Wrote %s PAS reference rows to %s", len(output_rows), out_tsv)
    except Exception:
        logging.exception("PAS reference construction failed")
        raise


if __name__ == "__main__":
    main()
