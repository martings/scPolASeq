#!/usr/bin/env python3

import argparse
import csv
import json
import re
from pathlib import Path


REQUIRED_COLUMNS = [
    "sample_id",
    "library_id",
    "protocol",
    "chemistry",
    "condition",
]

OPTIONAL_COLUMNS = [
    "replicate_id",
    "fastq_r1",
    "fastq_r2",
    "bam",
    "matrix_path",
    "barcode_whitelist",
]

ALLOWED_PROTOCOLS = {"10x_3p", "10x_5p", "bam_only"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--schema", required=False)
    parser.add_argument("--out", required=True)
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args()

    input_path = Path(args.input)
    rows = []
    errors = []
    canonical_library_ids = {}

    with input_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise SystemExit("Input samplesheet has no header row")

        missing = [column for column in REQUIRED_COLUMNS if column not in reader.fieldnames]
        if missing:
            raise SystemExit(f"Samplesheet missing required columns: {', '.join(missing)}")

        for index, row in enumerate(reader, start=2):
            clean = {key: (value or "").strip() for key, value in row.items()}
            for column in REQUIRED_COLUMNS:
                if not clean.get(column):
                    errors.append(f"Line {index}: missing required value for {column}")

            library_id = clean.get("library_id", "")
            if library_id:
                canonical_id = re.sub(r"[^A-Za-z0-9._-]+", "_", library_id)
                previous = canonical_library_ids.get(canonical_id)
                if previous is not None and previous != library_id:
                    errors.append(
                        f"Line {index}: library_id {library_id!r} collides with "
                        f"{previous!r} after filename sanitization ({canonical_id!r})"
                    )
                canonical_library_ids[canonical_id] = library_id

            protocol = clean.get("protocol", "")
            if protocol and protocol not in ALLOWED_PROTOCOLS:
                errors.append(f"Line {index}: unsupported protocol '{protocol}'")

            has_fastq = bool(clean.get("fastq_r1") and clean.get("fastq_r2"))
            has_bam = bool(clean.get("bam"))
            if not (has_fastq or has_bam):
                errors.append(f"Line {index}: requires fastq_r1/fastq_r2 or bam")

            normalized = {column: clean.get(column, "") for column in REQUIRED_COLUMNS + OPTIONAL_COLUMNS}
            rows.append(normalized)

    if errors:
        raise SystemExit("\n".join(errors))

    output_columns = REQUIRED_COLUMNS + OPTIONAL_COLUMNS
    with Path(args.out).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_columns)
        writer.writeheader()
        writer.writerows(rows)

    manifest = {
        "input": str(input_path.resolve()),
        "validated_rows": len(rows),
        "schema": args.schema,
        "required_columns": REQUIRED_COLUMNS,
        "optional_columns": OPTIONAL_COLUMNS,
    }
    Path(args.manifest).write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
