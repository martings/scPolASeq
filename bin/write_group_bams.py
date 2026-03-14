#!/usr/bin/env python3

import argparse
import csv
import re
from pathlib import Path


def sanitize(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--group-map", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--out-manifest", required=True)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    groups = set()
    with Path(args.group_map).open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("sample_id") == args.sample_id:
                groups.add((row.get("group_level", ""), row.get("group_id", "")))

    manifest_rows = []
    for group_level, group_id in sorted(groups):
        placeholder = out_dir / f"{sanitize(group_level)}__{sanitize(group_id)}.bam"
        placeholder.write_text("", encoding="utf-8")
        manifest_rows.append((group_level, group_id, str(placeholder.name)))

    with Path(args.out_manifest).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["group_level", "group_id", "placeholder_bam"])
        writer.writerows(manifest_rows)


if __name__ == "__main__":
    main()
