#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def count_rows(path: Path):
    if not path.exists() or path.stat().st_size == 0:
        return 0
    with path.open("r", encoding="utf-8") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-table", required=True)
    parser.add_argument("--manifest-name", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    input_path = Path(args.input_table)
    with Path(args.out).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["manifest_name", "input_table", "rows"])
        writer.writerow([args.manifest_name, str(input_path.resolve()), count_rows(input_path)])


if __name__ == "__main__":
    main()
