#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def count_rows(path: Path):
    if not path.exists() or path.stat().st_size == 0:
        return 0
    with path.open("r", encoding="utf-8") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def flatten_paths(items):
    flattened = []
    for item in items:
        path = Path(item)
        if path.is_dir():
            flattened.extend(sorted(str(child) for child in path.iterdir()))
        else:
            flattened.append(str(path))
    return flattened


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--site-catalog", required=True)
    parser.add_argument("--apa-usage", required=True)
    parser.add_argument("--apa-stats", required=True)
    parser.add_argument("--track-paths", nargs="*")
    parser.add_argument("--qc-paths", nargs="*")
    parser.add_argument("--out-html", required=True)
    parser.add_argument("--out-plots", required=True)
    parser.add_argument("--out-summary", required=True)
    args = parser.parse_args()

    track_paths = flatten_paths(args.track_paths or [])
    qc_paths = flatten_paths(args.qc_paths or [])

    plots_dir = Path(args.out_plots)
    plots_dir.mkdir(parents=True, exist_ok=True)
    (plots_dir / "summary.tsv").write_text(
        "metric\tvalue\n"
        f"site_catalog_rows\t{count_rows(Path(args.site_catalog))}\n"
        f"apa_usage_rows\t{count_rows(Path(args.apa_usage))}\n"
        f"apa_stats_rows\t{count_rows(Path(args.apa_stats))}\n",
        encoding="utf-8",
    )

    track_text = "".join(path + "\n" for path in track_paths)
    qc_text = "".join(path + "\n" for path in qc_paths)
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>scPolASeq report</title>
</head>
<body>
  <h1>scPolASeq report</h1>
  <ul>
    <li>site catalog rows: {count_rows(Path(args.site_catalog))}</li>
    <li>APA usage rows: {count_rows(Path(args.apa_usage))}</li>
    <li>APA stats rows: {count_rows(Path(args.apa_stats))}</li>
    <li>track entries: {len(track_paths)}</li>
    <li>QC entries: {len(qc_paths)}</li>
  </ul>
  <h2>Tracks</h2>
  <pre>{track_text}</pre>
  <h2>QC files</h2>
  <pre>{qc_text}</pre>
</body>
</html>
"""
    Path(args.out_html).write_text(html, encoding="utf-8")

    with Path(args.out_summary).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["artifact_type", "count"])
        writer.writerow(["track_paths", len(track_paths)])
        writer.writerow(["qc_paths", len(qc_paths)])


if __name__ == "__main__":
    main()
