#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict
from pathlib import Path


WINDOW = 50


def load_sites(path: Path):
    sites = defaultdict(list)
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = (row["chrom"], row["strand"])
            row["start"] = int(row["start"])
            row["end"] = int(row["end"])
            sites[key].append(row)
    return sites


def nearest_site(site_rows, position):
    best = None
    best_distance = None
    for site in site_rows:
        anchor = site["end"] if site["strand"] == "+" else site["start"]
        distance = abs(position - anchor)
        if distance <= WINDOW and (best_distance is None or distance < best_distance):
            best = site
            best_distance = distance
    return best


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coverage-files", nargs="+", required=True)
    parser.add_argument("--site-catalog", required=True)
    parser.add_argument("--group-map", required=True)
    parser.add_argument("--cell-annotations", required=True)
    parser.add_argument("--out-usage", required=True)
    parser.add_argument("--out-matrix", required=True)
    parser.add_argument("--out-gene-summary", required=True)
    args = parser.parse_args()

    sites = load_sites(Path(args.site_catalog))
    site_counts = defaultdict(lambda: {"read_count": 0, "umi_count": 0})

    for coverage_file in args.coverage_files:
        with Path(coverage_file).open("r", newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                key = (row["chrom"], row["strand"])
                position = int(row["position"])
                site = nearest_site(sites.get(key, []), position)
                if site is None:
                    continue
                count_key = (row["group_level"], row["group_id"], site["gene_id"], site["site_id"])
                site_counts[count_key]["read_count"] += int(row["read_count"])
                site_counts[count_key]["umi_count"] += int(row["umi_count"])

    gene_totals = defaultdict(int)
    for (group_level, group_id, gene_id, site_id), value in site_counts.items():
        gene_totals[(group_level, group_id, gene_id)] += value["umi_count"]

    usage_rows = []
    matrix_rows = defaultdict(dict)
    gene_summary = defaultdict(lambda: {"sites": 0, "umis": 0})
    for (group_level, group_id, gene_id, site_id), value in sorted(site_counts.items()):
        total_umis = gene_totals[(group_level, group_id, gene_id)] or 1
        usage_fraction = value["umi_count"] / total_umis
        usage_row = {
            "group_level": group_level,
            "group_id": group_id,
            "gene_id": gene_id,
            "site_id": site_id,
            "read_count": value["read_count"],
            "umi_count": value["umi_count"],
            "usage_fraction": f"{usage_fraction:.6f}",
            "pdui_like": f"{usage_fraction:.6f}",
        }
        usage_rows.append(usage_row)
        matrix_rows[site_id][group_id] = f"{usage_fraction:.6f}"
        gene_summary[(group_level, group_id, gene_id)]["sites"] += 1
        gene_summary[(group_level, group_id, gene_id)]["umis"] += value["umi_count"]

    with Path(args.out_usage).open("w", newline="", encoding="utf-8") as handle:
        fieldnames = ["group_level", "group_id", "gene_id", "site_id", "read_count", "umi_count", "usage_fraction", "pdui_like"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(usage_rows)

    group_ids = sorted({row["group_id"] for row in usage_rows})
    with Path(args.out_matrix).open("w", newline="", encoding="utf-8") as handle:
        fieldnames = ["site_id"] + group_ids
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for site_id in sorted(matrix_rows):
            row = {"site_id": site_id}
            row.update(matrix_rows[site_id])
            writer.writerow(row)

    with Path(args.out_gene_summary).open("w", newline="", encoding="utf-8") as handle:
        fieldnames = ["group_level", "group_id", "gene_id", "n_sites", "total_umi_count"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for (group_level, group_id, gene_id), summary in sorted(gene_summary.items()):
            writer.writerow(
                {
                    "group_level": group_level,
                    "group_id": group_id,
                    "gene_id": gene_id,
                    "n_sites": summary["sites"],
                    "total_umi_count": summary["umis"],
                }
            )


if __name__ == "__main__":
    main()
