#!/usr/bin/env python3
"""
apa_feature_extraction.py — Extract APA-relevant per-site features from
strand-aware grouped bedGraph coverage files.

For each candidate site in site_catalog, and for each group present in the
provided bedGraph directory, the following features are extracted:

  coverage_at_site       mean depth in ±window bp around the site
  upstream_cov           mean depth in [site-upstream, site-window)
  downstream_cov         mean depth in (site+window, site+window+upstream]
  read_end_density       mean depth in ±25 bp (proxy for 3' end density)
  proximal_distal_ratio  coverage_at_site / (downstream_cov + 0.01)
  dist_to_known_polya    distance (bp) to nearest known polyA site
  umi_support            per-site UMI count from site catalog (if present)
  cluster_support        number of groups with coverage ≥ 1 at this site

Output: apa_features.tsv
"""
import argparse
import glob
import logging
import os
import sys
from collections import defaultdict

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--site-catalog",     required=True)
    p.add_argument("--bedgraph-dir",     required=True)
    p.add_argument("--cell-annotations", required=True)
    p.add_argument("--known-polya",      required=True)
    p.add_argument("--out-features",     required=True)
    p.add_argument("--log",              required=True)
    p.add_argument("--window",           type=int, default=50)
    p.add_argument("--upstream",         type=int, default=200)
    return p.parse_args()


# ── bedGraph helpers ──────────────────────────────────────────────────────────

def load_bedgraph(path):
    """Parse bedGraph → dict {chrom: [(start, end, score), ...]}."""
    cov = defaultdict(list)
    try:
        with open(path) as fh:
            for line in fh:
                if line.startswith(("track", "#", "browser")):
                    continue
                parts = line.rstrip().split("\t")
                if len(parts) < 4:
                    continue
                cov[parts[0]].append((int(parts[1]), int(parts[2]), float(parts[3])))
    except Exception:
        pass
    return cov


def mean_coverage(cov_dict, chrom, start, end):
    """Mean coverage over half-open interval [start, end)."""
    if chrom not in cov_dict or start >= end:
        return 0.0
    span = end - start
    total = 0.0
    for seg_s, seg_e, score in cov_dict[chrom]:
        overlap = min(seg_e, end) - max(seg_s, start)
        if overlap > 0:
            total += score * overlap
    return total / span


# ── known polyA helpers ───────────────────────────────────────────────────────

def load_known_polya(path):
    """Load known polyA sites → list of (chrom, pos, strand)."""
    sites = []
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return sites
    try:
        df = pd.read_csv(path, sep="\t", header=None, comment="#")
        for _, row in df.iterrows():
            strand = str(row.iloc[2]) if df.shape[1] > 2 else "."
            sites.append((str(row.iloc[0]), int(row.iloc[1]), strand))
    except Exception:
        pass
    return sites


def nearest_polya_dist(chrom, pos, strand, known_sites):
    """Distance (bp) to the nearest known polyA site on the same chrom/strand."""
    min_dist = 10_000_000
    for k_chrom, k_pos, k_strand in known_sites:
        if k_chrom == chrom and (k_strand in (".", strand)):
            min_dist = min(min_dist, abs(pos - k_pos))
    return min_dist


# ── site parsing ──────────────────────────────────────────────────────────────

def parse_site_id(site_id):
    """Parse 'gene:chrom:pos:strand' → (chrom, pos, strand)."""
    parts = str(site_id).split(":")
    if len(parts) >= 4:
        try:
            return parts[1], int(parts[2]), parts[3]
        except ValueError:
            pass
    return "chr1", 0, "+"


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    logging.basicConfig(
        filename=args.log, level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s")
    log = logging.getLogger()
    log.info("Starting APA feature extraction")

    # Load site catalog
    try:
        catalog = pd.read_csv(args.site_catalog, sep="\t")
    except Exception as e:
        log.error(f"Cannot load site catalog: {e}. Writing empty output.")
        pd.DataFrame(columns=[
            "site_id", "group_level", "group_id", "coverage_at_site",
            "upstream_cov", "downstream_cov", "read_end_density",
            "proximal_distal_ratio", "dist_to_known_polya",
            "umi_support", "cluster_support",
        ]).to_csv(args.out_features, sep="\t", index=False)
        return

    if "site_id" not in catalog.columns:
        catalog["site_id"] = catalog.index.astype(str)
        log.warning("site_catalog has no 'site_id' column — using row index")

    # Load known polyA sites
    known_sites = load_known_polya(args.known_polya)
    log.info(f"Known polyA sites loaded: {len(known_sites)}")

    # Derive group levels from annotations
    try:
        anno = pd.read_csv(args.cell_annotations, sep="\t")
        group_levels = [c for c in ("cluster", "cell_type") if c in anno.columns]
    except Exception:
        group_levels = ["cluster"]
    log.info(f"Group levels: {group_levels}")

    # Discover bedGraph files
    fwd_files = glob.glob(os.path.join(args.bedgraph_dir, "*.fwd.bedGraph"))
    rev_files = glob.glob(os.path.join(args.bedgraph_dir, "*.rev.bedGraph"))
    log.info(f"bedGraph files — fwd: {len(fwd_files)}, rev: {len(rev_files)}")

    # Load into memory keyed by (group_id, strand)
    bg_data: dict = {}
    for path in fwd_files:
        gid = os.path.basename(path).replace(".fwd.bedGraph", "")
        bg_data[(gid, "+")] = load_bedgraph(path)
    for path in rev_files:
        gid = os.path.basename(path).replace(".rev.bedGraph", "")
        bg_data[(gid, "-")] = load_bedgraph(path)

    all_group_ids = sorted({k[0] for k in bg_data.keys()})
    log.info(f"Groups found: {len(all_group_ids)}")

    w = args.window
    u = args.upstream
    rows = []

    for _, site_row in catalog.iterrows():
        sid = site_row["site_id"]
        chrom, pos, strand = parse_site_id(sid)
        dist_known = nearest_polya_dist(chrom, pos, strand, known_sites)
        umi_support = int(site_row.get("umi_count", 0)) if "umi_count" in site_row else 0

        is_covered = 0
        for gid in all_group_ids:
            cov_dict = bg_data.get((gid, strand), {})

            at_site  = mean_coverage(cov_dict, chrom, pos - w, pos + w)
            up_cov   = mean_coverage(cov_dict, chrom, max(0, pos - u - w), pos - w)
            dn_cov   = mean_coverage(cov_dict, chrom, pos + w, pos + w + u)
            end_den  = mean_coverage(cov_dict, chrom, pos - 25, pos + 25)
            pd_ratio = at_site / (dn_cov + 0.01)

            if at_site > 0:
                is_covered += 1

            rows.append({
                "site_id":               sid,
                "group_level":           group_levels[0] if group_levels else "cluster",
                "group_id":              gid,
                "coverage_at_site":      round(at_site,  4),
                "upstream_cov":          round(up_cov,   4),
                "downstream_cov":        round(dn_cov,   4),
                "read_end_density":      round(end_den,  4),
                "proximal_distal_ratio": round(pd_ratio, 4),
                "dist_to_known_polya":   dist_known,
                "umi_support":           umi_support,
                "cluster_support":       0,  # filled below
            })

        # Backfill cluster_support for this site's rows
        for row in rows[-len(all_group_ids):]:
            row["cluster_support"] = is_covered

    df_out = pd.DataFrame(rows)
    df_out.to_csv(args.out_features, sep="\t", index=False)
    log.info(f"Wrote {len(df_out)} feature rows → {args.out_features}")


if __name__ == "__main__":
    main()
