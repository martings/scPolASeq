#!/usr/bin/env python3
"""
apa_feature_extraction.py — Extract APA-relevant per-site features from
strand-aware grouped bedGraph coverage files.

OPTIMIZED VERSION:
  • bedGraphs loaded once with pandas + converted to numpy arrays per chrom
  • Mean coverage over an interval computed in O(log N) via np.searchsorted
    on cumulative-coverage arrays (prefix sums), not a linear scan
  • Known polyA distances computed per (chrom, strand) with np.searchsorted
  • Output is streamed to disk in chunks so memory stays flat
  • Progress is logged every N sites so you know it's alive

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
import bisect
import glob
import logging
import os
import time
from collections import defaultdict

import numpy as np
import pandas as pd


# ── args ──────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--site-catalog",     required=True)
    p.add_argument("--bedgraph-dir",     required=True)
    p.add_argument("--bedgraph-glob",    default="*.bedGraph",
                   help="(kept for CLI compatibility; not used)")
    p.add_argument("--cell-annotations", required=True)
    p.add_argument("--known-polya",      required=True)
    p.add_argument("--out-features",     required=True)
    p.add_argument("--log",              required=True)
    p.add_argument("--window",           type=int, default=50)
    p.add_argument("--upstream",         type=int, default=200)
    p.add_argument("--chunk-size",       type=int, default=50_000,
                   help="Flush output every N sites (default 50000)")
    p.add_argument("--progress-every",   type=int, default=10_000,
                   help="Log progress every N sites (default 10000)")
    return p.parse_args()


# ── bedGraph: load as prefix-sum arrays ───────────────────────────────────────

def load_bedgraph_prefix(path):
    """
    Read a bedGraph into per-chromosome numpy arrays:
        chrom -> (starts, ends, scores, cum)
    where cum[i] = sum_{k<i} scores[k] * (ends[k]-starts[k])

    With this, the *weighted* sum of score over any segment index range [i,j)
    is cum[j] - cum[i]. For partial-overlap segments at the edges we handle
    them explicitly in mean_coverage().

    Assumes each chromosome's intervals are non-overlapping and sorted by
    start — which is the standard bedGraph convention from tools like
    bedtools genomecov / deepTools.
    """
    try:
        df = pd.read_csv(
            path, sep="\t", header=None, comment="#",
            names=["chrom", "start", "end", "score"],
            dtype={"chrom": str, "start": np.int64,
                   "end": np.int64, "score": np.float32},
            engine="c",
        )
    except Exception:
        return {}

    # drop any stray track/browser lines that slipped past comment filter
    df = df.dropna()
    if df.empty:
        return {}

    out = {}
    # groupby chrom, assume already sorted inside each chrom; sort defensively
    for chrom, sub in df.groupby("chrom", sort=False):
        sub = sub.sort_values("start", kind="mergesort")
        starts = sub["start"].to_numpy(np.int64, copy=False)
        ends   = sub["end"].to_numpy(np.int64,   copy=False)
        scores = sub["score"].to_numpy(np.float32, copy=False)
        # prefix sum of score * length, in float64 to avoid precision loss
        weighted = (scores.astype(np.float64) * (ends - starts))
        cum = np.concatenate(([0.0], np.cumsum(weighted)))
        out[chrom] = (starts, ends, scores, cum)
    return out


def mean_coverage(cov_entry, start, end):
    """
    Mean coverage on half-open [start, end) using prefix sums.

    cov_entry = (starts, ends, scores, cum) for one chromosome, or None.
    O(log N) per call thanks to searchsorted.
    """
    if cov_entry is None or start >= end:
        return 0.0
    starts, ends, scores, cum = cov_entry

    # Find the range of segments that could overlap [start, end):
    #   left:  first segment with end > start   → searchsorted on 'ends'
    #   right: first segment with start >= end  → searchsorted on 'starts'
    left  = np.searchsorted(ends,   start, side="right")
    right = np.searchsorted(starts, end,   side="left")
    if left >= right:
        return 0.0

    # Bulk contribution from fully-contained segments [left+1, right-1]
    # (we'll correct the edges separately)
    total = cum[right] - cum[left]

    # Subtract the portion of segment `left` that lies before `start`
    seg_s = starts[left]
    seg_e = ends[left]
    if seg_s < start:
        total -= float(scores[left]) * (start - seg_s)

    # Subtract the portion of segment `right-1` that lies at/after `end`
    last = right - 1
    seg_e_last = ends[last]
    if seg_e_last > end:
        total -= float(scores[last]) * (seg_e_last - end)

    return total / (end - start)


# ── known polyA: per-(chrom,strand) sorted arrays ─────────────────────────────

def load_known_polya_indexed(path):
    """Build {(chrom, strand): sorted np.ndarray of positions}.

    Handles two formats:
    - Simple 3-col: chrom\\tpos\\tstrand  (no header)
    - PolyASite2:   header present, columns chrom/chromStart/chromEnd/name/score/strand/...
                    uses chromEnd as the representative cleavage position.
                    Chromosomes without 'chr' prefix are normalised (e.g. '10' -> 'chr10').
    """
    index: dict = {}
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return index
    try:
        df = pd.read_csv(path, sep="\t", low_memory=False)
    except Exception:
        return index
    if df.empty:
        return index

    cols = list(df.columns)

    # PolyASite2 format: named columns
    if "chromEnd" in cols and "strand" in cols:
        chrom_col  = df["chrom"].astype(str)
        pos_col    = pd.to_numeric(df["chromEnd"], errors="coerce")
        strand_col = df["strand"].astype(str)
    else:
        # Fallback: first three columns as chrom, pos, strand
        chrom_col  = df.iloc[:, 0].astype(str)
        pos_col    = pd.to_numeric(df.iloc[:, 1], errors="coerce")
        strand_col = df.iloc[:, 2].astype(str) if len(cols) > 2 \
                     else pd.Series(["."] * len(df))

    tmp = defaultdict(list)
    for c, p, s in zip(chrom_col, pos_col, strand_col):
        if pd.isna(p):
            continue
        # Normalise: add 'chr' prefix if missing
        c = str(c)
        if not c.startswith("chr"):
            c = "chr" + c
        p = int(p)
        if s in ("+", "-"):
            tmp[(c, s)].append(p)
        else:
            tmp[(c, "+")].append(p)
            tmp[(c, "-")].append(p)

    for k, vals in tmp.items():
        index[k] = np.sort(np.asarray(vals, dtype=np.int64))
    return index


def nearest_polya_dist(polya_index, chrom, pos, strand):
    arr = polya_index.get((chrom, strand))
    if arr is None or arr.size == 0:
        return 10_000_000
    i = np.searchsorted(arr, pos)
    best = 10_000_000
    if i < arr.size:
        best = int(abs(arr[i] - pos))
    if i > 0:
        best = min(best, int(abs(arr[i - 1] - pos)))
    return best


# ── site parsing ──────────────────────────────────────────────────────────────

def parse_site_id(site_id):
    """Parse 'gene:chrom:pos_or_range:strand' → (chrom, pos, strand).
    Position can be an integer ('100627080') or a range ('100627080-100627121'),
    in which case the midpoint is used.
    """
    parts = str(site_id).split(":")
    if len(parts) >= 4:
        try:
            raw_pos = parts[2]
            if "-" in raw_pos:
                start, end = raw_pos.split("-", 1)
                pos = (int(start) + int(end)) // 2
            else:
                pos = int(raw_pos)
            chrom = parts[1]
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            return chrom, pos, parts[3]
        except ValueError:
            pass
    return "chr1", 0, "+"


# ── main ──────────────────────────────────────────────────────────────────────

OUT_COLS = [
    "site_id", "group_level", "group_id",
    "coverage_at_site", "upstream_cov", "downstream_cov",
    "read_end_density", "proximal_distal_ratio",
    "dist_to_known_polya", "umi_support", "cluster_support",
]


def main():
    args = parse_args()
    logging.basicConfig(
        filename=args.log, level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s")
    log = logging.getLogger()
    log.info("Starting APA feature extraction (optimized)")
    t0 = time.time()

    # 1) Load site catalog --------------------------------------------------
    try:
        catalog = pd.read_csv(args.site_catalog, sep="\t", low_memory=False)
    except Exception as e:
        log.error(f"Cannot load site catalog: {e}. Writing empty output.")
        pd.DataFrame(columns=OUT_COLS).to_csv(
            args.out_features, sep="\t", index=False)
        return

    if "site_id" not in catalog.columns:
        catalog["site_id"] = catalog.index.astype(str)
        log.warning("site_catalog has no 'site_id' column — using row index")

    n_sites = len(catalog)
    log.info(f"Site catalog rows: {n_sites}")

    # 2) Known polyA index --------------------------------------------------
    polya_index = load_known_polya_indexed(args.known_polya)
    log.info(f"Known polyA index buckets: {len(polya_index)}")

    # 3) Group level --------------------------------------------------------
    # Map annotation column names → canonical group_level labels.
    # cluster_id is the standard column produced by the labeling step.
    _LEVEL_MAP = [("cluster", "cluster"), ("cluster_id", "cluster"), ("cell_type", "cell_type")]
    try:
        anno = pd.read_csv(args.cell_annotations, sep="\t")
        group_level = next(
            (label for col, label in _LEVEL_MAP if col in anno.columns),
            "cluster"
        )
    except Exception:
        group_level = "cluster"
    log.info(f"Group level: {group_level}")

    # 4) Discover + load bedGraphs -----------------------------------------
    fwd_files = glob.glob(os.path.join(args.bedgraph_dir, "**", "*.fwd.bedGraph"),
                          recursive=True)
    rev_files = glob.glob(os.path.join(args.bedgraph_dir, "**", "*.rev.bedGraph"),
                          recursive=True)
    if not fwd_files:
        fwd_files = glob.glob(os.path.join(args.bedgraph_dir, "*.fwd.bedGraph"))
        rev_files = glob.glob(os.path.join(args.bedgraph_dir, "*.rev.bedGraph"))
    log.info(f"bedGraph files — fwd: {len(fwd_files)}, rev: {len(rev_files)}")

    bg_data: dict = {}
    t_load = time.time()
    for path in fwd_files:
        gid = os.path.basename(path).replace(".fwd.bedGraph", "")
        bg_data[(gid, "+")] = load_bedgraph_prefix(path)
    for path in rev_files:
        gid = os.path.basename(path).replace(".rev.bedGraph", "")
        bg_data[(gid, "-")] = load_bedgraph_prefix(path)
    log.info(f"Loaded bedGraphs in {time.time() - t_load:.1f}s")

    all_group_ids = sorted({k[0] for k in bg_data.keys()})
    n_groups = len(all_group_ids)
    log.info(f"Groups found: {n_groups}")

    w = args.window
    u = args.upstream

    # 5) Prepare output (stream in chunks) ---------------------------------
    out_path = args.out_features
    # Write header first, then append
    pd.DataFrame(columns=OUT_COLS).to_csv(out_path, sep="\t", index=False)
    first_chunk_written = False  # header already there

    # Pre-extract columns once (faster than .iterrows)
    site_ids  = catalog["site_id"].astype(str).to_numpy()
    if "umi_count" in catalog.columns:
        umi_vals = pd.to_numeric(catalog["umi_count"], errors="coerce") \
                     .fillna(0).astype(np.int64).to_numpy()
    else:
        umi_vals = np.zeros(n_sites, dtype=np.int64)

    rows_buf = []
    total_rows_written = 0
    t_loop = time.time()

    # Cache per-(group,strand,chrom) entry lookups to avoid dict.get overhead
    for i in range(n_sites):
        sid = site_ids[i]
        chrom, pos, strand = parse_site_id(sid)
        dist_known = nearest_polya_dist(polya_index, chrom, pos, strand)
        umi_support = int(umi_vals[i])

        # Resolve per-group chrom entries once for this site
        is_covered = 0
        site_rows_start = len(rows_buf)

        # Intervals used for this site
        at_s,  at_e  = pos - w,           pos + w
        up_s,  up_e  = max(0, pos - u - w), pos - w
        dn_s,  dn_e  = pos + w,           pos + w + u
        en_s,  en_e  = pos - 25,          pos + 25

        for gid in all_group_ids:
            cov_dict = bg_data.get((gid, strand))
            cov_entry = cov_dict.get(chrom) if cov_dict else None

            at_site = mean_coverage(cov_entry, at_s, at_e)
            up_cov  = mean_coverage(cov_entry, up_s, up_e)
            dn_cov  = mean_coverage(cov_entry, dn_s, dn_e)
            end_den = mean_coverage(cov_entry, en_s, en_e)
            pd_ratio = at_site / (dn_cov + 0.01)

            if at_site > 0:
                is_covered += 1

            rows_buf.append((
                sid, group_level, gid,
                round(at_site,  4),
                round(up_cov,   4),
                round(dn_cov,   4),
                round(end_den,  4),
                round(pd_ratio, 4),
                dist_known,
                umi_support,
                0,  # cluster_support placeholder
            ))

        # Backfill cluster_support for this site's rows
        for k in range(site_rows_start, len(rows_buf)):
            r = rows_buf[k]
            rows_buf[k] = r[:-1] + (is_covered,)

        # Flush chunk
        if (i + 1) % args.chunk_size == 0:
            df_chunk = pd.DataFrame(rows_buf, columns=OUT_COLS)
            df_chunk.to_csv(out_path, sep="\t", index=False,
                            header=False, mode="a")
            total_rows_written += len(df_chunk)
            rows_buf.clear()

        # Progress
        if (i + 1) % args.progress_every == 0:
            elapsed = time.time() - t_loop
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            eta = (n_sites - (i + 1)) / rate if rate > 0 else float("inf")
            log.info(
                f"Processed {i + 1}/{n_sites} sites "
                f"({100.0 * (i + 1) / n_sites:.1f}%) — "
                f"{rate:.0f} sites/s — ETA {eta / 60:.1f} min"
            )

    # Final flush
    if rows_buf:
        df_chunk = pd.DataFrame(rows_buf, columns=OUT_COLS)
        df_chunk.to_csv(out_path, sep="\t", index=False,
                        header=False, mode="a")
        total_rows_written += len(df_chunk)
        rows_buf.clear()

    log.info(
        f"Wrote {total_rows_written} feature rows → {out_path} "
        f"in {time.time() - t0:.1f}s total"
    )


if __name__ == "__main__":
    main()