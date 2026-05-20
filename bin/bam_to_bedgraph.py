#!/usr/bin/env python3
"""
bam_to_bedgraph.py — Strand-specific coverage in bedGraph format using pysam.

Replaces bedtools genomecov for the COVERAGE_TRACKS module.
Outputs 0-based half-open intervals [start, end) with integer coverage.

Optimizations vs. the pileup-based version:
- Uses bam.fetch() + NumPy arrays instead of bam.pileup() (much faster: one pass per read,
  not one pass per base).
- Splits large chromosomes into chunks for better load balancing across workers.
- Vectorized run-length encoding for bedGraph emission.
- Forward and reverse bedGraphs written in parallel.

Usage (single pass, both strands — preferred):
    python bam_to_bedgraph.py input.bam --fwd-output out.fwd.bedGraph --rev-output out.rev.bedGraph [--cpus 4]

Usage (single strand — kept for backward compatibility):
    python bam_to_bedgraph.py input.bam --strand + --output out.fwd.bedGraph [--cpus 4]
"""
import argparse
import gzip
import sys
import logging
from pathlib import Path
from multiprocessing import Pool, Process

import numpy as np

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger()

# Chunk size for splitting large chromosomes (in bp). 10 Mb is a good balance:
# - small enough to balance load across workers
# - large enough that per-chunk overhead (BAM open, fetch setup) is negligible
CHUNK_SIZE = 10_000_000


def _coverage_chunk(args):
    """Worker: compute strand-specific coverage for a chromosomal region.

    Returns (chrom, region_start, cov_fwd, cov_rev) where cov_fwd/cov_rev are
    NumPy int32 arrays of length (region_end - region_start), giving coverage
    over the region [region_start, region_end).
    """
    bam_path, chrom, region_start, region_end, min_mapq = args
    import pysam

    region_len = region_end - region_start
    cov_fwd = np.zeros(region_len, dtype=np.int32)
    cov_rev = np.zeros(region_len, dtype=np.int32)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(chrom, region_start, region_end):
            if read.is_unmapped or read.mapping_quality < min_mapq:
                continue
            target = cov_rev if read.is_reverse else cov_fwd
            # get_blocks() returns the aligned segments (handles N/D in CIGAR correctly)
            for block_start, block_end in read.get_blocks():
                # Clip block to the region we're processing
                s = max(block_start, region_start) - region_start
                e = min(block_end, region_end) - region_start
                if e > s:
                    target[s:e] += 1

    return chrom, region_start, cov_fwd, cov_rev


def _build_work_list(bam_path, min_mapq):
    """Open BAM header and build a list of (bam, chrom, start, end, min_mapq) work units.

    Large chromosomes are split into chunks of CHUNK_SIZE bp. Work is sorted by
    chunk size descending so the heaviest chunks start first (reduces tail latency).
    """
    import pysam

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        refs = list(bam.header.references)
        lens = list(bam.header.lengths)

    chrom_order = {ref: i for i, ref in enumerate(refs)}

    work = []
    for chrom, length in zip(refs, lens):
        if length <= CHUNK_SIZE:
            work.append((bam_path, chrom, 0, length, min_mapq))
        else:
            for start in range(0, length, CHUNK_SIZE):
                end = min(start + CHUNK_SIZE, length)
                work.append((bam_path, chrom, start, end, min_mapq))

    # Largest chunks first → better load balancing
    work.sort(key=lambda w: -(w[3] - w[2]))
    return work, chrom_order, dict(zip(refs, lens))


def _scatter_coverage(bam_path, min_mapq, cpus):
    """Scatter coverage computation across chromosomes/chunks, gather into full arrays.

    Returns:
        cov_fwd, cov_rev : dict[str, np.ndarray]  (one array per chromosome, full length)
        chrom_order      : dict[str, int]         (preserves BAM header order)
    """
    work, chrom_order, chrom_lens = _build_work_list(bam_path, min_mapq)
    log.info(f"Scattering {len(work)} chunk(s) across {cpus} worker(s) "
             f"(chunk size: {CHUNK_SIZE:,} bp)")

    # Pre-allocate full-length arrays per chromosome; chunks fill in their slice
    cov_fwd = {c: np.zeros(L, dtype=np.int32) for c, L in chrom_lens.items()}
    cov_rev = {c: np.zeros(L, dtype=np.int32) for c, L in chrom_lens.items()}

    if cpus == 1:
        # Skip multiprocessing overhead entirely
        for w in work:
            chrom, start, fwd, rev = _coverage_chunk(w)
            cov_fwd[chrom][start:start + len(fwd)] = fwd
            cov_rev[chrom][start:start + len(rev)] = rev
    else:
        with Pool(cpus) as pool:
            for chrom, start, fwd, rev in pool.imap_unordered(_coverage_chunk, work):
                cov_fwd[chrom][start:start + len(fwd)] = fwd
                cov_rev[chrom][start:start + len(rev)] = rev

    # Drop chromosomes with zero coverage to save memory before write step
    cov_fwd = {c: a for c, a in cov_fwd.items() if a.any()}
    cov_rev = {c: a for c, a in cov_rev.items() if a.any()}

    return cov_fwd, cov_rev, chrom_order


def _emit_bedgraph_intervals(chrom, cov_array):
    """Run-length encode a coverage array into bedGraph intervals (vectorized).

    Yields tuples (chrom, start, end, value) for nonzero runs.
    """
    # Find boundaries where coverage changes value
    # Sentinels of -1 at start and end guarantee that any nonzero run is bounded.
    diff_pts = np.where(np.diff(cov_array, prepend=-1, append=-1) != 0)[0]
    # diff_pts gives indices where a new run starts; consecutive pairs define [start, end)
    starts = diff_pts[:-1]
    ends   = diff_pts[1:]
    values = cov_array[starts]

    # Emit only nonzero runs (bedGraph omits zero-coverage intervals)
    mask = values > 0
    for s, e, v in zip(starts[mask], ends[mask], values[mask]):
        yield chrom, int(s), int(e), int(v)


def _open_output(path):
    """Open output path, transparently gzip-compressing if path ends in .gz."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "wt")
    return open(path, "w")


def _write_bedgraph(coverage, chrom_order, output_path, strand_label):
    """Write a strand's coverage dict to a bedGraph file."""
    n_chroms = len(coverage)
    if n_chroms == 0:
        log.warning(f"No coverage found for strand {strand_label} — writing empty bedGraph")
        Path(output_path).touch()
        return 0

    sorted_chroms = sorted(coverage.keys(), key=lambda c: chrom_order.get(c, 999_999))
    written = 0
    with _open_output(output_path) as out:
        for chrom in sorted_chroms:
            for c, s, e, v in _emit_bedgraph_intervals(chrom, coverage[chrom]):
                out.write(f"{c}\t{s}\t{e}\t{v}\n")
                written += 1

    log.info(f"Written {written:,} bedGraph intervals to {output_path} (strand {strand_label})")
    return written


def _write_bedgraph_worker(coverage, chrom_order, output_path, strand_label):
    """Wrapper for Process — re-runs _write_bedgraph in a child process."""
    _write_bedgraph(coverage, chrom_order, output_path, strand_label)


def bam_to_bedgraph_both(bam_path, fwd_output, rev_output, min_mapq=0, cpus=1):
    """Single logical pass (parallel by chunk) producing both strand bedGraphs.

    Forward and reverse outputs are written in parallel.
    """
    log.info(f"Reading {bam_path} (both strands, min_mapq={min_mapq}, cpus={cpus})")
    cov_fwd, cov_rev, chrom_order = _scatter_coverage(bam_path, min_mapq, cpus)

    log.info(f"Coverage built: + strand on {len(cov_fwd)} chrom(s), "
             f"- strand on {len(cov_rev)} chrom(s)")

    # Write both strands in parallel (each pickles its own dict to a child process).
    # For very large genomes this may double peak memory briefly; if that's an issue,
    # call _write_bedgraph sequentially instead.
    if cpus >= 2:
        p_fwd = Process(target=_write_bedgraph_worker,
                        args=(cov_fwd, chrom_order, fwd_output, "+"))
        p_rev = Process(target=_write_bedgraph_worker,
                        args=(cov_rev, chrom_order, rev_output, "-"))
        p_fwd.start(); p_rev.start()
        p_fwd.join();  p_rev.join()
        if p_fwd.exitcode != 0 or p_rev.exitcode != 0:
            log.error(f"bedGraph writer failed (fwd exit={p_fwd.exitcode}, "
                      f"rev exit={p_rev.exitcode})")
            sys.exit(1)
    else:
        _write_bedgraph(cov_fwd, chrom_order, fwd_output, "+")
        _write_bedgraph(cov_rev, chrom_order, rev_output, "-")


def bam_to_bedgraph(bam_path, strand, output_path, min_mapq=0, cpus=1):
    """Single-strand mode (backward compatibility)."""
    log.info(f"Reading {bam_path} (strand={strand}, min_mapq={min_mapq}, cpus={cpus})")
    cov_fwd, cov_rev, chrom_order = _scatter_coverage(bam_path, min_mapq, cpus)
    coverage = cov_fwd if strand == "+" else cov_rev
    _write_bedgraph(coverage, chrom_order, output_path, strand)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bam", help="Input BAM (must be sorted+indexed)")
    parser.add_argument("--min-mapq", type=int, default=0,
                        help="Minimum mapping quality (default: 0)")
    parser.add_argument("--cpus", type=int, default=1,
                        help="Parallel workers (default: 1)")

    # Both-strands mode (preferred — single pass)
    parser.add_argument("--fwd-output", help="Output bedGraph for forward strand "
                                              "(.gz extension → gzip-compressed)")
    parser.add_argument("--rev-output", help="Output bedGraph for reverse strand "
                                              "(.gz extension → gzip-compressed)")

    # Single-strand mode (backward compat)
    parser.add_argument("--strand", choices=["+", "-"],
                        help="Strand to compute coverage for (single-strand mode)")
    parser.add_argument("--output", help="Output bedGraph path (single-strand mode)")

    args = parser.parse_args()

    if not Path(args.bam).exists():
        log.error(f"BAM file not found: {args.bam}")
        sys.exit(1)

    # Verify index exists (fetch() requires it)
    bai = Path(str(args.bam) + ".bai")
    csi = Path(str(args.bam) + ".csi")
    if not (bai.exists() or csi.exists()):
        log.error(f"BAM index not found ({bai} or {csi}). Run `samtools index` first.")
        sys.exit(1)

    if args.fwd_output and args.rev_output:
        bam_to_bedgraph_both(args.bam, args.fwd_output, args.rev_output,
                              args.min_mapq, args.cpus)
    elif args.strand and args.output:
        bam_to_bedgraph(args.bam, args.strand, args.output,
                        args.min_mapq, args.cpus)
    else:
        log.error("Provide either --fwd-output + --rev-output, or --strand + --output")
        sys.exit(1)


if __name__ == "__main__":
    main()