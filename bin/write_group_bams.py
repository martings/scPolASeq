#!/usr/bin/env python3

import argparse
import csv
import logging
import re
from pathlib import Path


def sanitize(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value)


def grouped_bam_name(library_id: str | None, group_id: str) -> str:
    safe_group_id = sanitize(group_id)
    if library_id:
        return f"{sanitize(library_id)}.{safe_group_id}.grouped.bam"
    return f"{safe_group_id}.grouped.bam"


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    log = logging.getLogger()

    parser = argparse.ArgumentParser()
    parser.add_argument("--bam",           required=False, default=None,
                        help="Filtered BAM to split by group (optional; stub if absent/empty)")
    parser.add_argument("--group-map",     required=True)
    parser.add_argument("--sample-id",     required=True)
    parser.add_argument("--library-id",    required=False, default=None)
    parser.add_argument("--group-level",   required=False, default=None,
                        help="Only emit groups for this group_level")
    parser.add_argument("--out-dir",       required=True)
    parser.add_argument("--out-manifest",  required=True)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Collect group_ids for this sample/level from the group_map
    groups = []
    with Path(args.group_map).open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if args.sample_id and row.get("sample_id") != args.sample_id:
                continue
            row_library_id = (row.get("library_id") or "").strip()
            if args.library_id and row_library_id and row_library_id != args.library_id:
                continue
            level = row.get("group_level", "")
            if args.group_level and level != args.group_level:
                continue
            group_id = row.get("group_id", "")
            if level and group_id and (level, group_id) not in groups:
                groups.append((level, group_id))

    # Attempt real BAM splitting; fall back to placeholder output
    bam_valid = False
    if args.bam:
        bam_path = Path(args.bam)
        if bam_path.exists() and bam_path.stat().st_size > 0:
            try:
                import pysam
                _split_bam(args.bam, groups, out_dir, args.library_id, log)
                bam_valid = True
            except (ImportError, ValueError, OSError) as exc:
                log.warning(f"BAM split unavailable ({exc}); writing placeholder BAMs")

    if not bam_valid:
        for _level, group_id in groups:
            placeholder = out_dir / grouped_bam_name(args.library_id, group_id)
            placeholder.write_bytes(b"")
            log.info(f"Placeholder: {placeholder.name}")

    manifest_rows = []
    for level, group_id in groups:
        bam_name = grouped_bam_name(args.library_id, group_id)
        n_reads = 0
        manifest_rows.append({
            "group_level": level,
            "group_id": group_id,
            "bam_file": bam_name,
            "n_reads": n_reads,
        })

    with Path(args.out_manifest).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["group_level", "group_id", "bam_file", "n_reads"],
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    log.info(f"Groups written: {len(manifest_rows)}")


def _split_bam(bam_path, groups, out_dir, library_id, log):
    import pysam

    handles = {}
    tag_for_group = {
        "cluster": "CL",
        "cell_type": "CT",
        "cell": "CB",
    }
    with pysam.AlignmentFile(bam_path, "rb") as bam_in:
        header = bam_in.header
        for level, group_id in groups:
            out_bam = out_dir / grouped_bam_name(library_id, group_id)
            handles[(level, group_id)] = pysam.AlignmentFile(str(out_bam), "wb", header=header)

        for read in bam_in.fetch(until_eof=True):
            for level, group_id in groups:
                tag_name = tag_for_group.get(level, "CL")
                try:
                    group_tag = read.get_tag(tag_name)
                except KeyError:
                    continue
                if group_tag == group_id:
                    handles[(level, group_id)].write(read)

    for h in handles.values():
        h.close()


if __name__ == "__main__":
    main()
