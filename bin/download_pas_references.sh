#!/usr/bin/env bash
# download_pas_references.sh
#
# Downloads and normalises polyA site reference datasets for use with scPolASeq.
#
# Sources handled:
#   1. PolyASite v2.0  (GRCh38.96)          — 3'-seq bulk, already staged
#   2. PolyASite v3.0  (GRCh38.GENCODE_42)  — 10x single-cell 3'-end, 813 libs
#   3. PolyA_DB v4.1   Human                — manually curated PAS catalogue
#
# Each source is normalised to a tab-separated file with a header row containing
# at minimum the columns required by the pipeline's PAS reference builder:
#
#   gene_id  gene_name  chrom  start  end  strand  score  source
#
# Usage:
#   bash bin/download_pas_references.sh [--outdir <dir>] [--skip-existing]
#
# Defaults:
#   --outdir /scratch/data/polya_atlas
#
# The script also writes a summary report (<outdir>/pas_references_summary.tsv)
# and prints column headers + first 3 data rows of each downloaded file to
# help verify format compatibility before a pipeline run.
#
set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
OUTDIR="/scratch/data/polya_atlas"
SKIP_EXISTING=0

# ── Sources ───────────────────────────────────────────────────────────────────
POLYASITE2_URL="https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.bed.gz"
POLYASITE3_URL="https://polyasite.unibas.ch/download/atlas/3.0/GRCh38.GENCODE_42/atlas.clusters.3.0.GRCh38.GENCODE_42.bed.gz"
POLYADB41_URL="https://exon.apps.wistar.org/polya_db/v4/download/4.1/HumanPas.zip"

# ── Argument parsing ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --outdir)      OUTDIR="$2";      shift 2 ;;
        --skip-existing) SKIP_EXISTING=1; shift   ;;
        *) echo "Unknown option: $1" >&2; exit 1  ;;
    esac
done

mkdir -p "$OUTDIR"
SUMMARY="$OUTDIR/pas_references_summary.tsv"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
warn() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARN: $*" >&2; }

# ── Helper: normalise PolyASite BED to pipeline TSV ───────────────────────────
# PolyASite BED columns (with header in v2 tsv, headerless in raw .bed):
#   chrom  chromStart  chromEnd  name  score  strand  [rep  frac_samples
#   nr_prots  annotation  gene_name  gene_id  repSite_signals ...]
#
# Output columns: gene_id  gene_name  chrom  start  end  strand  score  source
normalise_polyasite() {
    local src="$1" out="$2" source_label="$3"
    log "Normalising PolyASite → $out"
    python3 - "$src" "$out" "$source_label" <<'PYEOF'
import csv, sys
from pathlib import Path

src, out, source_label = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]

BED_COLS = ["chrom","chromStart","chromEnd","name","score","strand",
            "rep","frac_samples","nr_prots","annotation","gene_name",
            "gene_id","repSite_signals"]

def fix_chrom(c):
    c = str(c).strip()
    if not c: return c
    if c in ("M","MT"): return "chrM"
    return c if c.startswith("chr") else f"chr{c}"

# Read ALL rows into memory while the file is open so no lazy iterators
# escape the with-block (that was the source of the closed-file error).
rows = []
with src.open(encoding="utf-8") as fh:
    first = fh.readline().rstrip("\n")
    if first.startswith("chrom\t") or first.startswith("#"):
        # File has a header row — let DictReader use it
        header = first.lstrip("#").split("\t")
        reader = csv.DictReader(fh, fieldnames=header, delimiter="\t")
        rows = list(reader)
    else:
        # Headerless BED — map positional columns
        parts = first.split("\t")
        rows.append(dict(zip(BED_COLS, parts)))
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts:
                rows.append(dict(zip(BED_COLS, parts)))

fieldnames_out = ["gene_id","gene_name","chrom","start","end","strand","score","source"]
written = 0
with out.open("w", newline="", encoding="utf-8") as fout:
    w = csv.DictWriter(fout, fieldnames=fieldnames_out, delimiter="\t",
                       extrasaction="ignore")
    w.writeheader()
    for row in rows:
        chrom = fix_chrom(row.get("chrom",""))
        start = row.get("start") or row.get("chromStart","0")
        end   = row.get("end")   or row.get("chromEnd", start)
        nan_vals = {"", "na", "nan", "none"}
        def clean_v(v, default=""):
            s = str(v).strip().lower() if v else ""
            return default if s in nan_vals else str(v).strip()
        gene_id   = clean_v(row.get("gene_id") or row.get("Gene"), "unknown_gene")
        gene_name = clean_v(row.get("gene_name") or row.get("Alias"))
        strand    = clean_v(row.get("strand"), "+")
        score     = clean_v(row.get("score"), "0")
        if not chrom:
            continue
        w.writerow({"gene_id": gene_id, "gene_name": gene_name,
                    "chrom": chrom, "start": start, "end": end,
                    "strand": strand, "score": score, "source": source_label})
        written += 1

print(f"  Wrote {written:,} rows to {out}")
PYEOF
}

# ── Helper: normalise PolyA_DB TSV to pipeline TSV ────────────────────────────
# PolyA_DB v4 columns (tab-separated, with header):
#   Gene  Alias  Species  ...  Chr  Strand  Position  ...
# 'Gene' = Entrez gene ID or symbol; 'Position' is a single-base coordinate.
normalise_polyadb() {
    local src="$1" out="$2"
    log "Normalising PolyA_DB → $out"
    python3 - "$src" "$out" <<'PYEOF'
import csv, sys
from pathlib import Path

src, out = Path(sys.argv[1]), Path(sys.argv[2])

def fix_chrom(c):
    c = str(c).strip()
    if not c: return c
    if c in ("M","MT"): return "chrM"
    return c if c.startswith("chr") else f"chr{c}"

def clean(v, default=""):
    """Return empty string for pandas-style NA/NaN placeholders."""
    v = str(v).strip() if v is not None else ""
    return default if v.lower() in ("", "na", "nan", "none") else v

fieldnames_out = ["gene_id","gene_name","chrom","start","end","strand","score","source"]

written = 0
with src.open(encoding="utf-8") as fh, out.open("w", newline="", encoding="utf-8") as fout:
    reader = csv.DictReader(fh, delimiter="\t")
    cols = reader.fieldnames or []
    print(f"  PolyA_DB columns (first 8): {cols[:8]}")
    w = csv.DictWriter(fout, fieldnames=fieldnames_out, delimiter="\t",
                       extrasaction="ignore")
    w.writeheader()

    # PolyA_DB v4.x: PAS_ID = "chr1:+:940182"  GeneSymbol = symbol
    # PolyA_DB v3.x: Chr / Strand / Position / Gene / Alias
    is_v4 = "PAS_ID" in cols

    for row in reader:
        if is_v4:
            pas_id = row.get("PAS_ID","").strip()   # chr1:+:940182
            parts = pas_id.split(":")
            if len(parts) < 3:
                continue
            chrom  = fix_chrom(parts[0])
            strand = parts[1] if parts[1] in ("+","-",".") else "+"
            pos    = parts[2]
            gene_name = clean(row.get("GeneSymbol",""))
            gene_id   = gene_name or "unknown_gene"
            score     = clean(row.get("PolyaStrength_percentile",""), "0")
        else:
            # v3.x style
            chrom     = fix_chrom(clean(row.get("Chr") or row.get("chrom","")))
            strand    = clean(row.get("Strand") or row.get("strand",""), "+")
            pos       = clean(row.get("Position") or row.get("position") or
                              row.get("start") or row.get("chromStart",""), "0")
            gene_id   = clean(row.get("Gene") or row.get("gene_id",""), "unknown_gene")
            gene_name = clean(row.get("Alias") or row.get("gene_name",""))
            score     = clean(row.get("Score") or row.get("score",""), "0")
        if not chrom:
            continue
        w.writerow({"gene_id": gene_id, "gene_name": gene_name,
                    "chrom": chrom, "start": pos, "end": pos,
                    "strand": strand, "score": score,
                    "source": "polyadb4.1"})
        written += 1

print(f"  Wrote {written:,} rows to {out}")
PYEOF
}

# ── Helper: show column summary ────────────────────────────────────────────────
show_summary() {
    local label="$1" file="$2"
    echo ""
    echo "┌── $label"
    echo "│   File : $file"
    if [[ ! -f "$file" ]]; then
        echo "│   STATUS: FILE NOT FOUND"
        echo "└──"
        return
    fi
    local rows
    rows=$(tail -n +2 "$file" | wc -l)
    echo "│   Rows : $rows"
    echo "│   Header:"
    head -1 "$file" | tr '\t' '\n' | nl -ba | sed 's/^/│     /'
    echo "│   First 3 data rows:"
    tail -n +2 "$file" | head -3 | cut -f1-8 | sed 's/^/│     /'
    echo "└──"
}

# ── Helper: safe download with progress ───────────────────────────────────────
download() {
    local url="$1" dest="$2"
    if [[ $SKIP_EXISTING -eq 1 && -f "$dest" ]]; then
        log "Skipping (exists): $dest"
        return 0
    fi
    log "Downloading: $url"
    curl -fSL --retry 3 --retry-delay 5 -o "$dest" "$url"
}

# ═══════════════════════════════════════════════════════════════════════════════
# 1. PolyASite v2.0  (already staged — just normalise if needed)
# ═══════════════════════════════════════════════════════════════════════════════
PA2_RAW="$OUTDIR/polyasite2_GRCh38.tsv"
PA2_NORM="$OUTDIR/polyasite2.0_GRCh38_norm.tsv"

if [[ -f "$PA2_RAW" ]]; then
    log "PolyASite v2.0 raw file already present: $PA2_RAW"
else
    log "PolyASite v2.0 not found — downloading..."
    PA2_GZ="$OUTDIR/polyasite2_GRCh38.bed.gz"
    download "$POLYASITE2_URL" "$PA2_GZ"
    log "Decompressing..."
    gunzip -c "$PA2_GZ" > "$PA2_RAW"
fi

if [[ ! -f "$PA2_NORM" || $SKIP_EXISTING -eq 0 ]]; then
    normalise_polyasite "$PA2_RAW" "$PA2_NORM" "polyasite2.0"
fi

# ═══════════════════════════════════════════════════════════════════════════════
# 2. PolyASite v3.0  (GRCh38.GENCODE_42, 10x single-cell)
# ═══════════════════════════════════════════════════════════════════════════════
PA3_GZ="$OUTDIR/polyasite3.0_GRCh38_10x.bed.gz"
PA3_RAW="$OUTDIR/polyasite3.0_GRCh38_10x.bed"
PA3_NORM="$OUTDIR/polyasite3.0_GRCh38_10x_norm.tsv"

if [[ ! -f "$PA3_GZ" || $SKIP_EXISTING -eq 0 ]]; then
    download "$POLYASITE3_URL" "$PA3_GZ"
fi

if [[ ! -f "$PA3_RAW" || $SKIP_EXISTING -eq 0 ]]; then
    log "Decompressing PolyASite v3.0..."
    gunzip -kf "$PA3_GZ"
    mv "${PA3_GZ%.gz}" "$PA3_RAW" 2>/dev/null || true
fi

if [[ ! -f "$PA3_NORM" || $SKIP_EXISTING -eq 0 ]]; then
    normalise_polyasite "$PA3_RAW" "$PA3_NORM" "polyasite3.0_10x"
fi

# ═══════════════════════════════════════════════════════════════════════════════
# 3. PolyA_DB v4.1 Human
# ═══════════════════════════════════════════════════════════════════════════════
PDB_ZIP="$OUTDIR/polyadb4.1_human.zip"
PDB_RAW_DIR="$OUTDIR/polyadb4.1_human_raw"
PDB_NORM="$OUTDIR/polyadb4.1_human_norm.tsv"

if [[ ! -f "$PDB_ZIP" || $SKIP_EXISTING -eq 0 ]]; then
    download "$POLYADB41_URL" "$PDB_ZIP"
fi

if [[ ! -d "$PDB_RAW_DIR" || $SKIP_EXISTING -eq 0 ]]; then
    log "Unzipping PolyA_DB v4.1..."
    mkdir -p "$PDB_RAW_DIR"
    python3 -c "import zipfile, sys; zipfile.ZipFile(sys.argv[1]).extractall(sys.argv[2])" \
        "$PDB_ZIP" "$PDB_RAW_DIR"
fi

# Find the main TSV/txt file inside the zip
PDB_TSV=$(find "$PDB_RAW_DIR" -maxdepth 2 \( -name "*.tsv" -o -name "*.txt" \) \
    ! -name "README*" ! -name "*.md" | sort | head -1)

if [[ -z "$PDB_TSV" ]]; then
    warn "Could not locate TSV/TXT inside PolyA_DB zip. Contents of $PDB_RAW_DIR:"
    find "$PDB_RAW_DIR" -type f | sed 's/^/  /'
else
    log "PolyA_DB raw file: $PDB_TSV"
    if [[ ! -f "$PDB_NORM" || $SKIP_EXISTING -eq 0 ]]; then
        normalise_polyadb "$PDB_TSV" "$PDB_NORM"
    fi
fi

# ═══════════════════════════════════════════════════════════════════════════════
# 4. Format summary report
# ═══════════════════════════════════════════════════════════════════════════════
log "Writing summary report..."

show_summary "PolyASite v2.0  (GRCh38.96, bulk 3'-seq)"  "$PA2_NORM"
show_summary "PolyASite v3.0  (GRCh38.GENCODE_42, 10x)"  "$PA3_NORM"
[[ -f "$PDB_NORM" ]] && show_summary "PolyA_DB v4.1   (Human curated)" "$PDB_NORM"

# Write machine-readable summary
{
    printf "source\tfile\trows\tnotes\n"
    for f_label_file in \
        "polyasite2.0|$PA2_NORM|bulk 3'-seq 221 libs" \
        "polyasite3.0_10x|$PA3_NORM|10x 3'-end 813 libs" \
        "polyadb4.1|$PDB_NORM|manually curated"
    do
        IFS="|" read -r lbl f notes <<< "$f_label_file"
        rows=0; [[ -f "$f" ]] && rows=$(tail -n +2 "$f" | wc -l)
        printf "%s\t%s\t%s\t%s\n" "$lbl" "$f" "$rows" "$notes"
    done
} > "$SUMMARY"

log "Summary written to: $SUMMARY"
echo ""
log "Done. Recommended --known_polya values:"
echo "  Single-cell optimised (recommended for 10x):  $PA3_NORM"
echo "  Bulk, broad coverage:                         $PA2_NORM"
[[ -f "$PDB_NORM" ]] && \
echo "  Curated (highest precision, lower recall):    $PDB_NORM"
echo ""
echo "  Update params.known_polya in nextflow.config or test_pbmc1k_deepthought.config"
echo "  to point to the desired file."
