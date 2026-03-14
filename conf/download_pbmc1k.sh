#!/usr/bin/env bash
# =============================================================================
# download_pbmc1k.sh — Download and stage PBMC 1k v3 dataset from 10x Genomics
#
# Usage:
#   bash conf/download_pbmc1k.sh [DEST_DIR]
#
# Default DEST_DIR: /cluster/datasets/pbmc1k
#
# Also downloads the 10x v3 barcode whitelist if not already present.
# =============================================================================
set -euo pipefail

DEST_DIR="${1:-/cluster/datasets/pbmc1k}"
WHITELIST_DIR="/cluster/datasets/10x_whitelists"

FASTQ_URL="https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar"
WHITELIST_URL="https://cf.10xgenomics.com/supp/cell-exp/whitelists/3M-february-2018.txt.gz"

echo "=== PBMC 1k v3 downloader ==="
echo "Destination : ${DEST_DIR}"
echo "Whitelist   : ${WHITELIST_DIR}"
echo

# ── FASTQs ────────────────────────────────────────────────────────────────────
mkdir -p "${DEST_DIR}/fastqs"

TAR_FILE="${DEST_DIR}/pbmc_1k_v3_fastqs.tar"

if ls "${DEST_DIR}/fastqs/"*.fastq.gz >/dev/null 2>&1; then
    echo "[SKIP] FASTQs already present in ${DEST_DIR}/fastqs/"
else
    echo "[1/3] Downloading PBMC 1k v3 FASTQs (~5 GB)..."
    curl --progress-bar -L "${FASTQ_URL}" -o "${TAR_FILE}"

    echo "[2/3] Extracting FASTQs..."
    tar -xf "${TAR_FILE}" -C "${DEST_DIR}/fastqs" --strip-components=1

    rm -f "${TAR_FILE}"
    echo "      Extracted:"
    ls "${DEST_DIR}/fastqs/"
fi

# ── Barcode whitelist ─────────────────────────────────────────────────────────
mkdir -p "${WHITELIST_DIR}"

WHITELIST_FILE="${WHITELIST_DIR}/3M-february-2018.txt"
if [[ -f "${WHITELIST_FILE}" ]]; then
    echo "[SKIP] Whitelist already present: ${WHITELIST_FILE}"
else
    echo "[3/3] Downloading 10x v3 barcode whitelist..."
    curl --progress-bar -L "${WHITELIST_URL}" -o "${WHITELIST_FILE}.gz"
    gunzip "${WHITELIST_FILE}.gz"
    echo "      Whitelist saved to ${WHITELIST_FILE}"
fi

echo
echo "=== Download complete ==="
echo
echo "Next steps:"
echo "  1. Build or download a GRCh38 STAR index:"
echo "     mkdir -p /cluster/datasets/reference/GRCh38"
echo "     # Download from GENCODE or use REFERENCE_PREPARE to build it"
echo
echo "  2. Launch the PBMC 1k test:"
echo "     nextflow run main.nf \\"
echo "       -profile test_pbmc1k,singularity \\"
echo "       -c /cluster/projects/nextflow.config \\"
echo "       --outdir /cluster/projects/scPolASeq_pbmc1k/results"
