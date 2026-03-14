#!/usr/bin/env bash
# =============================================================================
# setup_human_references.sh
# One-shot script to stage all human (GRCh38) reference data required by
# the scPolASeq PBMC 1k test on the QNAP NAS (/cluster).
#
# Usage:
#   bash conf/setup_human_references.sh [--skip-star-index]
#
# What it downloads / builds:
#   1. GRCh38 primary assembly genome     (GENCODE release 44, ~900 MB gz)
#   2. GENCODE v44 annotation GTF         (~50 MB gz)
#   3. PolyASite 2.0 human polyA atlas    (GRCh38 / Ensembl 96, ~30 MB gz)
#   4. 10x Genomics v3 barcode whitelist  (3M-february-2018, ~58 MB gz)
#   5. PBMC 1k v3 FASTQs                  (~5 GB tar)
#   6. STAR 2.7 genome index              (needs STAR ≥ 2.7, ~30 min, ~27 GB)
#
# Directory layout on /cluster after completion:
#   /cluster/datasets/reference/GRCh38/
#       genome.fa            (decompressed, ~3.1 GB)
#       genes.gtf            (decompressed, ~1.6 GB)
#       chrom.sizes
#       star_index/          (built from genome.fa + genes.gtf)
#   /cluster/datasets/polya_atlas/
#       polyasite2_GRCh38.tsv
#   /cluster/datasets/10x_whitelists/
#       3M-february-2018.txt
#   /cluster/datasets/pbmc1k/fastqs/
#       pbmc_1k_v3_S1_L00{1,2}_R{1,2}_001.fastq.gz
# =============================================================================
set -euo pipefail

SKIP_STAR=false
for arg in "$@"; do [[ "$arg" == "--skip-star-index" ]] && SKIP_STAR=true; done

REF_DIR="/cluster/datasets/reference/GRCh38"
POLYA_DIR="/cluster/datasets/polya_atlas"
WL_DIR="/cluster/datasets/10x_whitelists"
PBMC_DIR="/cluster/datasets/pbmc1k/fastqs"
NCPUS="${NCPUS:-8}"

# URLs
GENOME_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz"
POLYA_URL="https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.tsv.gz"
WL_URL="https://cf.10xgenomics.com/supp/cell-exp/whitelists/3M-february-2018.txt.gz"
PBMC_URL="https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar"

log() { echo -e "\033[1;32m[$(date '+%H:%M:%S')]\033[0m $*"; }
skip() { echo -e "\033[1;33m[SKIP]\033[0m $*"; }
err() { echo -e "\033[1;31m[ERROR]\033[0m $*" >&2; exit 1; }

# ── 0. Connectivity check ─────────────────────────────────────────────────────
log "Checking /cluster mount..."
[[ -d /cluster ]] || err "/cluster not mounted. Run: sudo mount -t nfs 192.168.0.224:/cluster /cluster"
[[ -w /cluster ]] || err "/cluster is not writable. Check NFS permissions."
log "/cluster is reachable and writable."

# ── 1. Genome (GRCh38 primary assembly) ───────────────────────────────────────
mkdir -p "${REF_DIR}"
if [[ -s "${REF_DIR}/genome.fa" ]]; then
    skip "Genome already present: ${REF_DIR}/genome.fa"
else
    log "Downloading GRCh38 primary assembly genome (~900 MB)..."
    curl --progress-bar -L "${GENOME_URL}" | pigz -dc > "${REF_DIR}/genome.fa" \
        || { gzip -dc < <(curl -sL "${GENOME_URL}") > "${REF_DIR}/genome.fa"; }
    log "Genome saved: ${REF_DIR}/genome.fa  ($(du -sh "${REF_DIR}/genome.fa" | cut -f1))"
fi

# ── 2. GTF annotation (GENCODE v44) ───────────────────────────────────────────
if [[ -s "${REF_DIR}/genes.gtf" ]]; then
    skip "GTF already present: ${REF_DIR}/genes.gtf"
else
    log "Downloading GENCODE v44 GTF (~50 MB)..."
    curl --progress-bar -L "${GTF_URL}" | pigz -dc > "${REF_DIR}/genes.gtf" \
        || { gzip -dc < <(curl -sL "${GTF_URL}") > "${REF_DIR}/genes.gtf"; }
    log "GTF saved: ${REF_DIR}/genes.gtf"
fi

# ── 3. Chromosome sizes ────────────────────────────────────────────────────────
if [[ -s "${REF_DIR}/chrom.sizes" ]]; then
    skip "chrom.sizes already present"
else
    log "Generating chrom.sizes from genome FASTA..."
    if command -v samtools >/dev/null 2>&1; then
        samtools faidx "${REF_DIR}/genome.fa"
        cut -f1,2 "${REF_DIR}/genome.fa.fai" > "${REF_DIR}/chrom.sizes"
    elif command -v faSize >/dev/null 2>&1; then
        faSize -detailed "${REF_DIR}/genome.fa" > "${REF_DIR}/chrom.sizes"
    else
        # Pure awk fallback from indexed genome
        awk '/^>/{name=substr($1,2); next} {len[name]+=length($0)} END{for(n in len) print n"\t"len[n]}' \
            "${REF_DIR}/genome.fa" | sort -k1,1V > "${REF_DIR}/chrom.sizes"
    fi
    log "chrom.sizes: $(wc -l < "${REF_DIR}/chrom.sizes") chromosomes/scaffolds"
fi

# ── 4. PolyASite 2.0 atlas ────────────────────────────────────────────────────
mkdir -p "${POLYA_DIR}"
if [[ -s "${POLYA_DIR}/polyasite2_GRCh38.tsv" ]]; then
    skip "PolyASite atlas already present"
else
    log "Downloading PolyASite 2.0 human atlas (~30 MB)..."
    curl --progress-bar -L "${POLYA_URL}" | pigz -dc > "${POLYA_DIR}/polyasite2_GRCh38.tsv" \
        || { gzip -dc < <(curl -sL "${POLYA_URL}") > "${POLYA_DIR}/polyasite2_GRCh38.tsv"; }
    log "Atlas saved: ${POLYA_DIR}/polyasite2_GRCh38.tsv"
fi

# ── 5. 10x v3 barcode whitelist ───────────────────────────────────────────────
mkdir -p "${WL_DIR}"
if [[ -f "${WL_DIR}/3M-february-2018.txt" ]]; then
    skip "Barcode whitelist already present"
else
    log "Downloading 10x v3 barcode whitelist..."
    curl --progress-bar -L "${WL_URL}" \
        | pigz -dc > "${WL_DIR}/3M-february-2018.txt" \
        || { gzip -dc < <(curl -sL "${WL_URL}") > "${WL_DIR}/3M-february-2018.txt"; }
    log "Whitelist: $(wc -l < "${WL_DIR}/3M-february-2018.txt") barcodes"
fi

# ── 6. PBMC 1k v3 FASTQs ─────────────────────────────────────────────────────
mkdir -p "${PBMC_DIR}"
if ls "${PBMC_DIR}"/*.fastq.gz >/dev/null 2>&1; then
    skip "FASTQs already present in ${PBMC_DIR}"
else
    log "Downloading PBMC 1k v3 FASTQs (~5 GB)..."
    TAR_TMP="${PBMC_DIR}/../pbmc_1k_v3_fastqs.tar"
    curl --progress-bar -L "${PBMC_URL}" -o "${TAR_TMP}"
    log "Extracting FASTQs..."
    tar -xf "${TAR_TMP}" -C "${PBMC_DIR}" --strip-components=1
    rm -f "${TAR_TMP}"
    log "FASTQs: $(ls "${PBMC_DIR}"/*.fastq.gz | wc -l) files"
fi

# ── 7. STAR genome index ──────────────────────────────────────────────────────
STAR_INDEX="${REF_DIR}/star_index"
if [[ "${SKIP_STAR}" == "true" ]]; then
    skip "STAR index build skipped (--skip-star-index)"
elif [[ -f "${STAR_INDEX}/SA" ]]; then
    skip "STAR index already present: ${STAR_INDEX}/"
elif ! command -v STAR >/dev/null 2>&1; then
    echo
    echo "  STAR not found in PATH. Install STAR and re-run, or run via Apptainer:"
    echo "  apptainer exec /cluster/containers/star_2.7.sif STAR --runMode genomeGenerate \\"
    echo "    --genomeDir ${STAR_INDEX} --genomeFastaFiles ${REF_DIR}/genome.fa \\"
    echo "    --sjdbGTFfile ${REF_DIR}/genes.gtf --runThreadN ${NCPUS}"
else
    log "Building STAR index (~30 min, ${NCPUS} CPUs)..."
    mkdir -p "${STAR_INDEX}"
    STAR \
        --runMode         genomeGenerate \
        --runThreadN      "${NCPUS}" \
        --genomeDir       "${STAR_INDEX}" \
        --genomeFastaFiles "${REF_DIR}/genome.fa" \
        --sjdbGTFfile     "${REF_DIR}/genes.gtf" \
        --outFileNamePrefix "${STAR_INDEX}/"
    log "STAR index complete: ${STAR_INDEX}/"
fi

# ── Summary ───────────────────────────────────────────────────────────────────
echo
log "=== Reference setup complete ==="
echo
echo "  Genome   : ${REF_DIR}/genome.fa"
echo "  GTF      : ${REF_DIR}/genes.gtf"
echo "  Sizes    : ${REF_DIR}/chrom.sizes"
echo "  STAR idx : ${STAR_INDEX}/"
echo "  polyA    : ${POLYA_DIR}/polyasite2_GRCh38.tsv"
echo "  whitelist: ${WL_DIR}/3M-february-2018.txt"
echo "  FASTQs   : ${PBMC_DIR}/"
echo
echo "Launch the pipeline:"
echo "  cd $(pwd)"
echo "  nextflow run main.nf \\"
echo "    -profile test_pbmc1k,local_qnap \\"
echo "    -c /cluster/projects/nextflow.config \\"
echo "    --outdir /cluster/projects/scPolASeq_pbmc1k/results"
