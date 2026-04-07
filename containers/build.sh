#!/usr/bin/env bash
# =============================================================================
# Build all scPolASeq Apptainer SIF images from conda environment specs.
#
# Usage:
#   cd /scratch/work/scPolASeq
#   bash containers/build.sh [--cache-dir /scratch/cache/apptainer/scpolaseq]
#
# Options:
#   --cache-dir DIR   Output directory for .sif files (default: /scratch/cache/apptainer/scpolaseq)
#   --env NAME        Build only one env (python|alignment|clustering|r)
#   --force           Rebuild even if .sif already exists
#
# Requirements: apptainer >= 1.0, mamba accessible in base conda env
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "${SCRIPT_DIR}")"
CACHE_DIR="/scratch/cache/apptainer/scpolaseq"
FORCE=0
ONLY_ENV=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cache-dir) CACHE_DIR="$2"; shift 2 ;;
        --env)       ONLY_ENV="$2";  shift 2 ;;
        --force)     FORCE=1;        shift 1 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

mkdir -p "${CACHE_DIR}"

ENVS=(python alignment clustering r)
if [[ -n "${ONLY_ENV}" ]]; then
    ENVS=("${ONLY_ENV}")
fi

for env in "${ENVS[@]}"; do
    def="${SCRIPT_DIR}/scpolaseq-${env}.def"
    sif="${CACHE_DIR}/scpolaseq-${env}.sif"

    if [[ ! -f "${def}" ]]; then
        echo "[SKIP] No .def file for env '${env}': ${def}"
        continue
    fi

    if [[ -f "${sif}" && "${FORCE}" -eq 0 ]]; then
        echo "[SKIP] SIF already exists: ${sif}  (use --force to rebuild)"
        continue
    fi

    echo "======================================================================"
    echo "[BUILD] ${env}  →  ${sif}"
    echo "======================================================================"

    # Build from repo root so %files paths resolve relative to project
    (cd "${REPO_ROOT}" && \
        apptainer build --fakeroot "${sif}" "${def}")

    echo "[OK] Built: ${sif}"
done

echo ""
echo "Done. SIF paths:"
for env in "${ENVS[@]}"; do
    sif="${CACHE_DIR}/scpolaseq-${env}.sif"
    [[ -f "${sif}" ]] && ls -lh "${sif}"
done
