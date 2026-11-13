#!/usr/bin/env bash
# =============================================================================
# Build all scPolASeq Apptainer SIF images from conda environment specs.
#
# Usage:
#   cd /path/to/scPolASeq
#   bash containers/build.sh [--cache-dir /path/to/sif-cache]
#
# Options:
#   --cache-dir DIR   Output directory for .sif files (default: ${XDG_CACHE_HOME:-$HOME/.cache}/scpolaseq/apptainer)
#   --env NAME        Build only one env (python|alignment|clustering|r|scapture)
#   --force           Rebuild even if .sif already exists
#
# Requirements: Apptainer >= 1.0 with fakeroot enabled for this user.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "${SCRIPT_DIR}")"
CACHE_DIR="${SCPOLASEQ_APPTAINER_CACHE:-${XDG_CACHE_HOME:-${HOME}/.cache}/scpolaseq/apptainer}"
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

ENVS=(python alignment clustering r scapture)
if [[ -n "${ONLY_ENV}" ]]; then
    ENVS=("${ONLY_ENV}")
fi

for env in "${ENVS[@]}"; do
    def="${SCRIPT_DIR}/scpolaseq-${env}.def"
    sif="${CACHE_DIR}/scpolaseq-${env}.sif"

    if [[ ! -f "${def}" ]]; then
        echo "[ERROR] No .def file for env '${env}': ${def}" >&2
        exit 1
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

    apptainer exec "${sif}" true
    echo "[OK] Built: ${sif}"
done

(
    cd "${CACHE_DIR}"
    shopt -s nullglob
    images=(scpolaseq-*.sif)
    if ((${#images[@]})); then
        sha256sum "${images[@]}" > SHA256SUMS
    else
        : > SHA256SUMS
    fi
)

echo ""
echo "Done. SIF paths:"
for env in "${ENVS[@]}"; do
    sif="${CACHE_DIR}/scpolaseq-${env}.sif"
    [[ -f "${sif}" ]] && ls -lh "${sif}"
done
