#!/usr/bin/env bash
# Smoke-check runtime SIFs before using them in Nextflow.
set -euo pipefail

SIF_DIR="${SCPOLASEQ_APPTAINER_CACHE:-${XDG_CACHE_HOME:-${HOME}/.cache}/scpolaseq/apptainer}"
ONLY_ENV=all
if [[ $# -gt 0 && "$1" != --* ]]; then
    SIF_DIR="$1"
    shift
fi
while [[ $# -gt 0 ]]; do
    case "$1" in
        --env) [[ $# -ge 2 ]] || { echo 'Missing --env value' >&2; exit 2; }
               ONLY_ENV="$2"; shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done
case "$ONLY_ENV" in
    all) ENVS=(python alignment clustering r scapture) ;;
    python|alignment|clustering|r|scapture) ENVS=("$ONLY_ENV") ;;
    *) echo "Unknown environment: $ONLY_ENV" >&2; exit 2 ;;
esac

for name in "${ENVS[@]}"; do
    image="${SIF_DIR}/scpolaseq-${name}.sif"
    [[ -s "${image}" ]] || { echo "Missing SIF: ${image}" >&2; exit 1; }
    apptainer exec "${image}" true
done

for name in "${ENVS[@]}"; do
    case "$name" in
        python)
            apptainer exec "${SIF_DIR}/scpolaseq-python.sif" python -c \
                'import numpy, pandas, pysam; print("python", numpy.__version__, pandas.__version__, pysam.__version__)'
            fixture_dir="$(mktemp -d "${TMPDIR:-/tmp}/scpolaseq-bigwig-sif.XXXXXX")"
            trap 'rm -rf -- "$fixture_dir"' EXIT
            printf 'chr1\t1000\n' > "${fixture_dir}/chrom.sizes"
            printf 'chr1\t100\t200\t3\n' > "${fixture_dir}/coverage.bedGraph"
            apptainer exec --cleanenv --bind "${fixture_dir}:${fixture_dir}" \
                "${SIF_DIR}/scpolaseq-python.sif" bedGraphToBigWig \
                "${fixture_dir}/coverage.bedGraph" "${fixture_dir}/chrom.sizes" \
                "${fixture_dir}/coverage.bw"
            [[ -s "${fixture_dir}/coverage.bw" ]] || {
                echo 'Python SIF produced an empty BigWig from synthetic coverage' >&2; exit 1;
            }
            magic="$(od -An -tx1 -N4 "${fixture_dir}/coverage.bw" | tr -d ' \n')"
            [[ "$magic" == 26fc8f88 ]] || {
                echo 'Python SIF produced an invalid BigWig header' >&2; exit 1;
            }
            rm -rf -- "$fixture_dir"
            trap - EXIT
            ;;
        alignment)
            apptainer exec "${SIF_DIR}/scpolaseq-alignment.sif" bash -lc \
                'STAR --version; samtools --version | head -1; python -c "import pysam; print(pysam.__version__)"'
            ;;
        clustering)
            apptainer exec "${SIF_DIR}/scpolaseq-clustering.sif" python -c \
                'import scanpy, anndata; print("scanpy", scanpy.__version__, "anndata", anndata.__version__)'
            ;;
        r)
            apptainer exec "${SIF_DIR}/scpolaseq-r.sif" Rscript -e \
                'library(Sierra); cat("Sierra", as.character(packageVersion("Sierra")), "OK\n")'
            ;;
        scapture)
            apptainer exec "${SIF_DIR}/scpolaseq-scapture.sif" python -c \
                'import tensorflow as tf, numpy, pandas, h5py; print("tensorflow", tf.__version__, "numpy", numpy.__version__, "pandas", pandas.__version__, "h5py", h5py.__version__)'
            ;;
    esac
done

echo "All Apptainer image smoke checks passed: ${SIF_DIR}"
