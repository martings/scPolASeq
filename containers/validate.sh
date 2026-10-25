#!/usr/bin/env bash
# Smoke-check the five runtime SIFs before using them in Nextflow.
set -euo pipefail

SIF_DIR="${1:-${SCPOLASEQ_APPTAINER_CACHE:-${XDG_CACHE_HOME:-${HOME}/.cache}/scpolaseq/apptainer}}"
for name in python alignment clustering r scapture; do
    image="${SIF_DIR}/scpolaseq-${name}.sif"
    [[ -s "${image}" ]] || { echo "Missing SIF: ${image}" >&2; exit 1; }
    apptainer exec "${image}" true
done

apptainer exec "${SIF_DIR}/scpolaseq-python.sif" python -c \
    'import numpy, pandas, pysam; print("python", numpy.__version__, pandas.__version__, pysam.__version__)'
apptainer exec "${SIF_DIR}/scpolaseq-python.sif" bash -lc \
    'command -v bedGraphToBigWig >/dev/null'
apptainer exec "${SIF_DIR}/scpolaseq-alignment.sif" bash -lc \
    'STAR --version; samtools --version | head -1; python -c "import pysam; print(pysam.__version__)"'
apptainer exec "${SIF_DIR}/scpolaseq-clustering.sif" python -c \
    'import scanpy, anndata; print("scanpy", scanpy.__version__, "anndata", anndata.__version__)'
apptainer exec "${SIF_DIR}/scpolaseq-r.sif" Rscript -e \
    'library(Sierra); cat("Sierra", as.character(packageVersion("Sierra")), "OK\n")'
apptainer exec "${SIF_DIR}/scpolaseq-scapture.sif" python -c \
    'import tensorflow as tf, numpy, pandas, h5py; print("tensorflow", tf.__version__, "numpy", numpy.__version__, "pandas", pandas.__version__, "h5py", h5py.__version__)'

echo "All Apptainer image smoke checks passed: ${SIF_DIR}"
