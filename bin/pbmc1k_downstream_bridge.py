#!/usr/bin/env python3
"""Run the unpublished downstream package on one scPolASeq PBMC1k result tree.

Only the first generated notebook cell is replaced. All analytical stages,
plots and enrichment continue to call scpolaseq_downstream itself.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sys


BASELINE = {
    "n_cells_qc_pass": 1113,
    "n_celltypes_curated": 16,
    "n_si_profiles": 22130,
    "n_differential_tests": 52750,
    "n_robust_apa_events": 942,
    "n_unique_robust_genes": 614,
}


def require_inputs(base: Path, package_src: Path) -> dict[str, str]:
    sample = "pbmc_1k_v3"
    required = {
        "clustered_h5ad": base / "labels" / f"{sample}.clustered.h5ad",
        "matrix": base / "alignment" / f"{sample}.Solo.out" / "Gene" / "filtered" / "matrix.mtx",
        "site_catalog": base / "apa" / "site_catalog.tsv",
        "apa_events": base / "apa_calls" / "apa_events.tsv",
        "pdui_matrix": base / "apa_calls" / "pdui_usage_matrix.tsv",
        "package_config": package_src / "scpolaseq_downstream" / "config.py",
        "package_template": package_src / "scpolaseq_downstream" / "templates.py",
    }
    missing = [name for name, path in required.items() if not path.is_file() or path.stat().st_size == 0]
    sierra = base / "sierra_quant" / "cluster"
    if not sierra.is_dir() or not any(p.stat().st_size > 0 for p in sierra.glob("*.sierra_quant.tsv")):
        missing.append("sierra_quant/cluster")
    if missing:
        raise ValueError(f"Missing or empty downstream inputs: {', '.join(missing)}")
    return {name: str(path.resolve()) for name, path in required.items()}


def source_fingerprint(package_src: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted((package_src / "scpolaseq_downstream").glob("*.py")):
        digest.update(path.name.encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def first_cell(package_src: Path, upstream: Path, work_root: Path) -> str:
    # Keep imports identical to the package-generated notebook while making
    # source and result locations explicit for this isolated invocation.
    return f'''import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

PACKAGE_SRC = Path({str(package_src)!r})
sys.path.insert(0, str(PACKAGE_SRC))
from scpolaseq_downstream.config import PipelineConfig
from scpolaseq_downstream.qc import run_qc_stage
from scpolaseq_downstream.celltyping import run_celltyping_stage
from scpolaseq_downstream.apa import run_apa_stage
from scpolaseq_downstream.enrichment import run_enrichment_stage
from scpolaseq_downstream.plotting import (
    plot_celltype_landscape, plot_marker_validation_heatmap,
    plot_apa_global_overview, plot_volcano_apa,
    plot_magnetization_concordance, plot_multilineage_gene_architecture,
    plot_reactome_and_targets,
)

config = PipelineConfig.for_pbmc1k(
    base_dir={str(upstream)!r},
    workspace_dir={str(work_root)!r},
)
subdirs = config.ensure_dirs()
print("Loaded PipelineConfig for:", config.dataset_name)
print("Output directory:", config.output_dir)
'''


def kernel_for_current_python(work_root: Path, package_src: Path) -> None:
    kernel_dir = work_root / ".jupyter" / "kernels" / "scpolaseq_bridge"
    kernel_dir.mkdir(parents=True)
    kernel = {
        "argv": [sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"],
        "display_name": "scPolASeq downstream bridge",
        "language": "python",
        "env": {"PYTHONPATH": os.pathsep.join([str(package_src), os.environ.get("PYTHONPATH", "")])},
    }
    (kernel_dir / "kernel.json").write_text(json.dumps(kernel), encoding="utf-8")
    os.environ["JUPYTER_PATH"] = os.pathsep.join(
        [str(work_root / ".jupyter"), os.environ.get("JUPYTER_PATH", "")]
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--upstream-results", type=Path, required=True)
    parser.add_argument("--package-src", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--generate-only", action="store_true")
    args = parser.parse_args()

    upstream = args.upstream_results.resolve()
    package_src = args.package_src.resolve()
    output = args.output_dir.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"Refusing to reuse a nonempty downstream output directory: {output}")
    inputs = require_inputs(upstream, package_src)
    output.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(package_src))
    import nbformat
    from scpolaseq_downstream.config import PipelineConfig
    from scpolaseq_downstream.templates import generate_template_notebook, execute_and_export_notebook

    config = PipelineConfig.for_pbmc1k(base_dir=upstream, workspace_dir=output)
    notebook = output / "notebooks" / "18_pbmc1k_template_downstream_analysis.ipynb"
    generate_template_notebook(config, notebook)
    nb = nbformat.read(notebook, as_version=4)
    if len(nb.cells) < 10 or nb.cells[1].cell_type != "code" or "PipelineConfig.for_pbmc1k" not in nb.cells[1].source:
        raise ValueError("The unpublished package template changed; review the bridge before execution")
    nb.cells[1].source = first_cell(package_src, upstream, output)
    nbformat.write(nb, notebook)
    kernel_for_current_python(output, package_src)

    manifest = {
        "status": "generated_only" if args.generate_only else "running",
        "dataset_id": "pbmc1k",
        "sample_id": "pbmc_1k_v3",
        "upstream_inputs": inputs,
        "downstream_source_sha256": source_fingerprint(package_src),
        "notebook": str(notebook.relative_to(output)),
        "baseline_metrics": BASELINE,
    }
    if not args.generate_only:
        result = execute_and_export_notebook(notebook, kernel_name="scpolaseq_bridge")
        summary_file = config.output_dir / "checkpoints" / "pipeline_execution_summary.json"
        summary = json.loads(summary_file.read_text(encoding="utf-8"))
        manifest["execution"] = {
            "notebook": str(notebook.relative_to(output)),
            "html": str(notebook.with_suffix('.html').relative_to(output)),
            "elapsed_seconds": result["elapsed_seconds"],
            "inline_plot_cells": result["inline_plot_cells"],
        }
        manifest["observed_metrics"] = {name: summary.get(name) for name in BASELINE}
        manifest["comparison"] = {
            name: {"baseline": old, "observed": summary.get(name), "matches": summary.get(name) == old}
            for name, old in BASELINE.items()
        }
        reactome_file = config.output_dir / "enrichment" / "reactome_pathway_enrichment.tsv"
        with reactome_file.open(newline="", encoding="utf-8") as handle:
            pathways = list(csv.DictReader(handle, delimiter="\t"))
        significant = sum(float(row["FDR"]) <= 0.05 for row in pathways)
        manifest["reactome_check"] = {
            "displayed_pathways": len(pathways),
            "pathways_fdr_le_0_05": significant,
            "package_checkpoint_n_reactome_pathways_fdr05": summary.get("n_reactome_pathways_fdr05"),
            "checkpoint_count_is_displayed_rows": summary.get("n_reactome_pathways_fdr05") == len(pathways),
        }
        manifest["status"] = "completed"
    (output / "bridge_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
