#!/usr/bin/env python3
"""Run the pinned downstream package on one scPolASeq PBMC1k result tree."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
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


def write_replication_report(output: Path, manifest: dict) -> None:
    """Summarize technical concordance without treating it as biological validation."""
    execution = manifest["execution"]
    artifacts = {
        "notebook_nonempty": (output / execution["notebook"]).is_file()
        and (output / execution["notebook"]).stat().st_size > 0,
        "html_nonempty": (output / execution["html"]).is_file()
        and (output / execution["html"]).stat().st_size > 0,
    }
    figures_dir = output / "outputs" / "pbmc1k" / "figures"
    figures = sorted(path for path in figures_dir.glob("*.png") if path.stat().st_size > 0)
    artifacts["nonempty_figure_count"] = len(figures)
    artifacts["expected_figure_count"] = 7
    artifacts["figures_match_expected_count"] = len(figures) == 7
    metrics_match = all(item["matches"] for item in manifest["comparison"].values())
    artifact_checks_pass = all(
        artifacts[key] for key in ("notebook_nonempty", "html_nonempty", "figures_match_expected_count")
    )
    manifest["technical_concordance"] = {
        "matches_reference_summary": metrics_match,
        "artifact_checks": artifacts,
        "complete": metrics_match and artifact_checks_pass,
        "scientific_review_status": "pending",
    }

    lines = [
        "# PBMC1k downstream replication check",
        "",
        f"Technical concordance: **{'complete' if metrics_match and artifact_checks_pass else 'incomplete'}**.",
        "Scientific review: **pending**. Matching counts and file presence do not verify labels, PAS assignments, or biological claims.",
        "",
        "| Summary metric | Notebook 18 reference | New run | Match |",
        "|---|---:|---:|:---:|",
    ]
    for name, item in manifest["comparison"].items():
        lines.append(
            f"| `{name}` | {item['baseline']} | {item['observed']} | {'yes' if item['matches'] else 'no'} |"
        )
    reactome = manifest["reactome_check"]
    checkpoint_count = reactome["package_checkpoint_n_reactome_pathways_fdr05"]
    count_note = (
        "matches the FDR-filtered count" if checkpoint_count == reactome["pathways_fdr_le_0_05"]
        else "differs from the FDR-filtered count; review package checkpoint"
    )
    lines.extend([
        "",
        "## Generated artifacts",
        "",
        f"- Executed notebook: {'present' if artifacts['notebook_nonempty'] else 'missing or empty'}.",
        f"- HTML export: {'present' if artifacts['html_nonempty'] else 'missing or empty'}.",
        f"- Nonempty PNG figures: {len(figures)} of 7 expected.",
        "",
        "## Reactome count audit",
        "",
        f"- Displayed pathways: {reactome['displayed_pathways']}.",
        f"- Pathways with FDR <= 0.05: {reactome['pathways_fdr_le_0_05']}.",
        f"- Package checkpoint field `n_reactome_pathways_fdr05`: {checkpoint_count} ({count_note}).",
        "",
        f"Downstream package source SHA-256: `{manifest['downstream_source_sha256']}`.",
        "",
    ])
    (output / "replication_check.md").write_text("\n".join(lines), encoding="utf-8")


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
    os.environ["PYTHONPATH"] = os.pathsep.join(filter(None, [str(package_src), os.environ.get("PYTHONPATH")]))
    import nbformat
    from scpolaseq_downstream.config import PipelineConfig
    from scpolaseq_downstream.templates import generate_template_notebook, execute_and_export_notebook
    if not args.generate_only:
        from scpolaseq_downstream.enrichment import REACTOME_SCRIPT, QUICKGO_SCRIPT, PUBMED_SCRIPT, UV_BIN
        missing = [str(path) for path in (REACTOME_SCRIPT, QUICKGO_SCRIPT, PUBMED_SCRIPT) if not path.is_file()]
        if not UV_BIN.is_file() and shutil.which("uv") is None:
            missing.append("uv executable")
        if missing:
            raise ValueError("Full enrichment requires external science tools: " + ", ".join(missing))

    config = PipelineConfig.for_pbmc1k(base_dir=upstream, workspace_dir=output)
    notebook = output / "notebooks" / "18_pbmc1k_template_downstream_analysis.ipynb"
    generate_template_notebook(config, notebook)
    nb = nbformat.read(notebook, as_version=4)
    if (len(nb.cells) < 10 or nb.cells[1].cell_type != "code"
            or str(upstream) not in nb.cells[1].source
            or str(output) not in nb.cells[1].source):
        raise ValueError("The package template did not use the requested input and output paths")

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
        result = execute_and_export_notebook(notebook)
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
        write_replication_report(output, manifest)
    (output / "bridge_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
