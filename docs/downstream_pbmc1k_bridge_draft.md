# PBMC1k downstream notebook bridge (draft)

This branch is an isolated experiment based on `main`. It adds one Nextflow process, `PBMC1K_DOWNSTREAM`, launched through `workflows/pbmc1k_downstream.nf`. It is separate from scientific workflow PR #19 and does not run as part of `main.nf`.

The downstream package is now published at [martings/scPolASeq_downstream](https://github.com/martings/scPolASeq_downstream). This branch pins commit `365c0722758f7d3a08e9a06050c974207374cfe0` as a Git submodule under `external/scPolASeq_downstream`. The process accepts an existing scPolASeq PBMC1k result directory and invokes the pinned package's `PipelineConfig`, `generate_template_notebook`, and `execute_and_export_notebook`. The package itself now generates a notebook with the requested input and output paths, so the bridge no longer edits its code cells. A fresh output directory prevents accidental reuse of historical checkpoints.

## Required inputs

- `labels/pbmc_1k_v3.clustered.h5ad`
- `alignment/pbmc_1k_v3.Solo.out/Gene/filtered/matrix.mtx`
- Nonempty Sierra TSVs in `sierra_quant/cluster/`
- `apa/site_catalog.tsv`
- `apa_calls/apa_events.tsv` and `pdui_usage_matrix.tsv`
- Initialized `external/scPolASeq_downstream` submodule, or an explicit `--downstream_package_src` override
- Python 3.11 environment with the downstream dependencies and `ipykernel`; for full enrichment, the external science skill scripts and `uv` are also needed

Run from the scPolASeq repository root after initializing the pinned submodule:

```bash
git submodule update --init external/scPolASeq_downstream
nextflow run workflows/pbmc1k_downstream.nf \
  --upstream_results /path/to/pbmc1k_results \
  --downstream_python /path/to/scpolaseq_downstream_env/bin/python \
  --downstream_outdir /path/to/new_downstream_output \
  -work-dir /path/to/new_nextflow_work
```

Use `--downstream_generate_only true` to validate the connection and generate the notebook without executing it. The output directory must not exist before launch; Nextflow rejects it to prevent publication over a previous run. Nextflow stages the existing input directory; this module does not rerun upstream alignment or APA tasks. The package source defaults to the pinned submodule; `--downstream_package_src` remains available for development overrides. A full run checks for the three external science scripts and `uv` before notebook execution. `bridge_manifest.json` records the source fingerprint, required inputs, run status, notebook paths relative to the published output, and comparison metrics. A full run also writes `replication_check.md`, with the six reference counts, observed counts, notebook and HTML presence, figure count, and Reactome FDR count. Its technical concordance status is a file and summary check; scientific review remains pending.

## Observed replication on the October PBMC1k output

The full module completed on the October 24–25 PBMC1k execution output with Nextflow 25.10.4 and the Python 3.11 downstream environment after switching to the published package submodule. On October 26, the notebook executed in 226.45 seconds, exported HTML, and produced seven figures. The package source SHA-256 fingerprint was `fbbe651e88c7a6666c2404ec422471dc6c4f5188f5a0e3b1998927364cbaa3ad` at this test. All six reported counts matched the historical notebook 18 checkpoint:

| Metric | Historical notebook 18 | Fresh October output |
|---|---:|---:|
| QC-pass cells | 1,113 | 1,113 |
| Curated cell labels | 16 | 16 |
| Gene × lineage SI profiles | 22,130 | 22,130 |
| Differential APA comparisons | 52,750 | 52,750 |
| Package-defined robust APA events | 942 | 942 |
| Genes among those events | 614 | 614 |

Matching summary counts support technical reproducibility for this package version and input pair. They do not independently validate cell labels, PAS assignments, differential tests, or therapeutic claims. The upstream October run's 18 BigWigs were empty; this downstream notebook uses Sierra counts, the PAS catalog, STARsolo matrices, and the clustered AnnData, so those BigWigs were not required for this replication.

**Enrichment correction:** the original package checkpoint wrote `n_reactome_pathways_fdr05 = len(reactome_df)`, which counted displayed pathways. The earlier fresh table displayed 20 rows and had **zero** pathways with FDR ≤ 0.05 (minimum FDR 0.923). The published package now applies the FDR threshold in the checkpoint. In the pinned full run, the table again displayed 20 rows and both the checkpoint and the bridge independently reported **0** pathways with FDR ≤ 0.05. The bridge's `replication_check.md` reports technical concordance as complete and scientific review as pending.

## Integration limits

- The package template now accepts the configured upstream and workspace locations. The bridge verifies that the generated notebook contains both paths before executing it.
- The bridge targets the PBMC1k result schema. Other datasets and cross-condition comparisons require explicit adapters and independent validation.
- The current package may query Reactome, QuickGO, and PubMed during a full run through external science skill scripts. Their responses and versions are outside the upstream Nextflow provenance contract. The scripts are not bundled in the package repository; set `SCPOLASEQ_SCIENCE_SKILLS_DIR` where needed.
- The bridge leaves upstream `main.nf` and PR #19 untouched. Publication or inclusion in the default workflow needs a separate review of package versioning, output contracts, and scientific interpretation.
