# PBMC1k downstream notebook bridge (draft)

This branch is an isolated experiment based on `main`. It adds one Nextflow process, `PBMC1K_DOWNSTREAM`, launched through `workflows/pbmc1k_downstream.nf`. It is separate from scientific workflow PR #19 and does not run as part of `main.nf`.

The process accepts an existing scPolASeq PBMC1k result directory and the source directory of the unpublished `scpolaseq_downstream` package. It invokes the package's `PipelineConfig`, `generate_template_notebook`, and `execute_and_export_notebook`. The bridge replaces the generated notebook's first code cell with explicit source and output paths; the analytical cells still call the package's QC, CellTypist, APA, enrichment, and plotting modules. A fresh output directory prevents accidental reuse of historical checkpoints.

## Required inputs

- `labels/pbmc_1k_v3.clustered.h5ad`
- `alignment/pbmc_1k_v3.Solo.out/Gene/filtered/matrix.mtx`
- Nonempty Sierra TSVs in `sierra_quant/cluster/`
- `apa/site_catalog.tsv`
- `apa_calls/apa_events.tsv` and `pdui_usage_matrix.tsv`
- Local, unpublished `scpolaseq_downstream/src` with its Python dependencies and `ipykernel` installed in the selected interpreter

Run from the scPolASeq repository root:

```bash
nextflow run workflows/pbmc1k_downstream.nf \
  --upstream_results /path/to/pbmc1k_results \
  --downstream_package_src /path/to/scPolASeq_downstream/src \
  --downstream_python /path/to/scpolaseq_downstream_env/bin/python \
  --downstream_outdir /path/to/new_downstream_output \
  -work-dir /path/to/new_nextflow_work
```

Use `--downstream_generate_only true` to validate the connection and generate the notebook without executing it. The output directory must be new. Nextflow stages the existing input directory; this module does not rerun upstream alignment or APA tasks. The package is imported from the supplied source directory; no published wheel is assumed. `bridge_manifest.json` records the source fingerprint, required inputs, run status, notebook paths relative to the published output, and comparison metrics. A full run also writes `replication_check.md`, with the six reference counts, observed counts, notebook and HTML presence, figure count, and a corrected Reactome FDR count. Its technical concordance status is a file and summary check; scientific review remains pending.

## Observed replication on the October PBMC1k output

The full module completed on the October 24–25 PBMC1k execution output with Nextflow 25.10.4 and the local Python 3.11 downstream environment. The notebook executed in 161.65 seconds, exported HTML, and produced seven figures. The package source SHA-256 fingerprint was `710ccc133cee3655ef9d19c91bddf36f54b5f180d7140f978faf2dc8f7062709` at this test. All six reported counts matched the historical notebook 18 checkpoint:

| Metric | Historical notebook 18 | Fresh October output |
|---|---:|---:|
| QC-pass cells | 1,113 | 1,113 |
| Curated cell labels | 16 | 16 |
| Gene × lineage SI profiles | 22,130 | 22,130 |
| Differential APA comparisons | 52,750 | 52,750 |
| Package-defined robust APA events | 942 | 942 |
| Genes among those events | 614 | 614 |

Matching summary counts support technical reproducibility for this package version and input pair. They do not independently validate cell labels, PAS assignments, differential tests, or therapeutic claims. The upstream October run's 18 BigWigs were empty; this downstream notebook uses Sierra counts, the PAS catalog, STARsolo matrices, and the clustered AnnData, so those BigWigs were not required for this replication.

**Enrichment correction:** the package checkpoint writes `n_reactome_pathways_fdr05 = len(reactome_df)`, which counts displayed pathways. The fresh table displays 20 rows and has **zero** pathways with FDR ≤ 0.05 (minimum FDR 0.923). The bridge's manifest computes the thresholded count from the table separately. The package owner should correct the checkpoint field before using it in scientific reports.

## Integration limits

- The unpublished package's template currently hard-codes its source and workspace paths. This bridge overrides only the first generated code cell and checks the template shape before execution. The package should expose these locations as supported parameters before this integration is promoted.
- The bridge targets the PBMC1k result schema. Other datasets and cross-condition comparisons require explicit adapters and independent validation.
- The current package may query Reactome, QuickGO, and PubMed during a full run. Their responses and versions are outside the upstream Nextflow provenance contract.
- The bridge leaves upstream `main.nf` and PR #19 untouched. Publication or inclusion in the default workflow needs a separate review of package versioning, output contracts, and scientific interpretation.
