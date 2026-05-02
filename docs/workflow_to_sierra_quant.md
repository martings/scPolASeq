# Workflow to Sierra Quant

## Purpose

This document explains the execution path from sample sheet ingestion to `SIERRA_QUANT` for the `pbmc1k` validation run, including the relevant processes, parameters, intermediate contracts, and improvement opportunities.

Run context described here matches the deepthought PBMC config and the current sample-merged execution strategy:

- Config: `conf/test_pbmc1k_deepthought.config`
- Results: `/scratch/results/scPolASeq/pbmc1k`
- Work directory: `/scratch/work/scPolASeq_runs/pbmc1k_issue10`

## Top-Level Workflow

The entrypoint is `workflow SCPOLASEQ` in `main.nf`.

High-level stages:

1. `INPUT_HARMONIZATION`
2. `REFERENCE_PREPARE`
3. `ALIGN_OR_IMPORT_SC`
4. `CELL_LABELING_SC`
5. `APA_CORE`
6. `REPORTING_SC`

The path into Sierra is:

`main.nf -> CELL_LABELING_SC -> APA_CORE -> BARCODE_FILTERING -> GROUPED_RECONSTRUCTION -> SIERRA_QUANT`

## Run Parameters Used for PBMC1K

From `conf/test_pbmc1k_deepthought.config`:

### Inputs

- `input = /scratch/data/pbmc1k/samplesheet.csv`
- `genome_fasta = /scratch/data/reference/GRCh38/genome.fa`
- `gtf = /scratch/data/reference/GRCh38/genes.gtf`
- `prebuilt_star_index = /scratch/results/scPolASeq/pbmc1k/reference/star_index`
- `known_polya = /scratch/data/polya_atlas/polyasite2_GRCh38.tsv`

### Output and execution mode

- `outdir = /scratch/results/scPolASeq/pbmc1k`
- `run_mode = full`
- `aligner = starsolo`
- `protocol_mode = 10x_3p`

### Clustering and grouping

- `enable_internal_clustering = true`
- `grouping_levels = cluster,cell_type`
- `apa_grouping_levels = cluster,cell_type`
- `apa_group_level = cluster`

### Coverage and APA

- `emit_group_bams = false`
- `emit_bigwigs = true`
- `apa_min_coverage = 5`
- `apa_min_pdui_delta = 0.2`
- `apa_model_type = random_forest`
- `enable_single_cell_apa_projection = false`
- `enable_pas_reference_build = true`
- `enable_sierra_quant = true`
- `enable_pas_scoring = false`

### Resources

- `max_cpus = 16`
- `max_memory = 120 GB`
- `max_time = 24 h`

## Stage-by-Stage Flow

## 1. Input Harmonization

Process family:

- `INPUT_HARMONIZATION`

Role:

- normalizes the sample sheet and optional external metadata inputs
- produces the sample descriptors consumed by alignment and labeling stages

Important contract:

- each sample/library travels downstream as `meta`
- `meta` is expected to carry at least `sample_id` and `library_id`

## 2. Reference Preparation

Process family:

- `REFERENCE_PREPARE`

Role:

- validates/stages genome FASTA and GTF
- exposes the reference bundle
- passes the known polyA atlas into later APA stages

Reference bundle consumed later by `APA_CORE` contains:

- STAR index
- GTF
- FASTA
- chromosome sizes
- terminal exons
- atlas
- blacklist

## 3. Alignment or Import

Process family:

- `ALIGN_OR_IMPORT_SC`
- for this run, the active aligner path is STARsolo

Relevant runtime parameters:

- `run_mode = full`
- `aligner = starsolo`
- `protocol_mode = 10x_3p`

Role:

- merge validated samplesheet rows by `sample_id` before STARsolo
- generate one alignment output set and one 10x-like expression matrix per sample
- preserve CB/UB tags in BAM
- emit:
  - `bam_bundle`
  - `matrix_bundle`

For the PBMC1K example, the input rows:

- `pbmc_1k_v3_L001`
- `pbmc_1k_v3_L002`

are collapsed into a single STARsolo job for:

- `pbmc_1k_v3`

Key outputs used downstream:

- `bam_bundle`: sample BAM + BAI
- `matrix_bundle`: sample STARsolo matrix directory with `barcodes.tsv`

## 4. Cell Labeling and Clustering

Subworkflow:

- `CELL_LABELING_SC`

Included modules:

- `SCANPY_CLUSTER_SC`
- `BUILD_GROUP_MAP`
- `WRITE_MANIFESTS`

### 4.1 `SCANPY_CLUSTER_SC`

Input tuple:

- `meta`
- `matrix_dir`
- `cell_metadata`
- `enable_internal_clustering`
- `celltypist_model`
- `reference_h5ad`
- `reference_label_col`

CLI arguments passed by the module:

- `--matrix-dir`
- `--cell-metadata`
- `--sample-id`
- `--library-id`
- `--enable-internal-clustering`
- optional:
  - `--celltypist-model`
  - `--reference-h5ad`
  - `--reference-label-col`
- outputs:
  - `--out-annotations`
  - `--out-report`
  - `--out-h5ad`
  - `--out-plot-prefix`

Role:

- load STARsolo matrix
- cluster cells internally
- attach canonical barcode fields
- emit per-library cell annotations
- emit clustered H5AD
- emit PCA and UMAP plots for inspection

Published outputs in `results/.../labels`:

- `<library>.cell_annotations.tsv`
- `<library>.cluster_report.tsv`
- `<library>.clustered.h5ad`
- `<library>.clustered.pca.png`
- `<library>.clustered.umap.png`

Canonical identity fields now expected in the H5AD and annotations:

- `sample_id`
- `library_id`
- `barcode_raw`
- `barcode_corrected`
- `cell_id`

### 4.2 `BUILD_GROUP_MAP`

Role:

- merges per-library annotation tables
- writes:
  - unified `cell_annotations.tsv`
  - unified `group_map.tsv`
  - group summary outputs

Important contract:

- `group_map.tsv` must retain `library_id`
- joins should be library-aware whenever the same 10x barcode appears in multiple libraries

## 5. APA Core Entry

Subworkflow:

- `APA_CORE`

Inputs into `APA_CORE` from `main.nf`:

- `reference_bundle`
- `bam_bundle`
- unified `cell_annotations.tsv`
- unified `group_map.tsv`
- `known_polya`
- `apa_grouping_levels`
- `apa_min_coverage`
- `apa_min_pdui_delta`
- `apa_model_type`
- `apa_model_group_level`

Feature toggles resolved at runtime:

- `enable_single_cell_apa_projection`
- `enable_pas_reference_build`
- `enable_sierra_quant`
- `enable_pas_scoring`

## 6. Barcode Filtering

Subworkflow:

- `BARCODE_FILTERING`

Included module:

- `BAM_BARCODE_FILTER`

Role:

- takes library BAMs from `bam_bundle`
- filters them against unified `cell_annotations.tsv`
- projects cell metadata back into BAM records through tags used later by grouping

Expected behavior after the issue #10 fix:

- annotation lookup is constrained by `sample_id` and `library_id`
- barcode mapping is driven by canonical barcode fields instead of reconstructed names

In sample-merged runs:

- `library_id` mirrors `sample_id`
- the filtered BAM already represents the whole sample

Outputs:

- filtered BAM + BAI
- filter stats

## 7. Grouped Reconstruction

Subworkflow:

- `GROUPED_RECONSTRUCTION`

Included module:

- `GROUPED_BAM_GENERATION`

Inputs:

- filtered BAM
- unified `group_map.tsv`
- grouping levels from `apa_grouping_levels`

For the PBMC run, active grouping levels are:

- `cluster`
- `cell_type`

### 7.1 `GROUPED_BAM_GENERATION`

CLI arguments passed:

- `--bam`
- `--group-map`
- `--sample-id`
- `--library-id`
- `--group-level`
- `--out-dir`
- `--out-manifest`

Role:

- split one filtered BAM into one BAM per `(library_id, group_level, group_id)`
- index every grouped BAM
- emit a grouping manifest

Current group tag behavior in `write_group_bams.py`:

- `cluster` uses BAM tag `CL`
- `cell_type` uses BAM tag `CT`
- `cell` uses BAM tag `CB`

Current naming behavior:

- grouped BAMs are prefixed with the active analysis identifier
- example: `pbmc_1k_v3.cluster_8.grouped.bam`

Why this matters:

- the filename prefix still prevents collisions in published outputs
- downstream code must not assume that the BAM basename equals the semantic `group_id`

## 8. Coverage and Feature Extraction

Parallel branches inside `APA_CORE`:

- `COVERAGE_GENERATION`
- `APA_FEATURE_PIPELINE`
- `MODEL_PIPELINE`
- `PAS_REFERENCE_BUILD`

These stages are not the direct cause of the current Sierra mismatch, but they depend on the same grouped BAM and annotation contracts.

Relevant outputs feeding reporting and APA interpretation:

- bedGraphs and bigWigs
- feature tables
- APA event tables
- PDUI matrix
- optional PAS reference

## 9. Sierra Quantification

Module:

- `SIERRA_QUANT`

Container:

- `${params.apptainer_cache_dir}/scpolaseq-r.sif`

For the current profile:

- `/scratch/cache/apptainer/scpolaseq/scpolaseq-r.sif`

Input tuple:

- `meta`
- `group_level`
- `group_id`
- `grouped_bam`
- `grouped_bai`
- `pas_reference`
- `cell_annotations`
- `gtf`

CLI arguments passed to `bin/sierra_quant.R`:

- `--bam`
- `--gtf`
- `--pas-reference`
- `--cell-annotations`
- `--library-id`
- `--group-level`
- `--group-id`
- `--output`
- `--log`
- `--ncores`

Role:

- build a group-scoped whitelist from annotations
- call `Sierra::CountPeaks`
- write:
  - `<library>.<group_level>.<group_id>.sierra_quant.tsv`
  - manifest
  - log

## Current Sierra Behavior

### Corrected behavior

The Sierra script now scopes annotations before counting.

Expected log signals for a healthy run:

- `Scoped annotations rows: N`
- `Whitelist barcodes after scoping: N`
- `Matrix loaded: ... x N cells`

Validated real example:

- workdir `/scratch/work/scPolASeq_runs/pbmc1k_issue10/2f/a88623ba24f3ed52fcac4ab74accb7`
- `Scoped annotations rows: 12`
- `Whitelist barcodes after scoping: 12`
- `Matrix loaded: 741003 peaks x 12 cells`

### Historical failure mode

`APA_CORE` currently derives `group_id` from the grouped BAM filename:

```groovy
def group_id = bam.baseName.replaceAll(/\.grouped(\.bam)?$/, '')
```

Because grouped BAMs are prefixed, this produced values like:

- `pbmc_1k_v3.cluster_8`

But `cluster_id` values in `cell_annotations.tsv` remain:

- `cluster_8`

So Sierra would receive a `group_id` that does not match the annotation schema, and scoping would collapse to zero rows.

Observed log pattern in failed rerun tasks on 2026-05-02 before the prefix-normalization patch:

- `Scoped annotations rows: 0`
- `Whitelist barcodes after scoping: 0`
- `WARNING: scoped whitelist is empty — writing empty TSV`

Result in published outputs before the patch:

- older non-prefixed cluster outputs remain populated
- newer prefixed outputs are header-only

## Recommendations

### Contract recommendations

- make `group_id` an explicit data field propagated in channels
- do not reconstruct semantic identifiers from filenames
- keep `library_id` and `group_id` orthogonal

### Workflow recommendations

1. Change `GROUPED_RECONSTRUCTION` or `GROUPED_BAM_GENERATION` to emit tuples containing explicit `group_id`.
2. Keep the current `APA_CORE` prefix normalization only as a compatibility bridge.
3. Keep filename prefixing for collision resistance, but treat it as presentation only.

### Validation recommendations

- add a pre-Sierra assertion that scoped annotations are non-empty for every emitted grouped BAM
- add a regression test where:
  - two libraries share barcode values
  - grouped BAM filenames are prefixed
  - Sierra cluster scoping still succeeds
- add a post-run report counting header-only Sierra outputs by `library_id` and `group_level`

### UX and debugging recommendations

- retain the H5AD PCA/UMAP publication in `labels/`
- add a compact `sierra_scope_summary.tsv` with:
  - `library_id`
  - `group_level`
  - `group_id`
  - `scoped_annotation_rows`
  - `sampled_cb_tags`
  - `output_rows`
- include a short warning in the final report when duplicate-style Sierra outputs exist for the same logical group

## Summary

The workflow up to `SIERRA_QUANT` is now much more explicit and sample-aware than before the issue #10 work. The main residual fragility is that `group_id` is still inferred from a filename; the new prefix-normalization patch makes that safe for now, but explicit identifier propagation would still be a cleaner contract.

Once that identifier propagation is made explicit, the sample-level STARsolo strategy should extend cleanly into Sierra for both `cluster` and `cell_type` groupings.
