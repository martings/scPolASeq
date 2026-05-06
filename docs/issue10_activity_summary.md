# Issue #10 Activity Summary

## Scope

This document summarizes the work completed for GitHub issue #10 on branch `10-barcode-mismatch-star`, with emphasis on the barcode identity fix, the downstream impact on `SIERRA_QUANT`, and the later shift to sample-scoped processing from `STARsolo` onward.

Reference dataset and rerun context:

- Results directory: `/scratch/results/scPolASeq/pbmc1k`
- Active work directory: `/scratch/work/scPolASeq_runs/pbmc1k_issue10`
- Cached rerun observed in `screen`: `issue10_cached_20260502_165823`

## Problem Statement

The original failure mode was a barcode identity mismatch across:

- STARsolo `Solo.out/Gene/filtered/barcodes.tsv`
- `barcode_registry.tsv`
- `cell_annotations.tsv`
- `*.clustered.h5ad`
- grouping / coverage / Sierra inputs

This was especially visible in the `pbmc1k` sample with technical libraries `L001` and `L002`, where `L002` failed to map consistently and produced incomplete QC.

## Changes Implemented

### 1. Canonical barcode schema

Barcode identity is now preserved explicitly instead of being reconstructed from `obs_names` or regexes.

Canonical fields introduced and propagated:

- `sample_id`
- `library_id`
- `barcode_raw`
- `barcode_corrected`
- `cell_id`

Canonical cell identity:

- `cell_id = "{sample_id}:{library_id}:{barcode_corrected}"`

Design choice:

- `barcode_corrected` remains the canonical STARsolo barcode
- `obs_names` may be made globally unique, but `barcode_corrected` remains the reversible source of truth

### 2. Clustering outputs hardened

Updated `bin/scanpy_cluster_sc.py` and its module wrapper so that:

- the H5AD stores explicit barcode columns in `adata.obs`
- `library_id` is retained throughout clustering outputs
- `cell_annotations.tsv` contains library-aware barcode fields
- dimensionality reduction plots are emitted from the final H5AD

Published outputs now include:

- `*.clustered.h5ad`
- `*.clustered.pca.png`
- `*.clustered.umap.png`

Confirmed in `/scratch/results/scPolASeq/pbmc1k/labels`:

- `pbmc_1k_v3_L001.clustered.h5ad`
- `pbmc_1k_v3_L001.clustered.pca.png`
- `pbmc_1k_v3_L001.clustered.umap.png`
- `pbmc_1k_v3_L002.clustered.h5ad`
- `pbmc_1k_v3_L002.clustered.pca.png`
- `pbmc_1k_v3_L002.clustered.umap.png`

### 3. Grouping and BAM projection made sample-aware

Updated downstream scripts so barcode-to-group joins remain stable after moving the analysis unit from lane-level `library_id` to merged `sample_id`.

Key changes:

- `group_map.tsv` still carries `library_id` for backward compatibility
- BAM filtering can still respect library scoping, but sample-merged STARsolo now emits `library_id == sample_id`
- grouped BAM generation uses the correct BAM tag by grouping mode:
  - `cluster -> CL`
  - `cell_type -> CT`
  - `cell -> CB`
- grouped BAM filenames are prefixed with the active analysis identifier, which now mirrors `sample_id`

### 4. Sierra whitelist scoping fixed

`bin/sierra_quant.R` now scopes `cell_annotations.tsv` before building the whitelist.

Current scoping logic:

- filter annotations by `library_id`
- then restrict by:
  - `cluster_id == group_id` for `group_level=cluster`
  - `cell_type == group_id` for `group_level=cell_type`
  - `cell_id == group_id` for `group_level=cell`

This removed the earlier bug where Sierra counted a grouped BAM with a full-library whitelist.

### 5. `STARsolo` now merges rows by `sample_id`

The pipeline now groups validated samplesheet rows by `sample_id` before entering `STARSOLO_ALIGN_SC`.

For `pbmc1k`, the original rows:

- `pbmc_1k_v3_L001`
- `pbmc_1k_v3_L002`

are now merged into one STARsolo execution keyed by:

- `sample_id = pbmc_1k_v3`

Downstream consequences:

- one STARsolo matrix per sample
- one alignment BAM per sample
- one H5AD per sample
- one unified barcode namespace per sample
- no downstream need to reconcile cells across lane-specific STARsolo outputs

## Container Validation

The Sierra module uses the Apptainer image configured by the pipeline:

- `container "${params.apptainer_cache_dir}/scpolaseq-r.sif"`

For the deepthought/apptainer profile used here, that resolves to:

- `/scratch/cache/apptainer/scpolaseq/scpolaseq-r.sif`

Validation completed inside that container:

- `Rscript` parse check for `bin/sierra_quant.R`: passed
- manual rerun of one Sierra workdir using the patched script: passed

Observed successful scoped run:

- Workdir: `/scratch/work/scPolASeq_runs/pbmc1k_issue10/2f/a88623ba24f3ed52fcac4ab74accb7`
- Log highlights:
  - `Scoped annotations rows: 12`
  - `Whitelist barcodes after scoping: 12`
  - `BAM CB tags sampled: 12 unique`
  - `Matrix loaded: 741003 peaks x 12 cells`
  - `Output written: recheck.cluster_8.sierra_quant.tsv`

## Current Rerun Status as of 2026-05-02

### What is corrected

- H5ADs are present and published
- dimensionality plots are present and published
- barcode identity is explicit in annotations and H5AD
- Sierra no longer defaults to a full-library whitelist when the inputs are scoped correctly

### What is still failing

The cached rerun on 2026-05-02 exposed a new mismatch in the APA chain.

Observed in recent Sierra logs under `/scratch/work/scPolASeq_runs/pbmc1k_issue10`:

- `Scoped annotations rows: 0`
- `Whitelist barcodes after scoping: 0`
- `WARNING: scoped whitelist is empty — writing empty TSV`

Examples:

- `.../72/de1b1aff255c209488be932929884b/pbmc_1k_v3_L002.cluster.pbmc_1k_v3_L002.cluster_1.sierra_quant.log`
- `.../ff/a1ce358d7945f5b9c47e850add047d/pbmc_1k_v3_L002.cluster.pbmc_1k_v3_L002.cluster_8.sierra_quant.log`
- `.../0b/a2bd3e4635cabbb252a358d0c8878f/pbmc_1k_v3_L001.cell_type.pbmc_1k_v3_L001.unlabeled.sierra_quant.log`

Published outputs show the same pattern:

- legacy files like `pbmc_1k_v3_L002.cluster.cluster_8.sierra_quant.tsv` remain populated
- newly prefixed files like `pbmc_1k_v3_L002.cluster.pbmc_1k_v3_L002.cluster_8.sierra_quant.tsv` are header-only

## Root Cause of the Remaining Sierra Issue

The problem is in `subworkflows/apa_core.nf`.

Current logic derives `group_id` from the grouped BAM filename:

```groovy
def group_id = bam.baseName.replaceAll(/\.grouped(\.bam)?$/, '')
```

After the grouped BAM naming change, filenames became analysis-prefixed, for example:

- `pbmc_1k_v3.cluster_8.grouped.bam`

That means `group_id` becomes:

- `pbmc_1k_v3.cluster_8`

But the annotation table still stores:

- `cluster_id = cluster_8`

So the Sierra scoping filter would compare:

- `cluster_id == pbmc_1k_v3.cluster_8`

instead of:

- `cluster_id == cluster_8`

That bug is now patched by stripping the active analysis prefix before calling Sierra.

Result before the patch:

- the whitelist becomes empty
- Sierra emits an empty TSV by design

## Recommendations

### Immediate fix

The active patch now strips the `{analysis_id}.` prefix before calling `SIERRA_QUANT`.

Preferred long-term option:

1. Pass `group_id` directly from the grouping manifest generated by `GROUPED_BAM_GENERATION`.

### Stability improvements

- Treat manifest metadata as authoritative and filenames as presentation only.
- Add a preflight assertion before Sierra:
  - `group_id` must match at least one row in scoped `cell_annotations.tsv`
- Add a regression test covering:
  - two libraries sharing the same barcode namespace
  - grouped BAM filenames prefixed by `library_id`
  - Sierra scoping for both `cluster` and `cell_type`
- Consider publishing grouped BAMs in per-sample subdirectories even when filenames are already prefixed.

### Observability improvements

- Keep `Scoped annotations rows` and `Whitelist barcodes after scoping` in Sierra logs
- Add a post-run QC summary that flags header-only Sierra outputs
- Add a report table mapping:
  - `library_id`
  - `group_level`
  - `group_id`
  - `n_cells_in_annotations`
  - `n_cb_tags_in_bam`
  - `n_rows_in_sierra_output`

## Files Touched During the Issue #10 Fix

Core logic:

- `bin/scanpy_cluster_sc.py`
- `bin/build_group_map.py`
- `bin/grouped_3prime_coverage.py`
- `bin/filter_bam_by_barcodes.py`
- `bin/write_group_bams.py`
- `bin/harmonize_cell_metadata.py`
- `bin/sierra_quant.R`

Workflow and module wiring:

- `modules/local/clustering/scanpy_cluster_sc/main.nf`
- `modules/local/bam/bam_barcode_filter/main.nf`
- `modules/local/grouping/grouped_bam_generation/main.nf`
- `subworkflows/apa_core.nf`

Support material:

- `docs/module_contracts.md`
- `docs/outputs.md`
- `tests/test_issue10_barcode_identity.py`

## Bottom Line

The original barcode identity problem from STARsolo through clustering has been substantially fixed, and the workflow is now being shifted to sample-scoped STARsolo processing so the biological unit stays consistent from alignment onward.

The remaining Sierra failure observed in the 2026-05-02 cached rerun was localized to `group_id` derivation from prefixed grouped BAM filenames inside `APA_CORE`. That normalization is now patched, and the clean verification run should confirm the end-to-end behavior under the new sample-level strategy.
