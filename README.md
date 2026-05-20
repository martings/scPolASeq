# scPolASeq

`scPolASeq` is a Nextflow DSL2 pipeline for single-cell alternative polyadenylation (APA) discovery and quantification. It processes 10x Chromium scRNA-seq data from raw FASTQs through alignment, barcode-aware BAM reconstruction, per-cluster/cell-type APA quantification, and ML-based APA event scoring.

## Key features

- **STARsolo** alignment with automatic `.gz` detection and per-library BAM output.
- **Barcode-aware grouping**: reconstructs per-cluster and per-cell-type BAMs from the filtered cell x UMI matrix.
- **PAS reference build**: integrates a site catalog (annotation-guided candidate sites), terminal exon catalog, and a known polyA atlas (PolyASite 2.0 / PolyA_DB 4.1) into a unified `pas_reference.tsv`.
- **Sierra::CountPeaks** per-group APA quantification (R, Bioconductor).
- **SCAPTURE** de novo PAS detection sidecar (optional, disabled by default).
- **APA calling + ML model**: PDUI-based differential usage detection and a random forest classifier trained on per-site coverage features.
- **Coverage tracks**: strand-specific bigWig / bedGraph per group level.
- Multi-sample, multi-condition support via the samplesheet.

## Pipeline overview

```
FASTQs / BAMs
    └─ VALIDATE_INPUTS / HARMONIZE_CELL_METADATA
    └─ PREPARE_REFERENCE_BUNDLE  ← genome.fa + GTF + polyA atlas → STAR index
    └─ STARSOLO_ALIGN_SC         ← per library
    └─ SCANPY_CLUSTER_SC         ← internal clustering (optional)
    └─ BAM_BARCODE_FILTER        ← keep only whitelisted barcodes
    └─ GROUPED_BAM_GENERATION    ← reconstruct BAMs per cluster / cell type
    └─ BAM_MERGE_SAMPLE          ← merge runs of the same sample (multi-run)
    └─ SCAPTURE_FILTER           ← de novo PAS calling (opt-in)
    └─ SIERRA_QUANT              ← Sierra::CountPeaks per group
    └─ COVERAGE_TRACKS           ← strand-specific bigWig/bedGraph
    └─ APA_FEATURE_EXTRACTION    ─┐
    └─ APA_CALLING               ─┤ → apa_events.tsv / pdui_usage_matrix.tsv
    └─ APA_MODEL_TRAIN/PREDICT   ─┘ → scored_apa_events.tsv
```

## Quick start

### 1. Requirements

- Nextflow ≥ 24.04
- Apptainer ≥ 1.0 (or Docker/Singularity)
- Pre-built container images (see `containers/build.sh`)

### 2. Download the PBMC 1k v3 test dataset

The **PBMC 1k v3** dataset (10x Genomics, healthy donor, 10x Chromium 3' v3) is the canonical test dataset for the pipeline.

```bash
# Downloads FASTQs (~5 GB) and the 10x v3 barcode whitelist.
# Pass any local directory as the destination.
bash conf/download_pbmc1k.sh /data/pbmc1k
```

Or manually:

```bash
# FASTQs (2 lanes)
mkdir -p /data/pbmc1k/fastqs
curl -L https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar \
  | tar -xf - -C /data/pbmc1k/fastqs --strip-components=1

# 10x v3 barcode whitelist
curl -L https://cf.10xgenomics.com/supp/cell-exp/whitelists/3M-february-2018.txt.gz \
  | gunzip > /data/pbmc1k/3M-february-2018.txt
```

Human reference data (GRCh38 + GENCODE v44 + PolyASite 2.0) is staged by:

```bash
bash conf/setup_human_references.sh /data/references
```

### 3. Build container images

```bash
bash containers/build.sh --cache-dir /data/apptainer/scpolaseq
```

Four SIF images are built: `scpolaseq-alignment`, `scpolaseq-clustering`, `scpolaseq-python`, `scpolaseq-r`.  
A fifth image `scpolaseq-scapture` is built only if `--with-scapture` is passed (required for `enable_scapture = true`).

### 4. Create a samplesheet

The samplesheet is a CSV with one row per sequencing library (one `library_id` may have multiple rows if sequenced across multiple lanes — they will be merged automatically):

```csv
sample_id,library_id,protocol,chemistry,condition,replicate_id,fastq_r1,fastq_r2,barcode_whitelist
pbmc_1k_v3,pbmc_1k_v3_L001,10x_3p,10xv3,healthy,rep1,/path/L001_R1.fastq.gz,/path/L001_R2.fastq.gz,/path/3M-february-2018.txt
pbmc_1k_v3,pbmc_1k_v3_L002,10x_3p,10xv3,healthy,rep1,/path/L002_R1.fastq.gz,/path/L002_R2.fastq.gz,/path/3M-february-2018.txt
```

| Column | Description |
|--------|-------------|
| `sample_id` | Biological sample identifier. Libraries with the same `sample_id` are merged after filtering. |
| `library_id` | Unique sequencing run ID. Used for BAM naming and provenance. |
| `protocol` | `10x_3p` (only supported protocol currently). |
| `chemistry` | `10xv2` or `10xv3` — controls STARsolo CB+UMI lengths. |
| `condition` | Experimental group (e.g. `healthy`, `DF`, `DHF`). Used for grouping levels. |
| `replicate_id` | Biological replicate label. |
| `fastq_r1` / `fastq_r2` | Absolute paths to R1 (barcode+UMI) and R2 (cDNA) FASTQs. Gzip optional. |
| `barcode_whitelist` | Path to 10x whitelist `.txt`. Leave empty to run without whitelist (knee-point filtering). |
| `scapture_group_id` | *(optional)* Groups libraries for SCAPTURE peak merging across runs. |

### 5. Run the pipeline

**Test run (PBMC 1k v3, Apptainer):**

```bash
nextflow run main.nf \
  -profile apptainer,test_pbmc1k \
  -work-dir /data/work/scPolASeq_pbmc1k \
  -resume
```

**Custom run:**

```bash
nextflow run main.nf \
  -profile apptainer \
  --input samplesheet.csv \
  --genome_fasta /path/to/genome.fa \
  --gtf /path/to/genes.gtf \
  --known_polya /path/to/polyasite2_GRCh38.tsv \
  --outdir /path/to/results \
  --prebuilt_star_index /path/to/star_index    # skip genome generation if index exists
```

## Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `prebuilt_star_index` | `null` | Path to an existing STAR genome directory. When set, genome generation is skipped (~10–20 min and 28 GB saved). |
| `grouping_levels` | `cluster,cell_type` | Comma-separated list of grouping levels for BAM reconstruction and coverage. |
| `apa_group_level` | `cluster` | Primary grouping level for APA feature extraction and ML training. |
| `enable_internal_clustering` | `true` | Run Scanpy clustering when no external labels are provided. |
| `enable_sierra_quant` | `true` | Run Sierra::CountPeaks APA quantification. |
| `enable_scapture` | `false` | Run SCAPTURE de novo PAS detection (requires `scpolaseq-scapture.sif`). |
| `enable_ml_model` | `true` | Train and apply random forest APA classifier. |
| `apa_min_coverage` | `5` | Minimum per-site read coverage for APA calling. |
| `apa_min_pdui_delta` | `0.2` | Minimum |ΔPDUI| threshold for APA significance. |
| `scapture_deeppass_threshold` | `0.5` | Minimum DeepPASS score for retaining a SCAPTURE site. |

## Output structure

```
<outdir>/
├── alignment/          STARsolo BAMs and Solo.out cell x gene matrices
├── filtered_bam/       Barcode-filtered BAMs (one per library)
├── grouped_bam/        BAMs reconstructed per cluster and cell type
├── labels/             cell_annotations.tsv, group_map.tsv, clustered h5ad
├── provenance/         barcode_registry.tsv, input_manifest.json
├── pas_reference/      pas_reference.tsv (unified PAS site catalog)
├── apa/                site_catalog.tsv/bed (annotation-guided candidates)
├── sierra_quant/       Per-group Sierra::CountPeaks TSVs
├── scapture/           SCAPTURE peak catalog and PASquant output (if enabled)
├── coverage/           bigWig/bedGraph tracks per group level
├── apa_features/       apa_features.tsv (per-site coverage features)
├── apa_calls/          apa_events.tsv, pdui_usage_matrix.tsv, scored_apa_events.tsv
├── models/             apa_model.pkl, model_metrics.tsv, feature_importance.tsv
├── qc/                 barcode_filter_stats per library
├── report/             scpolaseq_report.html
└── pipeline_info/      execution_trace.txt, execution_report.html, timeline
```

## Tested datasets

| Dataset | Organism | Chemistry | Conditions | Status |
|---------|----------|-----------|------------|--------|
| PBMC 1k v3 (10x Genomics) | Human (GRCh38) | 10x v3 | Healthy | ✅ Full pipeline |

## Further reading

- [Usage guide](docs/usage.md) — detailed run instructions and parameter reference
- [Output specification](docs/outputs.md) — file formats and field definitions
- [Module contracts](docs/module_contracts.md) — input/output contracts per DSL2 module
- [Design notes](docs/design.md) — architectural decisions
