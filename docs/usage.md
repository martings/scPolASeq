# scPolASeq usage

## Full run from FASTQ

```bash
nextflow run . -profile docker -params-file params/full.json --input samplesheet.csv --genome_fasta genome.fa --gtf genes.gtf
```

## Reuse aligned BAMs

```bash
nextflow run . -profile slurm,singularity --run_mode from_bam --input samplesheet.csv --cell_metadata cell_metadata.tsv
```

## Feedback mode after relabeling

```bash
nextflow run . -profile slurm,singularity --run_mode feedback --input samplesheet.csv --cell_metadata relabeled_cells.tsv
```

## pbmc1k MVP validation run (deepthought)

The PBMC 1k v3 (10x Genomics healthy donor) dataset is the canonical test for
the MVP stack on deepthought. Data, reference (GRCh38 + Gencode v44), and the
PolyASite 2.0 atlas are staged under `/scratch/data/`.

The STAR genome index from previous runs is cached; use `-resume` to skip the
~30 min genome generation step.

```bash
cd /scratch/work/scPolASeq

RUNDIR=/scratch/work/scPolASeq_runs/pbmc1k_$(date +%Y%m%d%H%M%S)

nextflow run main.nf \
  -profile deepthought,apptainer,test_pbmc1k_deepthought \
  -work-dir "$RUNDIR" \
  --outdir /scratch/results/scPolASeq/pbmc1k \
  -resume
```

### What to check after the run

| File | What to verify |
|---|---|
| `alignment/pbmc_1k_v3_L00*.Aligned.sortedByCoord.out.bam` | Non-zero size; `samtools flagstat` shows >80% aligned |
| `labels/group_map.tsv` | Multiple cluster/cell_type rows, not just one fake barcode |
| `coverage/cluster/*.bedGraph` | Non-empty; `wc -l` >1000 lines |
| `apa_calls/apa_events.tsv` | Header + at least some significant events |
| `apa_calls/pdui_usage_matrix.tsv` | Non-empty matrix |
| `pas_reference/pas_reference.tsv` | Should have many more rows than site_catalog |
| `report/scpolaseq_report.html` | Renders without error |

### Known non-fatal warnings

- `Failed to render execution report/timeline` — R not installed; HTML reports skipped
- `SIERRA_QUANT` and `PAS_SCORING` are disabled in this profile (scaffold stubs)

