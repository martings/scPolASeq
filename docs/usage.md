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

## Expected key outputs

- `provenance/alignment/*.manifest.tsv`
- `labels/group_map.tsv`
- `apa/site_catalog.tsv`
- `apa/apa_usage.tsv`
- `apa/apa_stats.tsv`
- `coverage/*/tracks/*.bedGraph`
- `report/scpolaseq_report.html`
