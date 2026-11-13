# scPolASeq usage

## Full run from FASTQ

Use a host profile that matches the machine and one supported runtime profile.
The generic PBMC profile contains no staged paths; provide all inputs and
references for the host where the run will execute.

```bash
nextflow run main.nf \
  -profile apptainer,test_pbmc1k \
  --input /data/pbmc1k/samplesheet.csv \
  --genome_fasta /data/references/GRCh38/genome.fa \
  --gtf /data/references/GRCh38/genes.gtf \
  --known_polya /data/references/polyasite2_GRCh38.tsv \
  --apptainer_cache_dir /data/apptainer/scpolaseq \
  --outdir /data/results/pbmc1k \
  -work-dir /data/work/scPolASeq_pbmc1k
```

## Shared reference root and STAR index reuse

Generic workflow defaults keep `reference_root`, `reference_build`, `genome_fasta`, `gtf`, `chrom_sizes`, `star_index`, and `known_polya` as `null` so no host-specific paths are required or embedded in the core workflow.

To point a run or server profile at a shared reference directory (such as a fast SSD mount), set `--reference_root` (and optionally `--reference_build`, e.g. `GRCh38` or `GRCm38`). Any relative path passed via `--genome_fasta`, `--gtf`, `--chrom_sizes`, `--star_index`, or `--known_polya` is resolved under `--reference_root`; when `--genome_fasta` and `--gtf` are omitted with `--reference_root` + `--reference_build`, standard paths (`genomes/<build>/genome.fa`, `genomes/<build>/genes.gtf`, optional `genomes/<build>/chrom.sizes`, optional `indexes/<build>/star_index`) are resolved automatically:

```bash
nextflow run main.nf \
  -profile apptainer \
  --input /data/pbmc1k/samplesheet.csv \
  --reference_root /mnt/ref \
  --reference_build GRCh38 \
  --known_polya databases/polya_atlas/polyasite2_GRCh38.tsv \
  --outdir /data/results/pbmc1k
```

`PREPARE_REFERENCE_BUNDLE` rejects missing or empty inputs, validates optional `--chrom_sizes` against the FASTA, and records SHA-256 digests and tool/format versions in `reference/reference_manifest.json`. For `--star_index`, STAR `2.7.11b` uses the same `versionGenome = 2.7.4a` binary index format introduced in STAR `2.7.4a`; before reusing any prebuilt index, `PREPARE_REFERENCE_BUNDLE` checks required files (`Genome`, `SA`, `SAindex`, `genomeParameters.txt`), verifies `versionGenome` compatibility, confirms `chrNameLength.txt` matches `genome_fasta`, and records `index_sha256` in `reference_manifest.json`. Pass `--star_index null` to regenerate an index from FASTA + GTF.

## Reuse aligned BAMs

Use `--run_mode from_bam` with the samplesheet and cell metadata required by
that mode. Add the same host and runtime profiles used for the full run.

```bash
nextflow run main.nf \
  -profile apptainer \
  --run_mode from_bam \
  --input /data/project/samplesheet.csv \
  --cell_metadata /data/project/cell_metadata.tsv \
  --outdir /data/results/from_bam
```

## Feedback mode after relabeling

```bash
nextflow run main.nf \
  -profile apptainer \
  --run_mode feedback \
  --input /data/project/samplesheet.csv \
  --cell_metadata /data/project/relabeled_cells.tsv \
  --outdir /data/results/feedback
```

Confirm required inputs for each mode in `nextflow_schema.json` and the
workflow's parameter validation before launching.

## PolyA reference downloads

The workflow accepts PolyA_DB and PolyASite references as local files, URLs,
or the value `download`. When both sources are set, it downloads missing
archives, normalizes them, and merges nearby sites according to
`--pas_reference_merge_distance`. A single configured source is normalized
and used directly. Downloaded archives are cached under `<outdir>/reference_cache`
by default, or under `--pas_reference_cache_dir` when supplied.

```bash
nextflow run main.nf \
  -profile apptainer \
  --input /data/project/samplesheet.csv \
  --genome_fasta /data/references/GRCh38/genome.fa \
  --gtf /data/references/GRCh38/genes.gtf \
  --polya_db download \
  --polyasite download \
  --pas_reference_cache_dir /data/cache/pas \
  --outdir /data/results/pas_references
```

For a single source, use `--polya_db download --polyasite skip`, or set one
parameter to a local path/URL and the other to `skip`.

## Standalone PolyA reference helper

The helper script can download and normalize the built-in PolyA_DB and
PolyASite references. It stores the cache under the user's XDG data directory
by default; use `--outdir` to select another destination. Set the resulting
file through `--known_polya` when launching the pipeline.

```bash
bash bin/download_pas_references.sh --help
bash bin/download_pas_references.sh --outdir /data/references/polya_atlas
```

## PBMC 1k staged profile on Deepthought

`test_pbmc1k_deepthought` and `deepthought` contain host-specific paths and
resource settings. Use them only on the Deepthought system with the staged PBMC
inputs. The Deepthought host profile sets `reference_root = '/mnt/ref'` (SSD layout `genomes/<build>/`, `indexes/<build>/star_index`, `databases/polya_atlas/`). Override `--reference_root /tmp/fixture` with the same profile to test against a temporary fixture tree. The work directory and result directory should be unique to a run.

```bash
RUN_ID=$(date +%Y%m%d%H%M%S)
RUNDIR="/scratch/work/scPolASeq_runs/pbmc1k_${RUN_ID}"
OUTDIR="/scratch/results/scPolASeq/pbmc1k_${RUN_ID}"

nextflow run main.nf \
  -profile deepthought,apptainer,test_pbmc1k_deepthought \
  -work-dir "$RUNDIR" \
  --outdir "$OUTDIR" \
  --apptainer_cache_dir /scratch/cache/apptainer
```

The corresponding staged inputs are in `conf/test_pbmc1k_deepthought.config`.
Inspect the task trace and required output files after completion; a successful
Nextflow exit alone does not establish scientific QC or reproduce published
APA findings.
