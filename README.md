# scPolASeq


`scPolASeq` is a modular Nextflow DSL2 scaffold for single-cell APA discovery and quantification with barcode-aware read reconstruction, grouped 3' coverage, and feedback re-analysis after clustering or cell-type labeling.

## Design goals

- `STARsolo`-first full alignment for supported 10x-like data.
- Explicit provenance from FASTQ/BAM through barcode registries, grouping, coverage, and APA statistics.
- Annotation-guided APA as the robust MVP path.
- Optional `Scanpy`-based clustering when external labels are not supplied.
- HPC- and container-friendly nf-core-like repository layout.

## Repository status

This scaffold includes:

- A coherent `main.nf` orchestration.
- Custom DSL2 modules and local subworkflows.
- Canonical input and output schemas.
- Python and R helper scripts for MVP logic.
- Test fixtures for `full`, `from_bam`, and `feedback` entry modes.

The repo intentionally keeps `modules/nf-core/` as a vendoring target and tracks preferred nf-core module equivalents in [`modules.json`](./modules.json), while remaining runnable from local custom modules.

## Quick start

See [`docs/usage.md`](./docs/usage.md) for example commands and [`docs/design.md`](./docs/design.md) for architecture notes.
