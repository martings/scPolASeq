# Design notes

## Defaults

- `STARsolo` is the default alignment path for supported 10x-style inputs.
- Annotation-guided APA site generation is the default detection path.
- Molecule collapse happens during site quantification, not as BAM-wide deduplication.
- `feedback` mode is intended to skip realignment and regenerate grouped evidence after updated labeling.

## Grouping strategy

`group_map.tsv` stores barcode-to-group assignments independently of the BAM. This keeps cell labels replaceable while preserving aligned read provenance.

## Provenance strategy

The scaffold emits manifest tables at alignment and labeling stages. These are meant to remain stable inputs to later `feedback` runs.

## APA orchestration boundary

`subworkflows/apa_core.nf` is the contract freezer for APA-stage execution. It
owns channel reshaping, toggle handling, and future sidecar module wiring. This
keeps `main.nf` focused on coarse orchestration and prevents later feature work
from leaking tuple-layout assumptions across the pipeline.

## Contract strategy

Current working outputs remain the primary public surface. New modules
(`PAS_REFERENCE_BUILD`, `SIERRA_QUANT`, `PAS_SCORING`) are scaffolded as
sidecars with explicit I/O so they can mature independently without breaking
existing consumers.
