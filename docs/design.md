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
