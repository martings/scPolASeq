# Outputs

## Reference

- `reference/reference_manifest.tsv`
- `reference/terminal_exons.tsv`

## Alignment

- `alignment/<sample_id>.Aligned.sortedByCoord.out.bam`
- `alignment/<sample_id>.barcode_registry.tsv`
- `alignment/<sample_id>.alignment_manifest.tsv`

`barcode_registry.tsv` preserves STARsolo barcodes verbatim in
`barcode_raw` and `barcode_corrected`, including 10x suffixes such as `-1`.
When multiple FASTQ rows share the same `sample_id`, they are merged before
`STARsolo`; downstream alignment outputs are therefore sample-scoped.

## Labels and grouping

- `labels/cell_annotations.tsv`
- `labels/group_map.tsv`
- `labels/group_summary.tsv`

`cell_annotations.tsv` canonical columns:
`sample_id`, `library_id`, `barcode_raw`, `barcode_corrected`, `cell_id`,
`cluster_id`, `cell_type`, `condition`, `batch`, `label_source`

`group_map.tsv` canonical columns:
`sample_id`, `library_id`, `barcode_corrected`, `group_level`, `group_id`

`<sample_id>.clustered.h5ad` stores explicit `obs` columns
`sample_id`, `library_id`, `barcode_raw`, `barcode_corrected`, and `cell_id`.
If the AnnData index is rewritten for uniqueness, `barcode_corrected` remains
the canonical STARsolo barcode and `cell_id` remains reversible.

## Coverage

- `coverage/<sample_id>.grouped_3prime_counts.tsv`
- `coverage/<sample_id>.coverage_summary.tsv`
- `coverage/<sample_id>/tracks/*.bedGraph`
- `coverage/<sample_id>/tracks/*.bigWig`

## APA

- `apa/site_catalog.tsv`
- `apa/apa_usage.tsv`
- `apa/site_usage_matrix.tsv`
- `apa/gene_summary.tsv`
- `apa/apa_stats.tsv`

## PAS sidecars

- `pas_reference/pas_reference.tsv`
- `pas_reference/pas_reference_build.manifest.tsv`
- `sierra_quant/<group_level>/<sample_id>.<group_level>.<group_id>.sierra_quant.tsv`
- `pas_scoring/pas_scored_events.tsv`

## Report

- `report/scpolaseq_report.html`
- `report/report_plots/`
