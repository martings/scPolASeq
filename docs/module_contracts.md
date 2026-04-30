# Module Contracts

This document freezes the DSL2 contracts for the APA backbone. The goal is to
keep orchestration stable even as biological logic evolves.

## Naming conventions

- `meta` is a Nextflow map with stable keys:
  `sample_id`, `library_id`, `protocol`, `chemistry`, `condition`, `replicate_id`
- `reference_bundle` is a tuple with fixed positional layout:
  `tuple val(ref_meta), path(star_index), path(gtf), path(fasta), path(chrom_sizes), path(terminal_exons), path(site_catalog), path(priming_blacklist)`
- Group-resolved files always carry both `group_level` and `group_id` in the channel contract, even if only one appears in the filename.
- New PAS sidecar stages must not replace existing published APA outputs until explicitly promoted.

## Subworkflow boundary

### `APA_CORE`

Purpose:
connect barcode filtering, grouped BAM generation, coverage, APA feature calling,
model scoring, and future PAS-oriented sidecar stages behind one orchestration
surface.

Input contract:

```nextflow
reference_bundle                   // tuple(ref_meta, star_index, gtf, fasta, chrom_sizes, terminal_exons, atlas, blacklist)
bam_bundle                         // tuple val(meta), path(bam), path(bai)
cell_annotations                   // path("cell_annotations.tsv")
group_map                          // path("group_map.tsv")
known_polya                        // path("*.tsv") or assets/NO_FILE
apa_grouping_levels                // val("cluster,cell_type,...")
apa_min_coverage                   // val(Integer)
apa_min_pdui_delta                 // val(Decimal)
apa_model_type                     // val(String)
apa_model_group_level              // val(String)
```

Primary emitted contracts:

```nextflow
feature_table         = path("apa_features.tsv")
apa_events            = path("apa_events.tsv")
pdui_matrix           = path("pdui_usage_matrix.tsv")
scored_events         = path("scored_apa_events.tsv")
track_bundle          = path("*.bw")
qc_bundle             = path("*.tsv")
pas_reference         = path("pas_reference.tsv") | fallback path(site_catalog)
sierra_quant          = tuple val(meta), val(group_level), val(group_id), path("*.sierra_quant.tsv")
pas_scored_events     = path("pas_scored_events.tsv")
```

## Stable stage contracts

### `BARCODE_FILTERING`

```nextflow
input:
tuple val(meta), path(bam), path(bai)
path("cell_annotations.tsv")

output:
tuple val(meta), path("*.bam"), path("*.bam.bai")
path("*filter*.tsv")
```

Notes:
`meta.library_id` remains the library-scoped identity anchor. Filtering must not
rewrite sample identities. Unified `cell_annotations.tsv` inputs must carry
explicit `sample_id`, `library_id`, `barcode_raw`, `barcode_corrected`, and a
reversible `cell_id = <sample_id>:<library_id>:<barcode_corrected>`.

### `GROUPED_RECONSTRUCTION`

```nextflow
input:
tuple val(meta), path(filtered_bam), path(filtered_bai)
path("group_map.tsv")
val(grouping_levels)

output:
tuple val(meta), val(group_level), path("*.grouped.bam"), path("*.grouped.bam.bai")
path("*.grouping_manifest.tsv")
```

Notes:
the output channel is library-scoped. One emission may contain multiple grouped
BAM paths for the same `(meta, group_level)` pair. When `group_map.tsv`
contains `library_id`, consumers must scope barcode/group lookups to the active
library before matching on `barcode_corrected`.

### `COVERAGE_GENERATION`

```nextflow
input:
tuple val(meta), val(group_level), path("*.grouped.bam"), path("*.grouped.bam.bai")
path("chrom.sizes")

output:
tuple val(meta), val(group_level), val(group_id), path("*.fwd.bedGraph"), path("*.rev.bedGraph")
tuple val(meta), val(group_level), val(group_id), path("*.fwd.bw"),       path("*.rev.bw")
path("*.coverage_stats.tsv")
```

Notes:
strand-separated outputs are mandatory. Forward and reverse tracks travel as one
logical contract unit.

### `APA_FEATURE_PIPELINE`

```nextflow
input:
path("site_catalog.tsv")
tuple val(meta), val(group_level), val(group_id), path("*.fwd.bedGraph"), path("*.rev.bedGraph")
path("cell_annotations.tsv")
path("known_polya.tsv") | path("NO_FILE")
val(min_coverage)
val(min_pdui_delta)

output:
path("apa_features.tsv")
path("apa_events.tsv")
path("pdui_usage_matrix.tsv")
```

Unified label tables:

- `cell_annotations.tsv` canonical columns:
  `sample_id`, `library_id`, `barcode_raw`, `barcode_corrected`, `cell_id`,
  `cluster_id`, `cell_type`, `condition`, `batch`, `label_source`
- `group_map.tsv` canonical columns:
  `sample_id`, `library_id`, `barcode_corrected`, `group_level`, `group_id`

### `MODEL_PIPELINE`

```nextflow
input:
path("apa_features.tsv")
path("apa_events.tsv")
val(model_type)
val(group_level)
val(enable_single_cell_projection)

output:
path("apa_model.pkl")
path("feature_importance.tsv")
path("model_metrics.tsv")
path("scored_apa_events.tsv")
```

## PAS sidecar modules

### `PAS_REFERENCE_BUILD`

```nextflow
input:
path("site_catalog.tsv")
path("terminal_exons.tsv")
path("known_polya.tsv") | path("NO_FILE")

output:
path("pas_reference.tsv")
path("pas_reference_build.manifest.tsv")
path("pas_reference_build.log")
```

Naming rule:
the emitted file is always exactly `pas_reference.tsv`, independent of the
source method.

Parameter control:
`params.pas_reference_merge_distance` controls evidence collapsing without
changing the channel contract.

### `SIERRA_QUANT`

```nextflow
input:
tuple val(meta), val(group_level), val(group_id), path("*.grouped.bam"), path("*.grouped.bam.bai"), path("pas_reference.tsv"), path("cell_annotations.tsv")

output:
tuple val(meta), val(group_level), val(group_id), path("*.sierra_quant.tsv")
path("*.sierra_quant.manifest.tsv")
path("*.sierra_quant.log")
```

Naming rule:
`<library_id>.<group_level>.<group_id>.sierra_quant.tsv`

### `PAS_SCORING`

```nextflow
input:
path("apa_features.tsv")
path("apa_events.tsv")
path("pdui_usage_matrix.tsv")
path("scored_apa_events.tsv")
path("pas_reference.tsv")

output:
path("pas_scored_events.tsv")
path("pas_scoring.manifest.tsv")
path("pas_scoring.log")
```

## Toggle behavior

- Toggle control is parameter-driven, not channel-driven:
  `params.enable_single_cell_apa_projection`,
  `params.enable_pas_reference_build`,
  `params.enable_sierra_quant`,
  `params.enable_pas_scoring`
- `enable_pas_reference_build = false`
  `APA_CORE` falls back to the existing annotation-guided site catalog as the
  effective PAS reference.
- `enable_sierra_quant = false`
  Sierra sidecar outputs emit an empty channel and do not affect legacy stages.
- `enable_pas_scoring = false`
  PAS scoring sidecar outputs emit an empty channel and do not affect reporting.
