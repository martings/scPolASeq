// Stage 4 subworkflow — Grouped BAM reconstruction.
// Produces one BAM per (group_level, group_id), enabling cluster/cell_type-level APA.
include { GROUPED_BAM_GENERATION } from '../../../modules/local/grouping/grouped_bam_generation/main'

workflow GROUPED_RECONSTRUCTION {
    take:
    filtered_bam    // channel: tuple(meta, bam, bai)
    group_map       // channel: path(group_map_tsv)  — single value
    grouping_levels // val: comma-separated string, e.g. 'cluster,cell_type'

    main:
    // Expand to one task per sampling level
    filtered_bam
        .combine(group_map)
        .combine(Channel.from(grouping_levels.tokenize(',')))
        .map { meta, bam, bai, gmap, glevel ->
            tuple(meta, bam, bai, gmap, glevel.trim())
        }
        .set { ch_group_input }

    GROUPED_BAM_GENERATION(ch_group_input)

    emit:
    grouped_bams      = GROUPED_BAM_GENERATION.out.grouped_bams      // tuple(meta, group_level, [bams], [bais])
    grouping_manifest = GROUPED_BAM_GENERATION.out.grouping_manifest  // path
}
