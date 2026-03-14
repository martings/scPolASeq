// Stage 3 subworkflow — Barcode-aware BAM filtering.
// Cell annotations from CELL_LABELING feed back into BAM to associate
// each read with its cluster, cell_type, and condition.
include { BAM_BARCODE_FILTER } from '../../../modules/local/bam/bam_barcode_filter/main'

workflow BARCODE_FILTERING {
    take:
    bam_bundle       // channel: tuple(meta, bam, bai)
    cell_annotations // path — single unified annotations TSV broadcast to every library

    main:
    // Broadcast the single annotations file to every BAM entry
    bam_bundle
        .combine(cell_annotations)
        .map { meta, bam, bai, annotations ->
            tuple(meta, bam, bai, annotations)
        }
        .set { ch_filter_input }

    BAM_BARCODE_FILTER(ch_filter_input)

    emit:
    filtered_bam = BAM_BARCODE_FILTER.out.filtered_bam   // tuple(meta, bam, bai)
    filter_stats = BAM_BARCODE_FILTER.out.filter_stats    // path
}
