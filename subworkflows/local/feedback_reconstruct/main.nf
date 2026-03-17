include { GROUPED_3PRIME_COVERAGE } from '../../../modules/local/grouping/grouped_3prime_coverage/main'
include { WRITE_GROUP_BAMS        } from '../../../modules/local/grouping/write_group_bams/main'

workflow FEEDBACK_RECONSTRUCT {
    take:
    bam_bundle
    barcode_registry
    cell_annotations
    group_map
    reference_bundle
    emit_group_bams
    emit_bigwigs

    main:
    bam_bundle
        .join(barcode_registry)
        .combine(group_map)
        .combine(reference_bundle)
        .map { meta, bam, bai, registry, group_map_file, reference_meta, star_index, gtf, fasta, chrom_sizes, terminal_exons, atlas, blacklist ->
            tuple(meta, bam, bai, registry, group_map_file, chrom_sizes, emit_bigwigs)
        }
        .set { ch_coverage_input }

    GROUPED_3PRIME_COVERAGE(ch_coverage_input)

    def ch_group_bams = Channel.empty()
    def ch_qc_bundle = GROUPED_3PRIME_COVERAGE.out.coverage_summary.map { meta, summary -> summary }

    if (emit_group_bams) {
        bam_bundle
            .combine(group_map)
            .map { meta, bam, bai, group_map_file ->
                tuple(meta, group_map_file)
            }
            .set { ch_group_bam_input }
        WRITE_GROUP_BAMS(ch_group_bam_input)
        ch_group_bams = WRITE_GROUP_BAMS.out.group_bam_dir.map { meta, group_dir -> group_dir }
        ch_qc_bundle = ch_qc_bundle.mix(WRITE_GROUP_BAMS.out.manifest.map { meta, manifest -> manifest })
    }

    emit:
    coverage_bundle = GROUPED_3PRIME_COVERAGE.out.coverage_bundle.map { meta, coverage_file -> coverage_file }
    track_bundle    = GROUPED_3PRIME_COVERAGE.out.track_bundle.map { meta, track_dir -> track_dir }
    qc_bundle       = ch_qc_bundle
    group_bams      = ch_group_bams
    group_map       = group_map
}
