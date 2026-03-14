// Stage 5 subworkflow — Strand-aware coverage tracks per grouped BAM.
include { COVERAGE_TRACKS } from '../../../modules/local/coverage/coverage_tracks/main'

workflow COVERAGE_GENERATION {
    take:
    grouped_bams  // channel: tuple(meta, group_level, [bam_files], [bai_files])
    chrom_sizes   // path

    main:
    // Flatten: one process call per (meta, group_level, group_id, bam, bai)
    grouped_bams
        .flatMap { meta, group_level, bams, bais ->
            def bam_list = (bams instanceof List) ? bams : [bams]
            def bai_list = (bais instanceof List) ? bais : [bais]
            [bam_list, bai_list].transpose().collect { bam, bai ->
                def gid = bam.baseName.replaceAll(/\.grouped(\.bam)?$/, '')
                tuple(meta, group_level, gid, bam, bai)
            }
        }
        .set { ch_flat_bams }

    COVERAGE_TRACKS(ch_flat_bams, chrom_sizes)

    emit:
    bedgraphs      = COVERAGE_TRACKS.out.bedgraphs      // tuple(meta, group_level, group_id, fwd_bg, rev_bg)
    bigwigs        = COVERAGE_TRACKS.out.bigwigs         // tuple(meta, group_level, group_id, fwd_bw, rev_bw)
    coverage_stats = COVERAGE_TRACKS.out.coverage_stats  // path
}
