// Stage 5 subworkflow — Strand-aware bedGraph coverage tracks and optional BigWig conversion.
include { COVERAGE_TRACKS    } from '../../../modules/local/coverage/coverage_tracks/main'
include { BEDGRAPH_TO_BIGWIG } from '../../../modules/local/coverage/bedgraph_to_bigwig/main'

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

    COVERAGE_TRACKS(ch_flat_bams, Channel.value(file("${projectDir}/bin/bam_to_bedgraph.py", checkIfExists: true)))

    def do_bigwigs = params.emit_bigwigs == null ? true : params.emit_bigwigs.toString().toBoolean()
    def ch_bigwigs = Channel.empty()
    def ch_bigwig_status = Channel.empty()
    def ch_bigwig_versions = Channel.empty()

    if (do_bigwigs) {
        BEDGRAPH_TO_BIGWIG(COVERAGE_TRACKS.out.bedgraphs, chrom_sizes)
        ch_bigwigs = BEDGRAPH_TO_BIGWIG.out.bigwigs
        ch_bigwig_status = BEDGRAPH_TO_BIGWIG.out.status
        ch_bigwig_versions = BEDGRAPH_TO_BIGWIG.out.versions
    }

    emit:
    bedgraphs       = COVERAGE_TRACKS.out.bedgraphs      // tuple(meta, group_level, group_id, fwd_bg, rev_bg)
    bigwigs         = ch_bigwigs                         // tuple(meta, group_level, group_id, fwd_bw, rev_bw)
    bigwig_status   = ch_bigwig_status                   // path
    bigwig_versions = ch_bigwig_versions                 // path
    coverage_stats  = COVERAGE_TRACKS.out.coverage_stats  // path
}
