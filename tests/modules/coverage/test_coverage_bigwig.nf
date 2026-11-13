nextflow.enable.dsl = 2

include { COVERAGE_TRACKS    } from '../../../modules/local/coverage/coverage_tracks/main'
include { BEDGRAPH_TO_BIGWIG } from '../../../modules/local/coverage/bedgraph_to_bigwig/main'

params.outdir = params.outdir ?: "${launchDir}/tests/results/coverage_bigwig"
params.publish_dir_mode = 'copy'
params.apptainer_cache_dir = params.apptainer_cache_dir ?: null
params.emit_bigwigs = params.emit_bigwigs == null ? true : params.emit_bigwigs.toString().toBoolean()

workflow {
    def meta = [sample_id: 'coverage_test', library_id: 'coverage_test']
    def bam_ch = Channel.of(tuple(
        meta,
        'cluster',
        'cluster_1',
        file(params.test_bam),
        file(params.test_bai)
    ))
    COVERAGE_TRACKS(
        bam_ch,
        Channel.value(file("${projectDir}/../../../bin/bam_to_bedgraph.py", checkIfExists: true))
    )
    if (params.emit_bigwigs) {
        BEDGRAPH_TO_BIGWIG(
            COVERAGE_TRACKS.out.bedgraphs,
            Channel.value(file(params.test_chrom_sizes))
        )
    }
}
