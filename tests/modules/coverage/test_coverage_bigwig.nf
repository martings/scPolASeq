nextflow.enable.dsl = 2

include { COVERAGE_TRACKS } from '../../../modules/local/coverage/coverage_tracks/main'

params.outdir = params.outdir ?: "${launchDir}/tests/results/coverage_bigwig"
params.publish_dir_mode = 'copy'
params.apptainer_cache_dir = null
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
    COVERAGE_TRACKS(bam_ch, Channel.value(file(params.test_chrom_sizes)))
}
