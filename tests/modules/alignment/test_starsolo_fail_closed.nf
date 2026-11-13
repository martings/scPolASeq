nextflow.enable.dsl = 2

include { STARSOLO_ALIGN_SC } from '../../../modules/local/alignment/starsolo_align_sc/main'

params.outdir = params.outdir ?: "${launchDir}/tests/results/fail_closed"
params.publish_dir_mode = 'copy'
params.apptainer_cache_dir = null

def meta = [sample_id: 'fail_closed', library_id: 'fail_closed', chemistry: '10xv3']

workflow {
    STARSOLO_ALIGN_SC(Channel.of(tuple(
        meta,
        [file("${launchDir}/tests/data/sample1_R1.fastq")],
        [file("${launchDir}/tests/data/sample1_R2.fastq")],
        file("${launchDir}/assets/NO_FILE"),
        file(params.test_star_index),
        file("${launchDir}/tests/data/reference.gtf"),
        '10x_3p'
    )))
}
