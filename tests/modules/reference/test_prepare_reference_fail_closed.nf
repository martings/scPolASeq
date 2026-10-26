nextflow.enable.dsl = 2

include { PREPARE_REFERENCE_BUNDLE } from '../../../modules/local/reference/prepare_reference_bundle/main'

params.outdir = params.outdir ?: "${launchDir}/tests/results/fail_closed"
params.publish_dir_mode = 'copy'
params.apptainer_cache_dir = null

workflow {
    PREPARE_REFERENCE_BUNDLE(
        Channel.value(file("${launchDir}/tests/data/reference.fa")),
        Channel.value(file("${launchDir}/tests/data/reference.gtf")),
        Channel.value(file("${launchDir}/assets/NO_FILE")),
        Channel.value(file("${launchDir}/assets/NO_FILE")),
        Channel.value(file("${launchDir}/assets/NO_FILE"))
    )
}
