nextflow.enable.dsl = 2

include { SCAPTURE_MERGE } from '../subworkflows/local/scapture_merge/main'

params.scapture_evaluated_dir = null
params.scapture_library_id = 'pbmc_1k_v3'
params.scapture_empty_extended = false

workflow {
    if (!params.scapture_evaluated_dir) {
        error 'Supply --scapture_evaluated_dir with three evaluated BEDs from a completed real run'
    }
    def prefix = params.scapture_library_id
    def source = file(params.scapture_evaluated_dir, checkIfExists: true)
    def bed = { category -> file("${source}/${prefix}.${category}.peaks.evaluated.bed", checkIfExists: true) }
    def extended = params.scapture_empty_extended
        ? file("${projectDir}/data/scapture_empty_evaluated.bed", checkIfExists: true)
        : bed('3primeExtended')
    SCAPTURE_MERGE(Channel.of(tuple([library_id: prefix], bed('exonic'), bed('intronic'), extended)))
    SCAPTURE_MERGE.out.reference.view { meta, reference -> "PAS reference: ${reference}" }
}
