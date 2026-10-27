nextflow.enable.dsl = 2

include { SCAPTURE_MERGE } from '../subworkflows/local/scapture_merge/main'

workflow {
    def empty = file("${projectDir}/data/scapture_empty_evaluated.bed", checkIfExists: true)
    SCAPTURE_MERGE(Channel.of(tuple([library_id: 'empty'], empty, empty, empty)))
}
