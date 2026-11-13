nextflow.enable.dsl = 2

include { SCAPTURE_STAGED } from '../subworkflows/local/scapture_staged/main'

workflow {
    // Parses and wires the complete graph; biological validation uses real
    // PBMC1k inputs and is documented separately.
    SCAPTURE_STAGED(Channel.empty())
}
