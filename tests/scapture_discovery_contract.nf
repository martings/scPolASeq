nextflow.enable.dsl = 2

include { SCAPTURE_DISCOVERY } from '../subworkflows/local/scapture_discovery/main'

workflow {
    // Compile the full module graph without manufacturing biological outputs.
    // Real-output equivalence is tracked separately in issue #29.
    SCAPTURE_DISCOVERY(Channel.empty())
}
