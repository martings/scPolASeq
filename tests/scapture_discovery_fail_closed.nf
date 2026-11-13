nextflow.enable.dsl = 2

include { SCAPTURE_DISCOVERY } from '../subworkflows/local/scapture_discovery/main'

workflow {
    def fixture = { name -> file("${projectDir}/data/${name}", checkIfExists: true) }
    SCAPTURE_DISCOVERY(Channel.of(tuple(
        [library_id: 'tiny'],
        fixture('placeholder_input.bam'),
        fixture('scapture_dummy.bai'),
        fixture('scapture_annotation_small.fa'),
        fixture('scapture_annotation_small.fa.fai'),
        fixture('scapture_annotation_small.gtf'),
        fixture('scapture_annotation_small.chrom.sizes'),
        file("${projectDir}/../assets/NO_FILE", checkIfExists: true),
        false
    )))
}
