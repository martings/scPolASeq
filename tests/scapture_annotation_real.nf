nextflow.enable.dsl = 2

include { SCAPTURE_PREPARE } from '../modules/local/apa/scapture_prepare/main'
include { SCAPTURE_ANNOTATION } from '../modules/local/apa/scapture_annotation/main'

workflow {
    def fixture = { name -> file("${projectDir}/data/${name}", checkIfExists: true) }
    SCAPTURE_PREPARE(Channel.of(tuple(
        [library_id: 'tiny'],
        fixture('scapture_annotation_small.gtf'),
        file("${projectDir}/../bin/prepare_scapture_gtf.py", checkIfExists: true)
    )))
    SCAPTURE_ANNOTATION(SCAPTURE_PREPARE.out.gtf.map { meta, gtf ->
        tuple(meta, gtf,
              fixture('scapture_annotation_small.fa'),
              fixture('scapture_annotation_small.fa.fai'),
              fixture('scapture_annotation_small.chrom.sizes'))
    })
    SCAPTURE_ANNOTATION.out.annotation.view { meta, files -> "annotation files: ${files}" }
}
