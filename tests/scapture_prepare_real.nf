nextflow.enable.dsl = 2

include { SCAPTURE_PREPARE } from '../modules/local/apa/scapture_prepare/main'

workflow {
    SCAPTURE_PREPARE(Channel.of(tuple(
        [library_id: 'tiny'],
        file("${projectDir}/data/scapture_duplicate_genes.gtf", checkIfExists: true),
        file("${projectDir}/../bin/prepare_scapture_gtf.py", checkIfExists: true)
    )))
    SCAPTURE_PREPARE.out.gene_name_map.view { meta, mapping -> "gene map: ${mapping}" }
}
