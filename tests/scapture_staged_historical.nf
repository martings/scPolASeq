nextflow.enable.dsl = 2

include { SCAPTURE_STAGED } from '../subworkflows/local/scapture_staged/main'

params.scapture_task_dir = null
params.scapture_library_id = 'pbmc_1k_v3'
params.scapture_test_fail_after_quant = false

process SCAPTURE_STAGED_TEST_GATE {
    input:
    tuple val(meta), path(matrix)
    val fail_after_quant

    output:
    path 'test_gate.ok'

    script:
    """
    test -s ${matrix}
    if [ '${fail_after_quant}' = 'true' ]; then
        echo 'Intentional failure after staged SCAPTURE quantification' >&2
        exit 42
    fi
    touch test_gate.ok
    """
}

workflow {
    if (!params.scapture_task_dir) {
        error 'Supply --scapture_task_dir with a completed frozen SCAPTURE task directory'
    }
    def root = file(params.scapture_task_dir, checkIfExists: true)
    def source = { name -> file("${root}/${name}", checkIfExists: true) }
    def id = params.scapture_library_id
    def meta = [library_id: id, sample_id: id]
    SCAPTURE_STAGED(Channel.of(tuple(meta,
        source("${id}.merged.bam"),
        source("${id}.merged.bam.bai"),
        source('reference.genome.fa'),
        source('reference.genome.fa.fai'),
        source('reference.annotation.gtf'),
        source('reference.chrom.sizes'),
        source('cell_annotations.tsv'),
        source('scapture_polya_db_input'),
        false)))
    SCAPTURE_STAGED_TEST_GATE(
        SCAPTURE_STAGED.out.quant,
        Channel.value(params.scapture_test_fail_after_quant)
    )
}
