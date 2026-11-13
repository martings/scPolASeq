nextflow.enable.dsl = 2

include { SCAPTURE_CATALOG } from '../modules/local/apa/scapture_catalog/main'
include { SCAPTURE_PASQUANT } from '../modules/local/apa/scapture_pasquant/main'

params.scapture_task_dir = null
params.scapture_library_id = 'pbmc_1k_v3'
params.scapture_test_catalog_only = false
params.scapture_test_fail_after_quant = false
params.scapture_test_sample_id = null

process SCAPTURE_TEST_GATE {
    input:
    tuple val(meta), path(matrix)
    val fail_after_quant

    output:
    path 'test_gate.ok'

    script:
    """
    test -s ${matrix}
    if [ '${fail_after_quant}' = 'true' ]; then
        echo 'Intentional test failure after real PASquant' >&2
        exit 42
    fi
    touch test_gate.ok
    """
}

workflow {
    if (!params.scapture_task_dir) {
        error 'Supply --scapture_task_dir with the completed historical PBMC1k task'
    }
    def prefix = params.scapture_library_id
    def root = file(params.scapture_task_dir, checkIfExists: true)
    def source = { name -> file("${root}/${name}", checkIfExists: true) }
    def meta = [library_id: prefix, sample_id: params.scapture_test_sample_id ?: prefix]
    SCAPTURE_CATALOG(Channel.of(tuple(meta,
        source("${prefix}.exonic.peaks.evaluated.bed"),
        source("${prefix}.intronic.peaks.evaluated.bed"),
        source("${prefix}.3primeExtended.peaks.evaluated.bed"),
        file("${projectDir}/../bin/run_scapture_filter.py", checkIfExists: true))))
    if (!params.scapture_test_catalog_only) {
        SCAPTURE_PASQUANT(Channel.of(tuple(meta,
            source("${prefix}.merged.bam"),
            source("${prefix}.merged.bam.bai"),
            source("${prefix}_merge.Integrated.bed"),
            source('cell_annotations.tsv'),
            file("${projectDir}/../bin/FilterBamByTag", checkIfExists: true))))
        SCAPTURE_TEST_GATE(SCAPTURE_PASQUANT.out.matrix, Channel.value(params.scapture_test_fail_after_quant))
    }
}
