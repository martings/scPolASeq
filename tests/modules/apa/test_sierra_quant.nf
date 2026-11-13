nextflow.enable.dsl = 2

include { SIERRA_QUANT } from '../../../subworkflows/local/sierra_quant/main'

params.outdir = params.outdir ?: "${launchDir}/tests/results/sierra_quant"
params.publish_dir_mode = 'copy'
params.apptainer_cache_dir = params.apptainer_cache_dir ?: null
params.sample_id = params.containsKey('sample_id') ? params.sample_id : 'pbmc_1k_v3'
params.library_id = params.containsKey('library_id') ? params.library_id : 'pbmc_1k_v3'
params.group_level = params.containsKey('group_level') ? params.group_level : 'cluster'
params.group_id = params.containsKey('group_id') ? params.group_id : 'cluster_5'
params.strandedness = params.containsKey('strandedness') ? params.strandedness : null
params.chemistry = params.containsKey('chemistry') ? params.chemistry : null

workflow {
    def meta = [sample_id: params.sample_id ?: 'pbmc_1k_v3',
                library_id: params.library_id ?: 'pbmc_1k_v3',
                strandedness: params.strandedness,
                chemistry: params.chemistry]
    SIERRA_QUANT(
        Channel.of(tuple(meta, params.group_level ?: 'cluster', params.group_id ?: 'cluster_5',
            file(params.test_bam, checkIfExists: true),
            file(params.test_bai, checkIfExists: true),
            file(params.test_pas_reference, checkIfExists: true),
            file(params.test_cell_annotations, checkIfExists: true),
            file(params.test_gtf, checkIfExists: true)))
    )
}
