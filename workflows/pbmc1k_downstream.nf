nextflow.enable.dsl = 2

include { PBMC1K_DOWNSTREAM } from '../modules/local/downstream/pbmc1k/main'

params.upstream_results = null
params.downstream_package_src = null
params.downstream_outdir = null
params.downstream_python = 'python3'
params.downstream_generate_only = false

workflow {
    if (!params.upstream_results || !params.downstream_package_src || !params.downstream_outdir) {
        error 'Set --upstream_results, --downstream_package_src, and --downstream_outdir'
    }
    PBMC1K_DOWNSTREAM(
        Channel.value(file(params.upstream_results, checkIfExists: true)),
        Channel.value(file(params.downstream_package_src, checkIfExists: true)),
        Channel.value(file("${projectDir}/../bin/pbmc1k_downstream_bridge.py", checkIfExists: true))
    )
}
