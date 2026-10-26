nextflow.enable.dsl = 2

include { PBMC1K_DOWNSTREAM } from '../modules/local/downstream/pbmc1k/main'

params.upstream_results = null
params.downstream_package_src = "${projectDir}/../external/scPolASeq_downstream/src"
params.downstream_outdir = null
params.downstream_python = 'python3'
params.downstream_generate_only = false

workflow {
    if (!params.upstream_results || !params.downstream_outdir) {
        error 'Set --upstream_results and --downstream_outdir'
    }
    def python_exec = params.downstream_python?.toString()
    if (!(python_exec == 'python3' ||
          (python_exec?.startsWith('/') && python_exec.matches('/[A-Za-z0-9_./-]+') &&
           !python_exec.tokenize('/').contains('..') && !python_exec.endsWith('/')))) {
        error 'Invalid --downstream_python: use python3 or an absolute executable path without shell characters'
    }
    try {
        java.nio.file.Files.createDirectory(file(params.downstream_outdir))
    } catch (java.nio.file.FileAlreadyExistsException ignored) {
        error "Downstream output already exists; choose a new --downstream_outdir: ${params.downstream_outdir}"
    }
    PBMC1K_DOWNSTREAM(
        Channel.value(file(params.upstream_results, checkIfExists: true)),
        Channel.value(file(params.downstream_package_src, checkIfExists: true)),
        Channel.value(file("${projectDir}/../bin/pbmc1k_downstream_bridge.py", checkIfExists: true)),
        Channel.value(params.downstream_generate_only)
    )
}
