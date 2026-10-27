process PBMC1K_DOWNSTREAM {
    tag 'pbmc1k_downstream'
    label 'process_single'
    publishDir params.downstream_outdir, mode: 'copy'

    input:
    path upstream_results, stageAs: 'upstream_results'
    path package_src, stageAs: 'package_src'
    path bridge_script, stageAs: 'bridge.py'
    val generate_only

    output:
    path 'downstream', emit: results

    script:
    def generation_flag = generate_only ? '--generate-only' : ''
    """
    ${params.downstream_python} bridge.py \\
      --upstream-results upstream_results \\
      --package-src package_src \\
      --output-dir downstream ${generation_flag}
    """
}
