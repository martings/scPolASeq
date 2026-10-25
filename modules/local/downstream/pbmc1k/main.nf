process PBMC1K_DOWNSTREAM {
    tag 'pbmc1k_downstream'
    label 'process_single'
    publishDir params.downstream_outdir, mode: 'copy'

    input:
    path upstream_results, stageAs: 'upstream_results'
    path package_src, stageAs: 'package_src'
    path bridge_script, stageAs: 'bridge.py'

    output:
    path 'downstream', emit: results

    script:
    def generate_only = params.containsKey('downstream_generate_only') && params.downstream_generate_only ? '--generate-only' : ''
    """
    ${params.downstream_python} bridge.py \\
      --upstream-results upstream_results \\
      --package-src package_src \\
      --output-dir downstream ${generate_only}
    """
}
