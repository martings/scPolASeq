process VALIDATE_INPUTS {
    tag "validate_inputs"
    label 'process_single'
    label 'process_python'

    conda "${projectDir}/envs/python.yml"
    publishDir "${params.outdir}/provenance", mode: params.publish_dir_mode, pattern: "samplesheet.validated.csv"
    publishDir "${params.outdir}/provenance", mode: params.publish_dir_mode, pattern: "input_manifest.json"

    input:
    path samplesheet
    path schema

    output:
    path "samplesheet.validated.csv", emit: samplesheet
    path "input_manifest.json", emit: manifest

    script:
    """
    python ${projectDir}/bin/validate_inputs.py \\
        --input ${samplesheet} \\
        --schema ${schema} \\
        --out samplesheet.validated.csv \\
        --manifest input_manifest.json
    """

    stub:
    """
    touch samplesheet.validated.csv
    echo "{}" > input_manifest.json
    """
}
