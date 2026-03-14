process WRITE_MANIFESTS {
    tag "${manifest_name}"
    label 'process_single'
    label 'process_python'

    conda "${projectDir}/envs/python.yml"
    publishDir "${params.outdir}/provenance", mode: params.publish_dir_mode

    input:
    path input_table
    val manifest_name

    output:
    path "${manifest_name}.manifest.tsv", emit: manifest

    script:
    """
    python ${projectDir}/bin/write_manifest.py \\
        --input-table ${input_table} \\
        --manifest-name ${manifest_name} \\
        --out ${manifest_name}.manifest.tsv
    """

    stub:
    """
    printf "manifest_name\\tinput_table\\trows\\n${manifest_name}\\t${input_table.getName()}\\t0\\n" > ${manifest_name}.manifest.tsv
    """
}
