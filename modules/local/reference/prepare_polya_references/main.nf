process PREPARE_POLYA_REFERENCES {
    tag "prepare_polya_references"
    label 'process_single'
    label 'process_python'

    conda "${projectDir}/envs/python.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-python.sif" : null
    publishDir "${params.outdir}/reference", mode: params.publish_dir_mode

    input:
    path known_polya, stageAs: 'known_polya_input'
    val polya_db
    path polya_db_file, stageAs: 'polya_db_file_input'
    val polyasite
    path polyasite_file, stageAs: 'polyasite_file_input'

    output:
    path "prepared_known_polya.tsv", emit: known_polya
    path "prepared_known_polya.manifest.tsv", emit: manifest

    script:
    def skip_values = ['skip', 'none', 'null', 'false', 'no', '0', '']
    def polya_db_file_arg = polya_db?.toString()?.trim()?.toLowerCase() in skip_values ? '' : polya_db_file
    def polyasite_file_arg = polyasite?.toString()?.trim()?.toLowerCase() in skip_values ? '' : polyasite_file
    def default_cache_dir = params.outdir.toString().startsWith('/') ? "${params.outdir}/reference_cache" : "${launchDir}/${params.outdir}/reference_cache"
    def cache_dir = params.pas_reference_cache_dir ?: default_cache_dir
    """
    python3 ${projectDir}/bin/merge_polya_references.py \\
        --known-polya ${known_polya} \\
        --polya-db '${polya_db ?: 'skip'}' \\
        --polya-db-file '${polya_db_file_arg}' \\
        --polyasite '${polyasite ?: 'skip'}' \\
        --polyasite-file '${polyasite_file_arg}' \\
        --cache-dir '${cache_dir}' \\
        --merge-distance ${params.pas_reference_merge_distance ?: 25} \\
        --out-tsv prepared_known_polya.tsv \\
        --out-manifest prepared_known_polya.manifest.tsv
    """

    stub:
    """
    printf "gene_id\\tgene_name\\tchrom\\tstart\\tend\\tstrand\\tscore\\tsource\\n" > prepared_known_polya.tsv
    printf "geneA\\tgeneA\\tchr1\\t120\\t120\\t+\\t0\\tstub\\n" >> prepared_known_polya.tsv
    printf "field\\tvalue\\nmode\\tstub\\n" > prepared_known_polya.manifest.tsv
    """
}
