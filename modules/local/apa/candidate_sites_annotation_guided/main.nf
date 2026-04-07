process CANDIDATE_SITES_ANNOTATION_GUIDED {
    tag "candidate_sites_annotation_guided"
    label 'process_single'
    label 'process_python'

    conda "${projectDir}/envs/python.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-python.sif" : null
    publishDir "${params.outdir}/apa", mode: params.publish_dir_mode, pattern: "site_catalog*"

    input:
    tuple path(terminal_exons), path(known_polya)

    output:
    path "site_catalog.tsv", emit: site_catalog
    path "site_catalog.bed", emit: site_catalog_bed

    script:
    """
    python ${projectDir}/bin/candidate_sites_annotation_guided.py \\
        --terminal-exons ${terminal_exons} \\
        --known-polya ${known_polya} \\
        --out-tsv site_catalog.tsv \\
        --out-bed site_catalog.bed
    """

    stub:
    """
    printf "site_id\\tgene_id\\tchrom\\tstart\\tend\\tstrand\\tsite_class\\tsite_source\\tpriming_flag\\ngeneA:chr1:100:+\\tgeneA\\tchr1\\t100\\t100\\t+\\tannotated_terminal_end\\tannotation\\tunchecked\\n" > site_catalog.tsv
    printf "chr1\\t99\\t100\\tgeneA:chr1:100:+\\t0\\t+\\n" > site_catalog.bed
    """
}
