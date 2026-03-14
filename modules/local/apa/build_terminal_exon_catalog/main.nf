process BUILD_TERMINAL_EXON_CATALOG {
    tag "build_terminal_exon_catalog"
    label 'process_single'
    label 'process_python'

    conda "${projectDir}/envs/python.yml"
    publishDir "${params.outdir}/reference", mode: params.publish_dir_mode

    input:
    path gtf

    output:
    path "terminal_exons.tsv", emit: terminal_exons
    path "terminal_exons.bed", emit: terminal_exons_bed

    script:
    """
    python ${projectDir}/bin/build_terminal_exons.py \\
        --gtf ${gtf} \\
        --out-tsv terminal_exons.tsv \\
        --out-bed terminal_exons.bed
    """

    stub:
    """
    printf "gene_id\\ttranscript_id\\tchrom\\tstart\\tend\\tstrand\\texon_number\\tterminal_exon_rank\\n" > terminal_exons.tsv
    touch terminal_exons.bed
    """
}
