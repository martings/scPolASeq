/*
 * Sierra quantification — wraps Sierra::CountPeaks.
 *
 * Contract:
 *   input:
 *     tuple val(meta), val(group_level), val(group_id), path(grouped_bam), path(grouped_bai), path(pas_reference), path(cell_annotations), path(gtf)
 *   output:
 *     tuple val(meta), val(group_level), val(group_id), path("<library>.<group_level>.<group_id>.sierra_quant.tsv")
 *     path "<library>.<group_level>.<group_id>.sierra_quant.manifest.tsv"
 *     path "<library>.<group_level>.<group_id>.sierra_quant.log"
 */
process SIERRA_QUANT {
    tag "${meta.library_id}:${group_level}:${group_id}"
    label 'process_medium'

    conda "${projectDir}/envs/r.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-r.sif" : null

    publishDir "${params.outdir}/sierra_quant/${group_level}", mode: params.publish_dir_mode

    input:
    tuple val(meta), val(group_level), val(group_id), path(grouped_bam), path(grouped_bai), path(pas_reference), path(cell_annotations), path(gtf)
    path quant_script

    output:
    tuple val(meta), val(group_level), val(group_id), path("${meta.library_id}.${group_level}.${group_id}.sierra_quant.tsv"), emit: quantified_pas
    path "${meta.library_id}.${group_level}.${group_id}.sierra_quant.manifest.tsv", emit: manifest
    path "${meta.library_id}.${group_level}.${group_id}.sierra_quant.log", emit: log

    script:
    def prefix = "${meta.library_id}.${group_level}.${group_id}"
    def ncores = task.cpus ?: 1
    """
    Rscript ${quant_script} \\
        --bam               ${grouped_bam} \\
        --gtf               ${gtf} \\
        --pas-reference     ${pas_reference} \\
        --cell-annotations  ${cell_annotations} \\
        --library-id        ${meta.library_id} \\
        --group-level       ${group_level} \\
        --group-id          ${group_id} \\
        --output            ${prefix}.sierra_quant.tsv \\
        --log               ${prefix}.sierra_quant.log \\
        --ncores            ${ncores}

    printf "field\\tvalue\\n"                                    > ${prefix}.sierra_quant.manifest.tsv
    printf "grouped_bam\\t${grouped_bam.name}\\n"              >> ${prefix}.sierra_quant.manifest.tsv
    printf "grouped_bai\\t${grouped_bai.name}\\n"              >> ${prefix}.sierra_quant.manifest.tsv
    printf "gtf\\t${gtf.name}\\n"                              >> ${prefix}.sierra_quant.manifest.tsv
    printf "pas_reference\\t${pas_reference.name}\\n"          >> ${prefix}.sierra_quant.manifest.tsv
    printf "cell_annotations\\t${cell_annotations.name}\\n"    >> ${prefix}.sierra_quant.manifest.tsv
    printf "quant_source\\tSierra::CountPeaks\\n"              >> ${prefix}.sierra_quant.manifest.tsv
    """

    stub:
    def prefix = "${meta.library_id}.${group_level}.${group_id}"
    """
    printf "library_id\\tgroup_level\\tgroup_id\\tsite_id\\tcell_barcode\\tumi_count\\n" > ${prefix}.sierra_quant.tsv
    printf "${meta.library_id}\\t${group_level}\\t${group_id}\\tstub_site_001\\tSTUBBARC001\\t1\\n" >> ${prefix}.sierra_quant.tsv
    printf "field\\tvalue\\nmode\\tstub\\n" > ${prefix}.sierra_quant.manifest.tsv
    printf "SIERRA_QUANT stub executed\\n" > ${prefix}.sierra_quant.log
    """
}
