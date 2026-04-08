/*
 * Scaffold-only Sierra quantification module.
 *
 * Contract:
 *   input:
 *     tuple val(meta), val(group_level), val(group_id), path(grouped_bam), path(grouped_bai), path(pas_reference), path(cell_annotations)
 *   output:
 *     tuple val(meta), val(group_level), val(group_id), path("<library>.<group_level>.<group_id>.sierra_quant.tsv")
 *     path "<library>.<group_level>.<group_id>.sierra_quant.manifest.tsv"
 *     path "<library>.<group_level>.<group_id>.sierra_quant.log"
 */
process SIERRA_QUANT {
    tag "${meta.library_id}:${group_level}:${group_id}"
    label 'process_single'

    publishDir "${params.outdir}/sierra_quant/${group_level}", mode: params.publish_dir_mode

    input:
    tuple val(meta), val(group_level), val(group_id), path(grouped_bam), path(grouped_bai), path(pas_reference), path(cell_annotations)

    output:
    tuple val(meta), val(group_level), val(group_id), path("${meta.library_id}.${group_level}.${group_id}.sierra_quant.tsv"), emit: quantified_pas
    path "${meta.library_id}.${group_level}.${group_id}.sierra_quant.manifest.tsv", emit: manifest
    path "${meta.library_id}.${group_level}.${group_id}.sierra_quant.log", emit: log

    script:
    """
    printf "library_id\\tgroup_level\\tgroup_id\\tpas_reference\\tquant_source\\n" > ${meta.library_id}.${group_level}.${group_id}.sierra_quant.tsv
    printf "${meta.library_id}\\t${group_level}\\t${group_id}\\t${pas_reference.name}\\tscaffold_only\\n" >> ${meta.library_id}.${group_level}.${group_id}.sierra_quant.tsv

    printf "field\\tvalue\\n" > ${meta.library_id}.${group_level}.${group_id}.sierra_quant.manifest.tsv
    printf "grouped_bam\\t%s\\n" "${grouped_bam.name}" >> ${meta.library_id}.${group_level}.${group_id}.sierra_quant.manifest.tsv
    printf "grouped_bai\\t%s\\n" "${grouped_bai.name}" >> ${meta.library_id}.${group_level}.${group_id}.sierra_quant.manifest.tsv
    printf "pas_reference\\t%s\\n" "${pas_reference.name}" >> ${meta.library_id}.${group_level}.${group_id}.sierra_quant.manifest.tsv
    printf "cell_annotations\\t%s\\n" "${cell_annotations.name}" >> ${meta.library_id}.${group_level}.${group_id}.sierra_quant.manifest.tsv

    printf "SIERRA_QUANT scaffold executed\\n" > ${meta.library_id}.${group_level}.${group_id}.sierra_quant.log
    """

    stub:
    """
    printf "library_id\\tgroup_level\\tgroup_id\\tpas_reference\\tquant_source\\n" > ${meta.library_id}.${group_level}.${group_id}.sierra_quant.tsv
    printf "${meta.library_id}\\t${group_level}\\t${group_id}\\tstub\\tstub\\n" >> ${meta.library_id}.${group_level}.${group_id}.sierra_quant.tsv
    printf "field\\tvalue\\nmode\\tstub\\n" > ${meta.library_id}.${group_level}.${group_id}.sierra_quant.manifest.tsv
    printf "SIERRA_QUANT stub executed\\n" > ${meta.library_id}.${group_level}.${group_id}.sierra_quant.log
    """
}
