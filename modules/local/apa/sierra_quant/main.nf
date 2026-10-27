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
 *     path "<library>.<group_level>.<group_id>.versions.yml" (Sierra and R runtime versions)
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
    tuple val(meta), val(group_level), val(group_id), path("${meta.library_id}.${group_level}.${group_id.toString().replaceAll(/[^A-Za-z0-9._-]+/, '_')}.sierra_quant.tsv"), emit: quantified_pas
    path "${meta.library_id}.${group_level}.${group_id.toString().replaceAll(/[^A-Za-z0-9._-]+/, '_')}.sierra_quant.manifest.tsv", emit: manifest
    path "${meta.library_id}.${group_level}.${group_id.toString().replaceAll(/[^A-Za-z0-9._-]+/, '_')}.sierra_quant.log", emit: log
    path "${meta.library_id}.${group_level}.${group_id.toString().replaceAll(/[^A-Za-z0-9._-]+/, '_')}.versions.yml", emit: versions

    script:
    def safeGroupId = group_id.toString().replaceAll(/[^A-Za-z0-9._-]+/, '_')
    def prefix = "${meta.library_id}.${group_level}.${safeGroupId}"
    def quotedGroupId = "'" + group_id.toString().replace("'", "'\"'\"'") + "'"
    def quotedLibraryId = "'" + meta.library_id.toString().replace("'", "'\"'\"'") + "'"
    def quotedGroupLevel = "'" + group_level.toString().replace("'", "'\"'\"'") + "'"
    def ncores = task.cpus ?: 1
    """
    test -s ${grouped_bam} && test -s ${grouped_bai} || { echo 'Sierra BAM/index missing or empty' >&2; exit 1; }
    test -s ${pas_reference} && test -s ${cell_annotations} && test -s ${gtf} || { echo 'Sierra reference/annotations missing or empty' >&2; exit 1; }
    Rscript ${quant_script} \\
        --bam               ${grouped_bam} \\
        --gtf               ${gtf} \\
        --pas-reference     ${pas_reference} \\
        --cell-annotations  ${cell_annotations} \\
        --library-id        ${quotedLibraryId} \\
        --group-level       ${quotedGroupLevel} \\
        --group-id          ${quotedGroupId} \\
        --output            ${prefix}.sierra_quant.tsv \\
        --log               ${prefix}.sierra_quant.log \\
        --versions          ${prefix}.versions.yml \\
        --ncores            ${ncores}

    printf "field\\tvalue\\n"                                    > ${prefix}.sierra_quant.manifest.tsv
    printf "grouped_bam\\t${grouped_bam.name}\\n"              >> ${prefix}.sierra_quant.manifest.tsv
    printf "grouped_bai\\t${grouped_bai.name}\\n"              >> ${prefix}.sierra_quant.manifest.tsv
    printf "gtf\\t${gtf.name}\\n"                              >> ${prefix}.sierra_quant.manifest.tsv
    printf "pas_reference\\t${pas_reference.name}\\n"          >> ${prefix}.sierra_quant.manifest.tsv
    printf "cell_annotations\\t${cell_annotations.name}\\n"    >> ${prefix}.sierra_quant.manifest.tsv
    printf "quant_source\\tSierra::CountPeaks\\n"              >> ${prefix}.sierra_quant.manifest.tsv
    printf 'group_id\\t%s\\n' ${quotedGroupId}                >> ${prefix}.sierra_quant.manifest.tsv
    test -s ${prefix}.sierra_quant.tsv && test -s ${prefix}.versions.yml || { echo 'Sierra outputs missing or empty' >&2; exit 1; }
    """

    stub:
    def safeGroupId = group_id.toString().replaceAll(/[^A-Za-z0-9._-]+/, '_')
    def prefix = "${meta.library_id}.${group_level}.${safeGroupId}"
    def quotedGroupId = "'" + group_id.toString().replace("'", "'\"'\"'") + "'"
    def quotedLibraryId = "'" + meta.library_id.toString().replace("'", "'\"'\"'") + "'"
    def quotedGroupLevel = "'" + group_level.toString().replace("'", "'\"'\"'") + "'"
    """
    printf "library_id\\tgroup_level\\tgroup_id\\tsite_id\\tcell_barcode\\tumi_count\\n" > ${prefix}.sierra_quant.tsv
    printf '%s\\t%s\\t%s\\tstub_site_001\\tSTUBBARC001\\t1\\n' ${quotedLibraryId} ${quotedGroupLevel} ${quotedGroupId} >> ${prefix}.sierra_quant.tsv
    printf "field\\tvalue\\nmode\\tstub\\n" > ${prefix}.sierra_quant.manifest.tsv
    printf "SIERRA_QUANT stub executed\\n" > ${prefix}.sierra_quant.log
    printf "SIERRA_QUANT:\\n  Sierra: stub\\n  R: stub\\n" > ${prefix}.versions.yml
    """
}
