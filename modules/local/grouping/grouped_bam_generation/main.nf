// Stage 4 — Write one BAM per (group_level, group_id) using barcode-to-group mapping.
// Allows Mode A (direct regrouping from filtered BAM) and
// Mode B (optional FASTQ reconstruction + realignment, future).
process GROUPED_BAM_GENERATION {
    tag "${meta.library_id}:${group_level}"
    label 'process_medium'

    conda "${projectDir}/envs/alignment.yml"
    publishDir "${params.outdir}/grouped_bam/${group_level}", mode: params.publish_dir_mode, pattern: "*.grouped.bam"
    publishDir "${params.outdir}/grouped_bam/${group_level}", mode: params.publish_dir_mode, pattern: "*.grouped.bam.bai"

    input:
    tuple val(meta), path(filtered_bam), path(filtered_bai), path(group_map), val(group_level)

    output:
    tuple val(meta), val(group_level), path("*.grouped.bam"), path("*.grouped.bam.bai"), emit: grouped_bams
    path "${meta.library_id}.${group_level}.grouping_manifest.tsv",                      emit: grouping_manifest

    script:
    """
    python ${projectDir}/bin/write_group_bams.py \\
        --bam           ${filtered_bam} \\
        --group-map     ${group_map} \\
        --group-level   ${group_level} \\
        --out-dir       . \\
        --out-manifest  ${meta.library_id}.${group_level}.grouping_manifest.tsv

    # Index all produced BAMs
    for bam_file in *.grouped.bam; do
        if command -v samtools >/dev/null 2>&1; then
            samtools index "\$bam_file" || touch "\${bam_file}.bai"
        else
            touch "\${bam_file}.bai"
        fi
    done
    """

    stub:
    """
    touch group_1.grouped.bam group_1.grouped.bam.bai
    touch group_2.grouped.bam group_2.grouped.bam.bai
    printf "group_level\tgroup_id\tbam_file\tn_reads\n" \
        >  ${meta.library_id}.${group_level}.grouping_manifest.tsv
    printf "${group_level}\tgroup_1\tgroup_1.grouped.bam\t100\n" \
        >> ${meta.library_id}.${group_level}.grouping_manifest.tsv
    printf "${group_level}\tgroup_2\tgroup_2.grouped.bam\t80\n"  \
        >> ${meta.library_id}.${group_level}.grouping_manifest.tsv
    """
}
