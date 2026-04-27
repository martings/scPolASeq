// Merge one or more per-library filtered BAMs into a single sample-level BAM.
// When the sample has only one library the BAM is passed through without copying.
process BAM_MERGE_SAMPLE {
    tag "${meta.library_id}"
    label 'process_medium'

    conda "${projectDir}/envs/alignment.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-alignment.sif" : null

    input:
    tuple val(meta), path(bams), path(bais)

    output:
    tuple val(meta), path("${meta.library_id}.merged.bam"), path("${meta.library_id}.merged.bam.bai"), emit: merged_bam

    script:
    def prefix  = meta.library_id
    def bam_list = (bams instanceof List ? bams : [bams])
    def n       = bam_list.size()
    """
    if [ "${n}" -eq 1 ]; then
        ln -sf \$(realpath ${bam_list[0]}) ${prefix}.merged.bam
        samtools index ${prefix}.merged.bam
    else
        samtools merge -@ ${task.cpus} -f ${prefix}.merged.bam ${bam_list.join(' ')}
        samtools index -@ ${task.cpus} ${prefix}.merged.bam
    fi
    """

    stub:
    """
    touch ${meta.library_id}.merged.bam ${meta.library_id}.merged.bam.bai
    """
}
