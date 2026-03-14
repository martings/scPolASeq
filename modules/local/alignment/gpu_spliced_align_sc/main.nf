process GPU_SPLICED_ALIGN_SC {
    tag "$meta.library_id"
    label 'process_high'

    input:
    tuple val(meta), path(reads), path(whitelist), path(star_index), path(gtf), val(protocol_mode)

    output:
    tuple val(meta), path("${meta.library_id}.gpu_placeholder.bam"), path("${meta.library_id}.gpu_placeholder.bam.bai"), emit: bam_bundle

    script:
    error "GPU_SPLICED_ALIGN_SC is a reserved extension point and is not implemented in the MVP scaffold."

    stub:
    """
    touch ${meta.library_id}.gpu_placeholder.bam
    touch ${meta.library_id}.gpu_placeholder.bam.bai
    """
}
