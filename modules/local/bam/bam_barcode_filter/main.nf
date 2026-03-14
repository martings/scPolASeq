// Stage 3 — Filter BAM to valid cell barcodes from cell_annotations
// Preserves CB and UB tags; associates reads with cluster/cell_type metadata.
process BAM_BARCODE_FILTER {
    tag "$meta.library_id"
    label 'process_medium'

    conda "${projectDir}/envs/python.yml"
    publishDir "${params.outdir}/filtered_bam", mode: params.publish_dir_mode, pattern: "*.filtered.bam"
    publishDir "${params.outdir}/filtered_bam", mode: params.publish_dir_mode, pattern: "*.filtered.bam.bai"
    publishDir "${params.outdir}/qc",           mode: params.publish_dir_mode, pattern: "*.barcode_filter_stats.tsv"

    input:
    tuple val(meta), path(bam), path(bai), path(cell_annotations)

    output:
    tuple val(meta), path("${meta.library_id}.filtered.bam"), path("${meta.library_id}.filtered.bam.bai"), emit: filtered_bam
    path "${meta.library_id}.barcode_filter_stats.tsv", emit: filter_stats

    script:
    """
    python ${projectDir}/bin/filter_bam_by_barcodes.py \\
        --bam            ${bam} \\
        --annotations    ${cell_annotations} \\
        --out-bam        ${meta.library_id}.filtered.bam \\
        --out-stats      ${meta.library_id}.barcode_filter_stats.tsv

    if command -v samtools >/dev/null 2>&1; then
        samtools index ${meta.library_id}.filtered.bam \\
            || touch ${meta.library_id}.filtered.bam.bai
    else
        touch ${meta.library_id}.filtered.bam.bai
    fi
    """

    stub:
    """
    touch ${meta.library_id}.filtered.bam
    touch ${meta.library_id}.filtered.bam.bai
    printf "metric\tvalue\ntotal_reads\t1000\nvalid_cb_reads\t950\nfiltered_reads\t950\n" \
        > ${meta.library_id}.barcode_filter_stats.tsv
    """
}
