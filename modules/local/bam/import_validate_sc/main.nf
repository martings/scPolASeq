process IMPORT_VALIDATE_SC {
    tag "$meta.library_id"
    label 'process_medium'

    conda "${projectDir}/envs/alignment.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-alignment.sif" : null
    publishDir "${params.outdir}/alignment", mode: params.publish_dir_mode, pattern: "*.bam"
    publishDir "${params.outdir}/alignment", mode: params.publish_dir_mode, pattern: "*.bai"
    publishDir "${params.outdir}/provenance", mode: params.publish_dir_mode, pattern: "*.barcode_registry.tsv"
    publishDir "${params.outdir}/provenance", mode: params.publish_dir_mode, pattern: "*.alignment_manifest.tsv"

    input:
    tuple val(meta), path(bam), path(matrix_bundle), path(whitelist)

    output:
    tuple val(meta), path("${meta.library_id}.Aligned.sortedByCoord.out.bam"), path("${meta.library_id}.Aligned.sortedByCoord.out.bam.bai"), emit: bam_bundle
    tuple val(meta), path("${meta.library_id}.matrix_bundle"), emit: matrix_bundle
    tuple val(meta), path("${meta.library_id}.barcode_registry.tsv"), emit: barcode_registry
    tuple val(meta), path("${meta.library_id}.alignment_manifest.tsv"), emit: alignment_manifest

    script:
    """
    cp ${bam} ${meta.library_id}.Aligned.sortedByCoord.out.bam
    mkdir -p ${meta.library_id}.matrix_bundle
    if [ "${matrix_bundle.getName()}" != "NO_FILE" ]; then
        cp -r ${matrix_bundle}/* ${meta.library_id}.matrix_bundle/ 2>/dev/null || true
    fi
    if command -v samtools >/dev/null 2>&1; then
        samtools index ${meta.library_id}.Aligned.sortedByCoord.out.bam || touch ${meta.library_id}.Aligned.sortedByCoord.out.bam.bai
    else
        touch ${meta.library_id}.Aligned.sortedByCoord.out.bam.bai
    fi
    python ${projectDir}/bin/extract_barcode_registry.py \\
        --bam ${meta.library_id}.Aligned.sortedByCoord.out.bam \\
        --sample ${meta.sample_id} \\
        --library ${meta.library_id} \\
        --out ${meta.library_id}.barcode_registry.tsv
    python ${projectDir}/bin/write_manifest.py \\
        --input-table ${meta.library_id}.barcode_registry.tsv \\
        --manifest-name ${meta.library_id}.alignment \\
        --out ${meta.library_id}.alignment_manifest.tsv
    """

    stub:
    """
    touch ${meta.library_id}.Aligned.sortedByCoord.out.bam
    touch ${meta.library_id}.Aligned.sortedByCoord.out.bam.bai
    mkdir -p ${meta.library_id}.matrix_bundle
    printf "sample_id\\tlibrary_id\\tbarcode_raw\\tbarcode_corrected\\tsource\\n${meta.sample_id}\\t${meta.library_id}\\t${meta.sample_id}-CELL001\\t${meta.sample_id}-CELL001\\tstub\\n" > ${meta.library_id}.barcode_registry.tsv
    printf "manifest_name\\tinput_table\\trows\\n${meta.library_id}.alignment\\t${meta.library_id}.barcode_registry.tsv\\t1\\n" > ${meta.library_id}.alignment_manifest.tsv
    """
}
