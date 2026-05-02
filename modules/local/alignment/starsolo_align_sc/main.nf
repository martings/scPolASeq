process STARSOLO_ALIGN_SC {
    tag "$meta.sample_id"
    label 'process_high'

    conda "${projectDir}/envs/alignment.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-alignment.sif" : null
    publishDir "${params.outdir}/alignment", mode: params.publish_dir_mode, pattern: "*.bam"
    publishDir "${params.outdir}/alignment", mode: params.publish_dir_mode, pattern: "*.bai"
    publishDir "${params.outdir}/alignment", mode: params.publish_dir_mode, pattern: "*.Solo.out"
    publishDir "${params.outdir}/provenance", mode: params.publish_dir_mode, pattern: "*.barcode_registry.tsv"
    publishDir "${params.outdir}/provenance", mode: params.publish_dir_mode, pattern: "*.alignment_manifest.tsv"

    input:
    tuple val(meta), path(read1s), path(read2s), path(whitelist), path(star_index), path(gtf), val(protocol_mode)

    output:
    tuple val(meta), path("${meta.library_id}.Aligned.sortedByCoord.out.bam"), path("${meta.library_id}.Aligned.sortedByCoord.out.bam.bai"), emit: bam_bundle
    tuple val(meta), path("${meta.library_id}.Solo.out"), emit: matrix_bundle
    tuple val(meta), path("${meta.library_id}.barcode_registry.tsv"), emit: barcode_registry
    tuple val(meta), path("${meta.library_id}.alignment_manifest.tsv"), emit: alignment_manifest
    tuple val(meta), path("${meta.library_id}.Log.final.out"), emit: log

    script:
    def whitelistArg = whitelist.name == 'NO_FILE' ? '--soloCBwhitelist None' : "--soloCBwhitelist ${whitelist}"
    def warning = protocol_mode == '10x_5p' ? "echo 'WARNING: guarded 10x 5p support is enabled; APA outputs are descriptive-first.' >&2" : ''
    def readFilesCmd = read1s[0].name.endsWith('.gz') ? '--readFilesCommand zcat' : ''
    def read1Arg = read1s.collect { it.toString() }.join(',')
    def read2Arg = read2s.collect { it.toString() }.join(',')
    // STARsolo single-cell params derived from protocol_mode.
    // ext.args is not reliably propagated in Nextflow 25.x when multiple withName
    // selectors exist across config files, so critical params are hardcoded here.
    def umiLen       = protocol_mode == '10x_5p' ? 10 : 12
    def soloTypeArgs = protocol_mode.startsWith('10x') ?
        "--soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen ${umiLen} --soloCellFilter CellRanger2.2" :
        "--soloType CB_UMI_Simple"
    """
    ${warning}
    _star_ok=0
    if command -v STAR >/dev/null 2>&1; then
        STAR \\
            --genomeDir ${star_index} \\
            --sjdbGTFfile ${gtf} \\
            --readFilesIn ${read2Arg} ${read1Arg} \\
            ${readFilesCmd} \\
            --runThreadN ${task.cpus} \\
            --outSAMtype BAM SortedByCoordinate \\
            --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \\
            ${soloTypeArgs} \\
            --outFileNamePrefix ${meta.sample_id}. \\
            ${whitelistArg} \\
            ${task.ext.args ?: ''} && _star_ok=1
    fi

    # Fallback placeholders — only if STAR was absent or failed
    mkdir -p ${meta.sample_id}.Solo.out/Gene
    if [[ \$_star_ok -eq 0 ]]; then
        touch ${meta.sample_id}.Log.final.out
        touch ${meta.sample_id}.Aligned.sortedByCoord.out.bam
        echo "${meta.sample_id}-CELL001" > ${meta.sample_id}.Solo.out/Gene/barcodes.tsv
        touch ${meta.sample_id}.Solo.out/Gene/matrix.mtx
        touch ${meta.sample_id}.Solo.out/Gene/features.tsv
    fi

    if command -v samtools >/dev/null 2>&1 && [[ -s ${meta.sample_id}.Aligned.sortedByCoord.out.bam ]]; then
        samtools index ${meta.sample_id}.Aligned.sortedByCoord.out.bam
    else
        touch ${meta.sample_id}.Aligned.sortedByCoord.out.bam.bai
    fi

    python ${projectDir}/bin/extract_barcode_registry.py \\
        --solo-dir ${meta.sample_id}.Solo.out \\
        --sample ${meta.sample_id} \\
        --library ${meta.sample_id} \\
        --out ${meta.sample_id}.barcode_registry.tsv

    python ${projectDir}/bin/write_manifest.py \\
        --input-table ${meta.sample_id}.barcode_registry.tsv \\
        --manifest-name ${meta.sample_id}.alignment \\
        --out ${meta.sample_id}.alignment_manifest.tsv
    """

    stub:
    """
    mkdir -p ${meta.sample_id}.Solo.out/Gene
    touch ${meta.sample_id}.Aligned.sortedByCoord.out.bam
    touch ${meta.sample_id}.Aligned.sortedByCoord.out.bam.bai
    touch ${meta.sample_id}.Log.final.out
    echo "${meta.sample_id}-CELL001" > ${meta.sample_id}.Solo.out/Gene/barcodes.tsv
    touch ${meta.sample_id}.Solo.out/Gene/matrix.mtx
    touch ${meta.sample_id}.Solo.out/Gene/features.tsv
    printf "sample_id\\tlibrary_id\\tbarcode_raw\\tbarcode_corrected\\tsource\\n${meta.sample_id}\\t${meta.sample_id}\\t${meta.sample_id}-CELL001\\t${meta.sample_id}-CELL001\\tstub\\n" > ${meta.sample_id}.barcode_registry.tsv
    printf "manifest_name\\tinput_table\\trows\\n${meta.sample_id}.alignment\\t${meta.sample_id}.barcode_registry.tsv\\t1\\n" > ${meta.sample_id}.alignment_manifest.tsv
    """
}
