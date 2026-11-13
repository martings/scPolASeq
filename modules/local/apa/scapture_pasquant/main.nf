process SCAPTURE_PASQUANT {
    tag "${meta.library_id}"
    label 'process_high'

    conda "${projectDir}/envs/scapture.yml"
    container "ghcr.io/martings/scapture:f056f02"

    input:
    tuple val(meta), path(bam), path(bai), path(pas_reference), path(cell_annotations), path(filter_helper, stageAs: 'FilterBamByTag.source')

    output:
    tuple val(meta), path("${prefix}.KeepCell.UMIs.tsv.gz"), emit: matrix
    tuple val(meta), path("${prefix}.PASquant.metrics.tsv"), emit: metrics
    tuple val(meta), path("${prefix}.PASquant.manifest.tsv"), emit: manifest
    tuple val(meta), path("${prefix}.PASquant.log"), emit: log
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix     = task.ext.prefix ?: "${meta.library_id?.toString()}"
    def sample = (meta.sample_id ?: meta.library_id)?.toString()
    if (!(prefix ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
        error "SCAPTURE_PASQUANT requires a safe library_id/prefix: ${prefix}"
    }
    if (!(sample ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
        error "SCAPTURE_PASQUANT requires a safe sample_id: ${sample}"
    }
    def ncores = task.cpus ?: 1
    def args   = task.ext.args != null ? task.ext.args.toString().trim() : '--strand forward'
    """
    set -euo pipefail
    export TMPDIR=\$(pwd)/tmp
    export MPLCONFIGDIR=\$(pwd)/tmp/mpl
    export XDG_CACHE_HOME=\$(pwd)/tmp/cache
    export SCAPTUREPATH=/opt/scapture/
    mkdir -p "\$TMPDIR" "\$MPLCONFIGDIR" "\$XDG_CACHE_HOME"
    test -s ${bam}
    test -s ${bai}
    test -s ${pas_reference}
    test -s ${cell_annotations}

    awk -v sid='${sample}' 'BEGIN{FS="\t"}
        NR==1 {
            for(i=1;i<=NF;i++) {
                if(\$i=="barcode_corrected") bc=i
                else if(\$i=="barcode_raw" && !bc) bc=i
                if(\$i=="sample_id") si=i
            }
            if (!bc) { print "ERROR: no barcode column" > "/dev/stderr"; exit 1 }
            next
        }
        si && \$si != sid { next }
        { print \$bc }
    ' ${cell_annotations} | sort -u > barcodes.txt
    test -s barcodes.txt || { echo 'ERROR: no barcodes for ${sample}' >&2; exit 1; }

    cp ${filter_helper} ./FilterBamByTag
    chmod +x ./FilterBamByTag
    export PATH=\$(pwd):\$PATH
    scapture -m PASquant \
        -b ${bam} \
        --celllist barcodes.txt \
        --pas ${pas_reference} \
        -o ${prefix} \
        -p ${ncores} \
        ${args} \
        &> ${prefix}.PASquant.log

    matrix=${prefix}.KeepCell.UMIs.tsv.gz
    test -s "\$matrix" || { echo 'ERROR: PASquant produced no matrix' >&2; exit 1; }
    gzip -t "\$matrix"
    zcat "\$matrix" | awk -F '\t' '
        NR==1 { width=NF; if (\$1!="gene" || width<2) exit 1; next }
        NF!=width || \$1=="" { exit 1 }
        END { if (NR<2) exit 1 }
    ' || { echo 'ERROR: malformed or empty PASquant matrix' >&2; exit 1; }
    printf 'metric\tvalue\nbarcodes\t%s\n' "\$(wc -l < barcodes.txt)" > ${prefix}.PASquant.metrics.tsv
    zcat "\$matrix" | awk -F '\t' 'BEGIN { OFS=FS } END { print "matrix_site_rows", NR-1; print "matrix_cell_columns", NF-1 }' >> ${prefix}.PASquant.metrics.tsv
    printf 'field\tvalue\nlibrary_id\t%s\nstatus\tcomplete\nreference_bed\t%s\nquant_matrix\t%s\n' \
        '${prefix}' '${pas_reference.name}' "\$matrix" > ${prefix}.PASquant.manifest.tsv
    awk 'NR>1 { print \$1 "\t" \$2 }' ${prefix}.PASquant.metrics.tsv >> ${prefix}.PASquant.manifest.tsv
    commit=\$(git -C /opt/scapture rev-parse HEAD 2>/dev/null || printf unknown)
    printf '"%s":\n  scapture_commit: "%s"\n  scapture_script_sha256: "%s"\n  filter_helper_sha256: "%s"\n' \
        "${task.process}" "\$commit" "\$(sha256sum "\$(command -v scapture)" | cut -d' ' -f1)" \
        "\$(sha256sum FilterBamByTag | cut -d' ' -f1)" > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.library_id?.toString()}"
    """
    printf 'gene\tAAACCCAAGAAACACT-1\npas_1\t5\n' | gzip -c > ${prefix}.KeepCell.UMIs.tsv.gz
    printf 'metric\tvalue\nbarcodes\t1\nmatrix_site_rows\t1\nmatrix_cell_columns\t1\n' > ${prefix}.PASquant.metrics.tsv
    printf 'field\tvalue\nlibrary_id\t%s\nstatus\tcomplete\nreference_bed\tstub.bed\nquant_matrix\t%s.KeepCell.UMIs.tsv.gz\nbarcodes\t1\nmatrix_site_rows\t1\nmatrix_cell_columns\t1\n' \
        '${prefix}' '${prefix}' > ${prefix}.PASquant.manifest.tsv
    printf 'stub scapture -m PASquant\n' > ${prefix}.PASquant.log
    printf '"%s":\n  scapture_commit: "f056f0202e64f801c52f885c069b3e31660e02f1"\n' "${task.process}" > versions.yml
    """
}
