process SCAPTURE_ANNOTATION {
    tag "${meta.library_id}"
    label 'process_high'

    conda "${projectDir}/envs/scapture.yml"
    container "ghcr.io/martings/scapture:f056f02"

    input:
    tuple val(meta), path(scapture_gtf), path(genome_fasta), path(genome_fai), path(chrom_sizes)

    output:
    tuple val(meta), path("${prefix}_annot.*"), emit: annotation
    tuple val(meta), path("${prefix}.annotation.log"), emit: log
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--extend 2000'
    prefix   = task.ext.prefix ?: "${meta.library_id}"
    """
    set -euo pipefail
    export TMPDIR=\$(pwd)/tmp
    mkdir -p "\$TMPDIR"
    test -s ${genome_fasta}
    test -s ${genome_fai}
    test -s ${chrom_sizes}
    # SCAPTURE uses the FASTA index via the conventional adjacent filename.
    # Nextflow stages the source as a symlink; use a local link rather than write
    # beside the source FASTA in a shared reference directory.
    ln -sfn "\$(readlink -f ${genome_fasta})" ${prefix}.scapture.genome.fa
    ln -sfn "\$(readlink -f ${genome_fai})" ${prefix}.scapture.genome.fa.fai
    scapture -m annotation \
        -o ${prefix}_annot \
        -g ${prefix}.scapture.genome.fa \
        --gtf ${scapture_gtf} \
        --cs ${chrom_sizes} \
        ${args} \
        &> ${prefix}.annotation.log
    test -s ${prefix}_annot.genetype.gtf
    test -s ${prefix}_annot.genetype.bed
    commit=\$(git -C /opt/scapture rev-parse HEAD 2>/dev/null || printf 'unknown')
    printf '"%s":\n  scapture_commit: "%s"\n  scapture_script_sha256: "%s"\n' \
        "${task.process}" "\$commit" "\$(sha256sum "\$(command -v scapture)" | cut -d' ' -f1)" > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.library_id}"
    """
    touch ${prefix}_annot.genetype.gtf ${prefix}_annot.genetype.bed
    printf 'stub scapture -m annotation\n' > ${prefix}.annotation.log
    printf '"%s":\n  scapture_commit: "f056f0202e64f801c52f885c069b3e31660e02f1"\n' "${task.process}" > versions.yml
    """
}
