process PREPARE_REFERENCE_BUNDLE {
    tag "prepare_reference_bundle"
    label 'process_medium'
    label 'process_python'

    conda "${projectDir}/envs/alignment.yml"
    publishDir "${params.outdir}/reference", mode: params.publish_dir_mode

    input:
    path genome_fasta
    path gtf
    path known_polya
    path priming_blacklist

    output:
    tuple val("primary"), path("star_index"), path("reference.annotation.gtf"), path("reference.genome.fa"), path("reference.chrom.sizes"), path("known_polya.reference.tsv"), path("priming_blacklist.reference.bed"), emit: reference_meta
    path "reference_manifest.json", emit: manifest
    path "reference.annotation.gtf", emit: gtf
    path "reference.genome.fa", emit: fasta

    script:
    """
    python ${projectDir}/bin/prepare_reference_bundle.py \\
        --genome-fasta ${genome_fasta} \\
        --gtf ${gtf} \\
        --known-polya ${known_polya} \\
        --priming-blacklist ${priming_blacklist} \\
        --out-fasta reference.genome.fa \\
        --out-gtf reference.annotation.gtf \\
        --out-fai reference.genome.fa.fai \\
        --out-dict reference.genome.dict \\
        --out-sizes reference.chrom.sizes \\
        --out-known-polya known_polya.reference.tsv \\
        --out-blacklist priming_blacklist.reference.bed \\
        --out-manifest reference_manifest.json

    mkdir -p star_index
    if command -v STAR >/dev/null 2>&1; then
        STAR \\
            --runMode genomeGenerate \\
            --runThreadN ${task.cpus} \\
            --genomeDir star_index \\
            --genomeFastaFiles reference.genome.fa \\
            --sjdbGTFfile reference.annotation.gtf || true
    fi
    touch star_index/SA
    touch star_index/Genome
    """

    stub:
    """
    cp ${genome_fasta} reference.genome.fa
    cp ${gtf} reference.annotation.gtf
    touch reference.genome.fa.fai
    touch reference.genome.dict
    touch reference.chrom.sizes
    touch known_polya.reference.tsv
    touch priming_blacklist.reference.bed
    echo "{}" > reference_manifest.json
    mkdir -p star_index
    touch star_index/SA
    touch star_index/Genome
    """
}
