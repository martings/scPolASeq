process PREPARE_REFERENCE_BUNDLE {
    tag "prepare_reference_bundle"
    label 'process_medium'
    label 'process_python'

    conda "${projectDir}/envs/alignment.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-alignment.sif" : null
    publishDir "${params.outdir}/reference", mode: params.publish_dir_mode

    input:
    path genome_fasta
    path gtf
    path known_polya
    path priming_blacklist
    path prebuilt_star_index   // directory or NO_FILE sentinel

    output:
    tuple val("primary"), path("star_index"), path("reference.annotation.gtf"), path("reference.genome.fa"), path("reference.genome.fa.fai"), path("reference.chrom.sizes"), path("known_polya.reference.tsv"), path("priming_blacklist.reference.bed"), emit: reference_meta
    path "reference_manifest.json", emit: manifest
    path "reference.annotation.gtf", emit: gtf
    path "reference.genome.fa", emit: fasta
    path "reference.genome.fa.fai", emit: fasta_fai

    script:
    def use_prebuilt = prebuilt_star_index.name != 'NO_FILE'
    """
    python3 ${projectDir}/bin/prepare_reference_bundle.py \\
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

    if ${use_prebuilt}; then
        cp -rL star_index .star_index_staged
        rm -rf star_index
        mv .star_index_staged star_index
    else
        mkdir -p star_index
        if command -v STAR >/dev/null 2>&1; then
            STAR \\
                --runMode genomeGenerate \\
                --runThreadN ${task.cpus} \\
                --genomeDir star_index \\
                --genomeFastaFiles reference.genome.fa \\
                --sjdbGTFfile reference.annotation.gtf || true
        fi
        [ -f star_index/SA ]     || touch star_index/SA
        [ -f star_index/Genome ] || touch star_index/Genome
    fi
    """

    stub:
    """
    touch reference.genome.fa
    touch reference.annotation.gtf
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
