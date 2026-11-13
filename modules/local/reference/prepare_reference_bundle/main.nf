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
    path prebuilt_star_index, stageAs: 'prebuilt_star_index_in'   // directory or NO_FILE sentinel
    path chrom_sizes                                             // optional chrom.sizes or NO_FILE sentinel

    output:
    tuple val("primary"), path("star_index"), path("reference.annotation.gtf"), path("reference.genome.fa"), path("reference.genome.fa.fai"), path("reference.chrom.sizes"), path("known_polya.reference.tsv"), path("priming_blacklist.reference.bed"), emit: reference_meta
    path "reference_manifest.json", emit: manifest
    path "reference.annotation.gtf", emit: gtf
    path "reference.genome.fa", emit: fasta
    path "reference.genome.fa.fai", emit: fasta_fai

    script:
    def use_prebuilt = !prebuilt_star_index.name.startsWith('NO_FILE')
    def staged_prebuilt = use_prebuilt ? 'prebuilt_star_index_in' : prebuilt_star_index
    def ref_root_arg = params.reference_root ? "--reference-root \"${params.reference_root}\"" : ''
    def ref_build_arg = params.reference_build ? "--reference-build \"${params.reference_build}\"" : ''
    def star_mode = use_prebuilt ? 'prebuilt' : 'generated'
    def bundle_script = file("${moduleDir}/../../../../bin/prepare_reference_bundle.py")
    """
    STAR_BIN_VERSION="unknown"
    if ${use_prebuilt} && command -v STAR >/dev/null 2>&1; then
        STAR_BIN_VERSION="\$(STAR --version 2>/dev/null | head -n 1 || echo unknown)"
    fi

    python3 ${bundle_script} \\
        --genome-fasta ${genome_fasta} \\
        --gtf ${gtf} \\
        --known-polya ${known_polya} \\
        --priming-blacklist ${priming_blacklist} \\
        --chrom-sizes ${chrom_sizes} \\
        --prebuilt-star-index ${staged_prebuilt} \\
        --star-binary-version "\${STAR_BIN_VERSION}" \\
        ${ref_root_arg} \\
        ${ref_build_arg} \\
        --out-fasta reference.genome.fa \\
        --out-gtf reference.annotation.gtf \\
        --out-fai reference.genome.fa.fai \\
        --out-dict reference.genome.dict \\
        --out-sizes reference.chrom.sizes \\
        --out-known-polya known_polya.reference.tsv \\
        --out-blacklist priming_blacklist.reference.bed \\
        --out-manifest reference_manifest.json

    if ${use_prebuilt}; then
        mkdir -p star_index
        # Try hard-linking first (zero disk space, instantaneous, and produces
        # content-identical files so downstream cache hashes remain stable).
        # Hard links fail when the source and destination are on different
        # filesystems, so fall back to a regular copy in that case.
        if cp -rl ${staged_prebuilt}/. star_index/ 2>/dev/null; then
            echo "Prebuilt STAR index staged via hard-links."
        else
            echo "Hard-link staging failed (likely cross-device); falling back to regular copy." >&2
            cp -r ${staged_prebuilt}/. star_index/ || \
                { echo "ERROR: Failed to stage prebuilt STAR index." >&2; exit 1; }
        fi
    else
        mkdir -p star_index
        command -v STAR >/dev/null 2>&1 || { echo 'STAR executable unavailable' >&2; exit 127; }
        STAR \\
            --runMode genomeGenerate \\
            --runThreadN ${task.cpus} \\
            --genomeDir star_index \\
            --genomeFastaFiles reference.genome.fa \\
            --sjdbGTFfile reference.annotation.gtf ||
            { echo 'STAR index generation failed' >&2; exit 1; }
        STAR_BIN_VERSION="\$(STAR --version 2>/dev/null | head -n 1 || echo unknown)"
    fi
    test -s star_index/SA && test -s star_index/Genome ||
        { echo 'STAR index is incomplete' >&2; exit 1; }

    python3 ${bundle_script} \\
        --finalize-star-index star_index \\
        --star-index-mode ${star_mode} \\
        --star-index-source ${prebuilt_star_index} \\
        --star-binary-version "\${STAR_BIN_VERSION}" \\
        --out-manifest reference_manifest.json
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
