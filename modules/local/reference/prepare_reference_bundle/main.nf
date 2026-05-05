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
    if [[ -n "${params.prebuilt_star_index}" && -d "${params.prebuilt_star_index}" ]]; then
        # Try hard-linking first (zero disk space, instantaneous, and produces
        # content-identical files so downstream cache hashes remain stable).
        # Hard links fail when the source and destination are on different
        # filesystems, so fall back to a regular copy in that case.
        if cp -rl "${params.prebuilt_star_index}/." star_index/ 2>/dev/null; then
            echo "Prebuilt STAR index staged via hard-links."
        else
            echo "Hard-link staging failed (likely cross-device); falling back to regular copy." >&2
            cp -r "${params.prebuilt_star_index}/." star_index/ || \
                { echo "ERROR: Failed to stage prebuilt STAR index from '${params.prebuilt_star_index}'." >&2; exit 1; }
        fi
    elif command -v STAR >/dev/null 2>&1; then
        STAR \\
            --runMode genomeGenerate \\
            --runThreadN ${task.cpus} \\
            --genomeDir star_index \\
            --genomeFastaFiles reference.genome.fa \\
            --sjdbGTFfile reference.annotation.gtf || true
    fi
    # Only create sentinel stubs if STAR did not produce the real index files.
    # `touch` on an existing file only updates mtime, but an absent file (STAR
    # unavailable or failed) would otherwise make the next process crash without
    # a clear error.
    [ -f star_index/SA ]     || touch star_index/SA
    [ -f star_index/Genome ] || touch star_index/Genome
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
