process FETCH_POLYA_REFERENCES {
    tag "fetch_polya_references"
    label 'process_low'
    label 'process_python'

    conda "${projectDir}/envs/alignment.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-alignment.sif" : null
    publishDir "${params.outdir}/reference", mode: params.publish_dir_mode

    input:
    path known_polya      // legacy fallback (may be NO_FILE sentinel)
    path polya_db_file    // staged local PolyA_DB file (may be NO_FILE sentinel)
    path polyasite_file   // staged local PolyASite file (may be NO_FILE sentinel)

    output:
    path "polya_references.tsv",          emit: polya
    path "polya_references.manifest.tsv", emit: manifest

    script:
    def db_spec    = params.polya_db  ?: "skip"
    def site_spec  = params.polyasite ?: "skip"
    def merge_dist = params.pas_reference_merge_distance ?: 25
    def cache_dir  = params.pas_reference_cache_dir ?: "${params.outdir}/reference/polya_cache"
    """
    python ${projectDir}/bin/merge_polya_references.py \\
        --known-polya    ${known_polya}    \\
        --polya-db       "${db_spec}"      \\
        --polya-db-file  ${polya_db_file}  \\
        --polyasite      "${site_spec}"    \\
        --polyasite-file ${polyasite_file} \\
        --cache-dir      ${cache_dir}      \\
        --merge-distance ${merge_dist}     \\
        --out-tsv        polya_references.tsv                   \\
        --out-manifest   polya_references.manifest.tsv
    """

    stub:
    """
    printf "gene_id\tgene_name\tchrom\tstart\tend\tstrand\tscore\tsource\n" > polya_references.tsv
    printf "field\tvalue\nmode\tstub\n" > polya_references.manifest.tsv
    """
}
