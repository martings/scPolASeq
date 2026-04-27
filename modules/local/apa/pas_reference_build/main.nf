/*
 * PAS reference builder.
 *
 * Contract:
 *   input:
 *     path site_catalog
 *     path terminal_exons
 *     path known_polya
 *   output:
 *     path "pas_reference.tsv"
 *     path "pas_reference_build.manifest.tsv"
 *     path "pas_reference_build.log"
 */
process PAS_REFERENCE_BUILD {
    tag "pas_reference_build"
    label 'process_single'
    label 'process_python'

    conda "${projectDir}/envs/python.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-python.sif" : null
    publishDir "${params.outdir}/pas_reference", mode: params.publish_dir_mode

    input:
    path site_catalog
    path terminal_exons
    path known_polya

    output:
    path "pas_reference.tsv",                emit: pas_reference
    path "pas_reference_build.manifest.tsv", emit: manifest
    path "pas_reference_build.log",          emit: log

    script:
    """
    python3 ${projectDir}/bin/build_pas_reference.py \\
        --site-catalog ${site_catalog} \\
        --terminal-exons ${terminal_exons} \\
        --known-polya ${known_polya} \\
        --out-tsv pas_reference.tsv \\
        --out-manifest pas_reference_build.manifest.tsv \\
        --log pas_reference_build.log \\
        --merge-distance ${params.pas_reference_merge_distance ?: 25}
    """

    stub:
    """
    printf "pas_reference_id\\tsite_id\\tgene_id\\tchrom\\tstart\\tend\\tstrand\\treference_source\\n" > pas_reference.tsv
    printf "pas_ref_stub_001\\tNA\\tNA\\tNA\\t0\\t0\\t+\\tstub\\n" >> pas_reference.tsv

    printf "field\\tvalue\\n" > pas_reference_build.manifest.tsv
    printf "site_catalog\\t${site_catalog.name}\\n"   >> pas_reference_build.manifest.tsv
    printf "terminal_exons\\t${terminal_exons.name}\\n" >> pas_reference_build.manifest.tsv
    printf "known_polya\\t${known_polya.name}\\n"     >> pas_reference_build.manifest.tsv
    printf "mode\\tstub\\n"                           >> pas_reference_build.manifest.tsv

    printf "PAS_REFERENCE_BUILD stub executed\\n" > pas_reference_build.log
    """
}
