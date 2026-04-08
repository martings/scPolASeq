/*
 * Scaffold-only PAS reference builder.
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
    printf "pas_reference_id\\tsite_id\\tgene_id\\tchrom\\tstart\\tend\\tstrand\\treference_source\\n" > pas_reference.tsv
    printf "pas_ref_stub_001\\tNA\\tNA\\tNA\\t0\\t0\\t+\\tcontract_scaffold\\n" >> pas_reference.tsv

    printf "field\\tvalue\\n" > pas_reference_build.manifest.tsv
    printf "site_catalog\\t%s\\n" "${site_catalog.name}" >> pas_reference_build.manifest.tsv
    printf "terminal_exons\\t%s\\n" "${terminal_exons.name}" >> pas_reference_build.manifest.tsv
    printf "known_polya\\t%s\\n" "${known_polya.name}" >> pas_reference_build.manifest.tsv
    printf "mode\\tscaffold_only\\n" >> pas_reference_build.manifest.tsv

    printf "PAS_REFERENCE_BUILD scaffold executed\\n" > pas_reference_build.log
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
