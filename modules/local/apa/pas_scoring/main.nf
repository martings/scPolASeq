/*
 * Scaffold-only PAS scoring module.
 *
 * Contract:
 *   input:
 *     path feature_table
 *     path apa_events
 *     path pdui_matrix
 *     path scored_events
 *     path pas_reference
 *   output:
 *     path "pas_scored_events.tsv"
 *     path "pas_scoring.manifest.tsv"
 *     path "pas_scoring.log"
 */
process PAS_SCORING {
    tag "pas_scoring"
    label 'process_single'

    publishDir "${params.outdir}/pas_scoring", mode: params.publish_dir_mode

    input:
    path feature_table
    path apa_events
    path pdui_matrix
    path scored_events
    path pas_reference

    output:
    path "pas_scored_events.tsv",    emit: scored_pas
    path "pas_scoring.manifest.tsv", emit: manifest
    path "pas_scoring.log",          emit: log

    script:
    """
    printf "event_id\\tscore_namespace\\tscore_source\\treference_input\\n" > pas_scored_events.tsv
    printf "stub_event_001\\tpas_scoring\\tscaffold_only\\t${pas_reference.name}\\n" >> pas_scored_events.tsv

    printf "field\\tvalue\\n" > pas_scoring.manifest.tsv
    printf "feature_table\\t%s\\n" "${feature_table.name}" >> pas_scoring.manifest.tsv
    printf "apa_events\\t%s\\n" "${apa_events.name}" >> pas_scoring.manifest.tsv
    printf "pdui_matrix\\t%s\\n" "${pdui_matrix.name}" >> pas_scoring.manifest.tsv
    printf "scored_events\\t%s\\n" "${scored_events.name}" >> pas_scoring.manifest.tsv
    printf "pas_reference\\t%s\\n" "${pas_reference.name}" >> pas_scoring.manifest.tsv

    printf "PAS_SCORING scaffold executed\\n" > pas_scoring.log
    """

    stub:
    """
    printf "event_id\\tscore_namespace\\tscore_source\\treference_input\\n" > pas_scored_events.tsv
    printf "stub_event_001\\tpas_scoring\\tstub\\tstub\\n" >> pas_scored_events.tsv
    printf "field\\tvalue\\nmode\\tstub\\n" > pas_scoring.manifest.tsv
    printf "PAS_SCORING stub executed\\n" > pas_scoring.log
    """
}
