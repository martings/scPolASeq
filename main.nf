nextflow.enable.dsl = 2

// ── Existing pre-processing subworkflows ──────────────────────────────────────
include { INPUT_HARMONIZATION    } from './subworkflows/local/input_harmonization/main'
include { REFERENCE_PREPARE      } from './subworkflows/local/reference_prepare/main'
include { ALIGN_OR_IMPORT_SC     } from './subworkflows/local/alignment_or_import/main'
include { CELL_LABELING_SC       } from './subworkflows/local/cell_labeling/main'

// ── APA orchestration boundary ────────────────────────────────────────────────
include { APA_CORE               } from './subworkflows/apa_core'

// ── Reporting ─────────────────────────────────────────────────────────────────
include { REPORTING_SC           } from './subworkflows/local/reporting/main'

workflow SCPOLASEQ {

    // ── Parameter validation ──────────────────────────────────────────────────
    if (!params.input)        error "Missing required parameter: --input"
    if (!params.genome_fasta) error "Missing required parameter: --genome_fasta"
    if (!params.gtf)          error "Missing required parameter: --gtf"

    // ── Stage 0 — Input harmonization + Reference preparation ─────────────────
    // Returns true only if value is a real path (not null, not keyword like "null"/"skip"/"download")
    def _NULL_STRINGS = ['', 'null', 'none', 'skip', 'false', 'no', '0', 'no_file'] as Set
    def param_is_set = { value ->
        def text = value?.toString()?.trim()
        return text && !(text.toLowerCase() in _NULL_STRINGS)
    }
    def is_local_polya_source = { value ->
        if (!param_is_set(value)) return false
        def lower = value.toString().trim().toLowerCase()
        def non_local = ['download', 'auto', 'true', 'yes', '1', 'http://', 'https://', 'ftp://']
        return !non_local.any { lower == it || lower.startsWith(it) }
    }
    def polya_db_source = params.polya_db ?: 'skip'
    def polyasite_source = params.polyasite ?: 'skip'
    def polya_db_file = is_local_polya_source(polya_db_source) ? file(polya_db_source, checkIfExists: true) : file("${projectDir}/assets/NO_FILE_POLYA_DB")
    def polyasite_file = is_local_polya_source(polyasite_source) ? file(polyasite_source, checkIfExists: true) : file("${projectDir}/assets/NO_FILE_POLYASITE")

    INPUT_HARMONIZATION(
        file(params.input, checkIfExists: true),
        param_is_set(params.cell_metadata)       ? file(params.cell_metadata,       checkIfExists: true) : file("${projectDir}/assets/NO_FILE_CELL_METADATA"),
        param_is_set(params.cluster_assignments) ? file(params.cluster_assignments, checkIfExists: true) : file("${projectDir}/assets/NO_FILE_CLUSTER_ASSIGNMENTS"),
        param_is_set(params.cell_type_labels)    ? file(params.cell_type_labels,    checkIfExists: true) : file("${projectDir}/assets/NO_FILE_CELL_TYPE_LABELS")
    )

    REFERENCE_PREPARE(
        file(params.genome_fasta, checkIfExists: true),
        file(params.gtf,          checkIfExists: true),
        param_is_set(params.known_polya)       ? file(params.known_polya,       checkIfExists: true) : file("${projectDir}/assets/NO_FILE"),
        param_is_set(params.priming_blacklist) ? file(params.priming_blacklist, checkIfExists: true) : file("${projectDir}/assets/NO_FILE_PRIMING_BLACKLIST"),
        polya_db_source,
        polya_db_file,
        polyasite_source,
        polyasite_file
    )

    // ── Stage 1 — Initial alignment (STARsolo, CB/UB tag preservation) ────────
    ALIGN_OR_IMPORT_SC(
        INPUT_HARMONIZATION.out.samplesheet,
        REFERENCE_PREPARE.out.reference_bundle,
        params.run_mode,
        params.aligner,
        params.protocol_mode
    )

    // ── Stage 2 — Expression quantification + cell labeling ───────────────────
    def ch_celltypist_model = params.celltypist_model
        ? Channel.value(file(params.celltypist_model, checkIfExists: true))
        : Channel.value(file("${projectDir}/assets/NO_FILE"))
    def ch_cell_type_ref_h5ad = params.cell_type_ref_h5ad
        ? Channel.value(file(params.cell_type_ref_h5ad, checkIfExists: true))
        : Channel.value(file("${projectDir}/assets/NO_FILE"))

    CELL_LABELING_SC(
        ALIGN_OR_IMPORT_SC.out.matrix_bundle,
        INPUT_HARMONIZATION.out.cell_metadata,
        params.enable_internal_clustering,
        params.grouping_levels,
        ch_celltypist_model,
        ch_cell_type_ref_h5ad,
        params.cell_type_ref_label_col ?: 'cell_type'
    )

    // ── Stages 3–9 — APA core orchestration ───────────────────────────────────
    // All APA-stage channel wiring is isolated behind APA_CORE so future
    // modules can be added without widening main.nf.
    APA_CORE(
        REFERENCE_PREPARE.out.reference_bundle,
        ALIGN_OR_IMPORT_SC.out.bam_bundle,
        CELL_LABELING_SC.out.cell_annotations,
        CELL_LABELING_SC.out.group_map,
        REFERENCE_PREPARE.out.known_polya,
        params.apa_grouping_levels ?: 'cluster,cell_type',
        params.apa_min_coverage ?: 5,
        params.apa_min_pdui_delta ?: 0.2,
        params.apa_model_type ?: 'random_forest',
        params.apa_group_level ?: 'cluster'
    )

    // ── Reporting ─────────────────────────────────────────────────────────────
    REPORTING_SC(
        APA_CORE.out.feature_table,   // site_catalog proxy
        APA_CORE.out.pdui_matrix,     // apa_usage
        APA_CORE.out.scored_events,   // apa_stats
        APA_CORE.out.track_bundle,
        APA_CORE.out.qc_bundle
    )
}

workflow {
    SCPOLASEQ()
}
