nextflow.enable.dsl = 2

include { INPUT_HARMONIZATION } from './subworkflows/local/input_harmonization/main'
include { REFERENCE_PREPARE    } from './subworkflows/local/reference_prepare/main'
include { ALIGN_OR_IMPORT_SC   } from './subworkflows/local/alignment_or_import/main'
include { CELL_LABELING_SC     } from './subworkflows/local/cell_labeling/main'
include { FEEDBACK_RECONSTRUCT } from './subworkflows/local/feedback_reconstruct/main'
include { APA_ANALYSIS_SC      } from './subworkflows/local/apa_analysis/main'
include { REPORTING_SC         } from './subworkflows/local/reporting/main'

workflow SCPOLASEQ {
    if (!params.input) {
        error "Missing required parameter: --input"
    }
    if (!params.genome_fasta) {
        error "Missing required parameter: --genome_fasta"
    }
    if (!params.gtf) {
        error "Missing required parameter: --gtf"
    }

    INPUT_HARMONIZATION(
        file(params.input, checkIfExists: true),
        params.cell_metadata ? file(params.cell_metadata, checkIfExists: true) : file("${projectDir}/assets/NO_FILE"),
        params.cluster_assignments ? file(params.cluster_assignments, checkIfExists: true) : file("${projectDir}/assets/NO_FILE"),
        params.cell_type_labels ? file(params.cell_type_labels, checkIfExists: true) : file("${projectDir}/assets/NO_FILE")
    )

    REFERENCE_PREPARE(
        file(params.genome_fasta, checkIfExists: true),
        file(params.gtf, checkIfExists: true),
        params.known_polya ? file(params.known_polya, checkIfExists: true) : file("${projectDir}/assets/NO_FILE"),
        params.priming_blacklist ? file(params.priming_blacklist, checkIfExists: true) : file("${projectDir}/assets/NO_FILE")
    )

    ALIGN_OR_IMPORT_SC(
        INPUT_HARMONIZATION.out.samplesheet,
        REFERENCE_PREPARE.out.reference_bundle,
        params.run_mode,
        params.aligner,
        params.protocol_mode
    )

    CELL_LABELING_SC(
        ALIGN_OR_IMPORT_SC.out.matrix_bundle,
        INPUT_HARMONIZATION.out.cell_metadata,
        params.enable_internal_clustering,
        params.grouping_levels
    )

    FEEDBACK_RECONSTRUCT(
        ALIGN_OR_IMPORT_SC.out.bam_bundle,
        ALIGN_OR_IMPORT_SC.out.barcode_registry,
        CELL_LABELING_SC.out.cell_annotations,
        CELL_LABELING_SC.out.group_map,
        REFERENCE_PREPARE.out.reference_bundle,
        params.emit_group_bams,
        params.emit_bigwigs
    )

    APA_ANALYSIS_SC(
        FEEDBACK_RECONSTRUCT.out.coverage_bundle,
        FEEDBACK_RECONSTRUCT.out.group_map,
        CELL_LABELING_SC.out.cell_annotations,
        REFERENCE_PREPARE.out.reference_bundle
    )

    REPORTING_SC(
        APA_ANALYSIS_SC.out.site_catalog,
        APA_ANALYSIS_SC.out.apa_usage,
        APA_ANALYSIS_SC.out.apa_stats,
        FEEDBACK_RECONSTRUCT.out.track_bundle,
        FEEDBACK_RECONSTRUCT.out.qc_bundle
    )
}

workflow {
    SCPOLASEQ()
}
