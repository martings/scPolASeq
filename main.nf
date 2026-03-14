nextflow.enable.dsl = 2

// ── Existing pre-processing subworkflows ──────────────────────────────────────
include { INPUT_HARMONIZATION    } from './subworkflows/local/input_harmonization/main'
include { REFERENCE_PREPARE      } from './subworkflows/local/reference_prepare/main'
include { ALIGN_OR_IMPORT_SC     } from './subworkflows/local/alignment_or_import/main'
include { CELL_LABELING_SC       } from './subworkflows/local/cell_labeling/main'

// ── New 9-stage APA subworkflows ──────────────────────────────────────────────
include { BARCODE_FILTERING      } from './subworkflows/local/barcode_filtering/main'
include { GROUPED_RECONSTRUCTION } from './subworkflows/local/grouped_reconstruction/main'
include { COVERAGE_GENERATION    } from './subworkflows/local/coverage_generation/main'
include { APA_FEATURE_PIPELINE   } from './subworkflows/local/apa_feature_pipeline/main'
include { MODEL_PIPELINE         } from './subworkflows/local/model_pipeline/main'

// ── Reporting ─────────────────────────────────────────────────────────────────
include { REPORTING_SC           } from './subworkflows/local/reporting/main'

workflow SCPOLASEQ {

    // ── Parameter validation ──────────────────────────────────────────────────
    if (!params.input)        error "Missing required parameter: --input"
    if (!params.genome_fasta) error "Missing required parameter: --genome_fasta"
    if (!params.gtf)          error "Missing required parameter: --gtf"

    def nofile = file("${projectDir}/assets/NO_FILE")

    // ── Stage 0 — Input harmonization + Reference preparation ─────────────────
    INPUT_HARMONIZATION(
        file(params.input, checkIfExists: true),
        params.cell_metadata        ? file(params.cell_metadata,        checkIfExists: true) : nofile,
        params.cluster_assignments  ? file(params.cluster_assignments,  checkIfExists: true) : nofile,
        params.cell_type_labels     ? file(params.cell_type_labels,     checkIfExists: true) : nofile
    )

    REFERENCE_PREPARE(
        file(params.genome_fasta, checkIfExists: true),
        file(params.gtf,          checkIfExists: true),
        params.known_polya       ? file(params.known_polya,       checkIfExists: true) : nofile,
        params.priming_blacklist ? file(params.priming_blacklist, checkIfExists: true) : nofile
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
    CELL_LABELING_SC(
        ALIGN_OR_IMPORT_SC.out.matrix_bundle,
        INPUT_HARMONIZATION.out.cell_metadata,
        params.enable_internal_clustering,
        params.grouping_levels
    )

    // ── Stage 3 — Barcode-aware BAM filtering ─────────────────────────────────
    // Reads are subset to valid cell barcodes; each read is tagged with its
    // cluster / cell_type derived from CELL_LABELING_SC.
    BARCODE_FILTERING(
        ALIGN_OR_IMPORT_SC.out.bam_bundle,
        CELL_LABELING_SC.out.cell_annotations          // plain path (unified TSV)
    )

    // ── Stage 4 — Grouped BAM reconstruction (one BAM per group_level/group_id)
    GROUPED_RECONSTRUCTION(
        BARCODE_FILTERING.out.filtered_bam,
        CELL_LABELING_SC.out.group_map,                // plain path
        params.apa_grouping_levels ?: 'cluster,cell_type'
    )

    // ── Stage 5 — Strand-aware coverage tracks (bedGraph + bigWig) ───────────
    // Extract chrom_sizes from the reference bundle (field index 4)
    REFERENCE_PREPARE.out.reference_bundle
        .map { ref_meta, star_index, gtf, fasta, chrom_sizes,
               terminal_exons, atlas, blacklist -> chrom_sizes }
        .first()
        .set { ch_chrom_sizes }

    COVERAGE_GENERATION(
        GROUPED_RECONSTRUCTION.out.grouped_bams,
        ch_chrom_sizes
    )

    // ── Stages 6 + 7 — APA feature extraction + statistical APA calling ───────
    // site_catalog / atlas from reference bundle (field index 6)
    REFERENCE_PREPARE.out.reference_bundle
        .map { ref_meta, star_index, gtf, fasta, chrom_sizes,
               terminal_exons, atlas, blacklist -> atlas }
        .first()
        .set { ch_site_catalog }

    def ch_known_polya = params.known_polya
        ? Channel.value(file(params.known_polya, checkIfExists: true))
        : Channel.value(nofile)

    APA_FEATURE_PIPELINE(
        ch_site_catalog,
        COVERAGE_GENERATION.out.bedgraphs,
        CELL_LABELING_SC.out.cell_annotations,
        ch_known_polya,
        params.apa_min_coverage   ?: 5,
        params.apa_min_pdui_delta ?: 0.2
    )

    // ── Stages 8 + 9 — ML model training + APA event scoring ─────────────────
    MODEL_PIPELINE(
        APA_FEATURE_PIPELINE.out.feature_table,
        APA_FEATURE_PIPELINE.out.apa_events,
        params.apa_model_type  ?: 'random_forest',
        params.apa_group_level ?: 'cluster',
        params.enable_single_cell_apa_projection ?: false
    )

    // ── Reporting ─────────────────────────────────────────────────────────────
    // Map bigwig tuples to flat paths for track_bundle; mix QC files for qc_bundle.
    def ch_track_bundle = COVERAGE_GENERATION.out.bigwigs
        .flatMap { meta, group_level, group_id, fwd_bw, rev_bw -> [fwd_bw, rev_bw] }

    def ch_qc_bundle = BARCODE_FILTERING.out.filter_stats
        .mix(MODEL_PIPELINE.out.model_metrics)
        .mix(GROUPED_RECONSTRUCTION.out.grouping_manifest)

    REPORTING_SC(
        APA_FEATURE_PIPELINE.out.feature_table,   // site_catalog proxy
        APA_FEATURE_PIPELINE.out.pdui_matrix,     // apa_usage
        MODEL_PIPELINE.out.scored_events,          // apa_stats
        ch_track_bundle,
        ch_qc_bundle
    )
}

workflow {
    SCPOLASEQ()
}
