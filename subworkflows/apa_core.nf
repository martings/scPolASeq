/*
 * APA core orchestration boundary.
 *
 * This subworkflow freezes the stage-3+ contract surface so new APA methods can
 * plug in as sidecar modules without changing main.nf or the legacy outputs
 * that downstream users already consume.
 */
include { BARCODE_FILTERING      } from './local/barcode_filtering/main'
include { GROUPED_RECONSTRUCTION } from './local/grouped_reconstruction/main'
include { COVERAGE_GENERATION    } from './local/coverage_generation/main'
include { APA_FEATURE_PIPELINE   } from './local/apa_feature_pipeline/main'
include { MODEL_PIPELINE         } from './local/model_pipeline/main'
include { PAS_REFERENCE_BUILD    } from '../modules/local/apa/pas_reference_build/main'
include { SIERRA_QUANT           } from '../modules/local/apa/sierra_quant/main'
include { PAS_SCORING            } from '../modules/local/apa/pas_scoring/main'

workflow APA_CORE {
    take:
    reference_bundle                   // tuple(ref_meta, star_index, gtf, fasta, chrom_sizes, terminal_exons, atlas, blacklist)
    bam_bundle                         // tuple(meta, bam, bai)
    cell_annotations                   // path: unified cell annotations TSV
    group_map                          // path: barcode-to-group map TSV
    known_polya                        // path: known PAS table or assets/NO_FILE sentinel
    apa_grouping_levels                // val: comma-separated grouping levels
    apa_min_coverage                   // val: minimum coverage threshold for APA calling
    apa_min_pdui_delta                 // val: minimum |delta PDUI| threshold
    apa_model_type                     // val: model type selector
    apa_model_group_level              // val: group level used by training and prediction

    main:
    def do_single_cell_projection = params.enable_single_cell_apa_projection.toString().toBoolean()
    def do_pas_reference_build = params.enable_pas_reference_build.toString().toBoolean()
    def do_sierra_quant = params.enable_sierra_quant.toString().toBoolean()
    def do_pas_scoring = params.enable_pas_scoring.toString().toBoolean()

    // Freeze the reference tuple layout here so downstream modules never index
    // into the bundle ad hoc.
    reference_bundle
        .map { ref_meta, star_index, gtf, fasta, chrom_sizes, terminal_exons, atlas, blacklist -> chrom_sizes }
        .first()
        .set { ch_chrom_sizes }

    reference_bundle
        .map { ref_meta, star_index, gtf, fasta, chrom_sizes, terminal_exons, atlas, blacklist -> terminal_exons }
        .first()
        .set { ch_terminal_exons }

    reference_bundle
        .map { ref_meta, star_index, gtf, fasta, chrom_sizes, terminal_exons, atlas, blacklist -> atlas }
        .first()
        .set { ch_site_catalog }

    BARCODE_FILTERING(
        bam_bundle,
        cell_annotations
    )

    GROUPED_RECONSTRUCTION(
        BARCODE_FILTERING.out.filtered_bam,
        group_map,
        apa_grouping_levels
    )

    COVERAGE_GENERATION(
        GROUPED_RECONSTRUCTION.out.grouped_bams,
        ch_chrom_sizes
    )

    APA_FEATURE_PIPELINE(
        ch_site_catalog,
        COVERAGE_GENERATION.out.bedgraphs,
        cell_annotations,
        known_polya,
        apa_min_coverage,
        apa_min_pdui_delta
    )

    MODEL_PIPELINE(
        APA_FEATURE_PIPELINE.out.feature_table,
        APA_FEATURE_PIPELINE.out.apa_events,
        apa_model_type,
        apa_model_group_level,
        do_single_cell_projection
    )

    def ch_pas_reference
    def ch_pas_reference_manifest
    if (do_pas_reference_build) {
        PAS_REFERENCE_BUILD(
            ch_site_catalog,
            ch_terminal_exons,
            known_polya
        )
        ch_pas_reference = PAS_REFERENCE_BUILD.out.pas_reference
        ch_pas_reference_manifest = PAS_REFERENCE_BUILD.out.manifest
    } else {
        // Disabled mode still emits a stable effective PAS reference by falling
        // back to the annotation-guided site catalog already used today.
        ch_pas_reference = ch_site_catalog
        ch_pas_reference_manifest = Channel.empty()
    }

    def ch_sierra_quant
    def ch_sierra_manifest
    if (do_sierra_quant) {
        GROUPED_RECONSTRUCTION.out.grouped_bams
            .flatMap { meta, group_level, bams, bais ->
                def bam_list = (bams instanceof List) ? bams : [bams]
                def bai_list = (bais instanceof List) ? bais : [bais]
                [bam_list, bai_list].transpose().collect { bam, bai ->
                    def group_id = bam.baseName.replaceAll(/\.grouped(\.bam)?$/, '')
                    tuple(meta, group_level, group_id, bam, bai)
                }
            }
            .combine(ch_pas_reference)
            .combine(cell_annotations)
            .map { meta, group_level, group_id, bam, bai, pas_reference, annotations ->
                tuple(meta, group_level, group_id, bam, bai, pas_reference, annotations)
            }
            .set { ch_sierra_input }

        SIERRA_QUANT(ch_sierra_input)
        ch_sierra_quant = SIERRA_QUANT.out.quantified_pas
        ch_sierra_manifest = SIERRA_QUANT.out.manifest
    } else {
        ch_sierra_quant = Channel.empty()
        ch_sierra_manifest = Channel.empty()
    }

    def ch_pas_scored
    def ch_pas_scoring_manifest
    if (do_pas_scoring) {
        PAS_SCORING(
            APA_FEATURE_PIPELINE.out.feature_table,
            APA_FEATURE_PIPELINE.out.apa_events,
            APA_FEATURE_PIPELINE.out.pdui_matrix,
            MODEL_PIPELINE.out.scored_events,
            ch_pas_reference
        )
        ch_pas_scored = PAS_SCORING.out.scored_pas
        ch_pas_scoring_manifest = PAS_SCORING.out.manifest
    } else {
        ch_pas_scored = Channel.empty()
        ch_pas_scoring_manifest = Channel.empty()
    }

    def ch_track_bundle = COVERAGE_GENERATION.out.bigwigs
        .flatMap { meta, group_level, group_id, fwd_bw, rev_bw -> [fwd_bw, rev_bw] }

    // Keep the reporting QC surface backward compatible; new sidecar manifests
    // are emitted separately so they can be adopted intentionally later.
    def ch_qc_bundle = BARCODE_FILTERING.out.filter_stats
        .mix(MODEL_PIPELINE.out.model_metrics)
        .mix(GROUPED_RECONSTRUCTION.out.grouping_manifest)
        .mix(ch_pas_reference_manifest)

    emit:
    filtered_bam                = BARCODE_FILTERING.out.filtered_bam
    filter_stats                = BARCODE_FILTERING.out.filter_stats
    grouped_bams                = GROUPED_RECONSTRUCTION.out.grouped_bams
    grouping_manifest           = GROUPED_RECONSTRUCTION.out.grouping_manifest
    bedgraphs                   = COVERAGE_GENERATION.out.bedgraphs
    bigwigs                     = COVERAGE_GENERATION.out.bigwigs
    coverage_stats              = COVERAGE_GENERATION.out.coverage_stats
    feature_table               = APA_FEATURE_PIPELINE.out.feature_table
    apa_events                  = APA_FEATURE_PIPELINE.out.apa_events
    pdui_matrix                 = APA_FEATURE_PIPELINE.out.pdui_matrix
    trained_model               = MODEL_PIPELINE.out.trained_model
    feature_importance          = MODEL_PIPELINE.out.feature_importance
    model_metrics               = MODEL_PIPELINE.out.model_metrics
    scored_events               = MODEL_PIPELINE.out.scored_events
    pas_reference               = ch_pas_reference
    pas_reference_manifest      = ch_pas_reference_manifest
    sierra_quant                = ch_sierra_quant
    sierra_manifest             = ch_sierra_manifest
    pas_scored_events           = ch_pas_scored
    pas_scoring_manifest        = ch_pas_scoring_manifest
    track_bundle                = ch_track_bundle
    qc_bundle                   = ch_qc_bundle
}
