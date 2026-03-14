// Stages 6 + 7 subworkflow — Feature extraction and statistical APA calling.
include { APA_FEATURE_EXTRACTION } from '../../../modules/local/apa/apa_feature_extraction/main'
include { APA_CALLING            } from '../../../modules/local/apa/apa_calling/main'

workflow APA_FEATURE_PIPELINE {
    take:
    site_catalog      // path
    bedgraphs         // channel: tuple(meta, group_level, group_id, fwd_bg, rev_bg)
    cell_annotations  // path — unified cell annotations TSV
    known_polya       // path (may be NO_FILE)
    min_coverage      // val
    min_pdui_delta    // val

    main:
    // Collect all bedGraph files (fwd + rev) into one flat list for the extraction script
    bedgraphs
        .flatMap { meta, group_level, group_id, fwd_bg, rev_bg -> [fwd_bg, rev_bg] }
        .collect()
        .set { ch_all_bedgraphs }

    // cell_annotations is a plain path — use directly
    APA_FEATURE_EXTRACTION(
        site_catalog,
        ch_all_bedgraphs,
        cell_annotations,
        known_polya
    )

    APA_CALLING(
        APA_FEATURE_EXTRACTION.out.feature_table,
        site_catalog,
        cell_annotations,
        min_coverage,
        min_pdui_delta
    )

    emit:
    feature_table = APA_FEATURE_EXTRACTION.out.feature_table  // path: apa_features.tsv
    apa_events    = APA_CALLING.out.apa_events                 // path: apa_events.tsv
    pdui_matrix   = APA_CALLING.out.pdui_matrix                // path: pdui_usage_matrix.tsv
}
