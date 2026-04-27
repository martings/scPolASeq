include { CANDIDATE_SITES_ANNOTATION_GUIDED } from '../../../modules/local/apa/candidate_sites_annotation_guided/main'
include { INTERNAL_PRIMING_FILTER           } from '../../../modules/local/apa/internal_priming_filter/main'
include { APA_QUANTIFY_SITE_USAGE           } from '../../../modules/local/apa/quantify_site_usage/main'
include { APA_TEST_DIFFERENTIAL_USAGE       } from '../../../modules/local/apa/test_differential_usage/main'

workflow APA_ANALYSIS_SC {
    take:
    coverage_bundle
    group_map
    cell_annotations
    reference_bundle

    main:
    reference_bundle
        .map { reference_meta, star_index, gtf, fasta, fasta_fai, chrom_sizes, terminal_exons, atlas, blacklist ->
            tuple(terminal_exons, atlas)
        }
        .set { ch_candidate_input }

    CANDIDATE_SITES_ANNOTATION_GUIDED(ch_candidate_input)

    CANDIDATE_SITES_ANNOTATION_GUIDED.out.site_catalog
        .combine(reference_bundle)
        .map { site_catalog, reference_meta, star_index, gtf, fasta, fasta_fai, chrom_sizes, terminal_exons, atlas, blacklist ->
            tuple(site_catalog, fasta, blacklist)
        }
        .set { ch_priming_input }

    INTERNAL_PRIMING_FILTER(ch_priming_input)
    APA_QUANTIFY_SITE_USAGE(coverage_bundle.collect(), INTERNAL_PRIMING_FILTER.out.site_catalog, group_map, cell_annotations)
    APA_TEST_DIFFERENTIAL_USAGE(APA_QUANTIFY_SITE_USAGE.out.apa_usage, cell_annotations)

    emit:
    site_catalog = INTERNAL_PRIMING_FILTER.out.site_catalog
    apa_usage    = APA_QUANTIFY_SITE_USAGE.out.apa_usage
    apa_stats    = APA_TEST_DIFFERENTIAL_USAGE.out.apa_stats
    qc_bundle    = APA_QUANTIFY_SITE_USAGE.out.gene_summary.mix(INTERNAL_PRIMING_FILTER.out.metrics).mix(APA_TEST_DIFFERENTIAL_USAGE.out.model_summary)
}
