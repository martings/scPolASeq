include { PREPARE_REFERENCE_BUNDLE } from '../../../modules/local/reference/prepare_reference_bundle/main'
include { BUILD_TERMINAL_EXON_CATALOG } from '../../../modules/local/apa/build_terminal_exon_catalog/main'

workflow REFERENCE_PREPARE {
    take:
    genome_fasta
    gtf
    known_polya
    priming_blacklist

    main:
    PREPARE_REFERENCE_BUNDLE(genome_fasta, gtf, known_polya, priming_blacklist)
    BUILD_TERMINAL_EXON_CATALOG(PREPARE_REFERENCE_BUNDLE.out.gtf)

    PREPARE_REFERENCE_BUNDLE.out.reference_meta
        .combine(BUILD_TERMINAL_EXON_CATALOG.out.terminal_exons)
        .map { reference_row, terminal_exons ->
            def (meta, star_index, normalized_gtf, fasta, chrom_sizes, atlas, blacklist) = reference_row
            tuple(meta, star_index, normalized_gtf, fasta, chrom_sizes, terminal_exons, atlas, blacklist)
        }
        .set { ch_reference_bundle }

    emit:
    reference_bundle  = ch_reference_bundle
    reference_manifest = PREPARE_REFERENCE_BUNDLE.out.manifest
    terminal_exons    = BUILD_TERMINAL_EXON_CATALOG.out.terminal_exons
}
