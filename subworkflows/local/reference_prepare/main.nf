include { PREPARE_POLYA_REFERENCES } from '../../../modules/local/reference/prepare_polya_references/main'
include { PREPARE_REFERENCE_BUNDLE } from '../../../modules/local/reference/prepare_reference_bundle/main'
include { BUILD_TERMINAL_EXON_CATALOG } from '../../../modules/local/apa/build_terminal_exon_catalog/main'

workflow REFERENCE_PREPARE {
    take:
    genome_fasta
    gtf
    known_polya
    priming_blacklist
    polya_db
    polya_db_file
    polyasite
    polyasite_file

    main:
    def polya_skip_values = ['skip', 'none', 'null', 'false', 'no', '0', '']
    def polya_db_source = polya_db?.toString()?.trim() ?: 'skip'
    def polyasite_source = polyasite?.toString()?.trim() ?: 'skip'
    def wants_polya_db = !(polya_db_source.toLowerCase() in polya_skip_values)
    def wants_polyasite = !(polyasite_source.toLowerCase() in polya_skip_values)

    def ch_effective_known_polya
    def ch_polya_reference_manifest
    if (wants_polya_db || wants_polyasite) {
        PREPARE_POLYA_REFERENCES(
            known_polya,
            polya_db,
            polya_db_file,
            polyasite,
            polyasite_file
        )
        ch_effective_known_polya = PREPARE_POLYA_REFERENCES.out.known_polya
        ch_polya_reference_manifest = PREPARE_POLYA_REFERENCES.out.manifest
    } else {
        ch_effective_known_polya = known_polya
        ch_polya_reference_manifest = Channel.empty()
    }

    PREPARE_REFERENCE_BUNDLE(genome_fasta, gtf, ch_effective_known_polya, priming_blacklist)
    BUILD_TERMINAL_EXON_CATALOG(PREPARE_REFERENCE_BUNDLE.out.gtf)

    PREPARE_REFERENCE_BUNDLE.out.reference_meta
        .combine(BUILD_TERMINAL_EXON_CATALOG.out.terminal_exons)
        .map { meta, star_index, normalized_gtf, fasta, chrom_sizes, atlas, blacklist, terminal_exons ->
            tuple(meta, star_index, normalized_gtf, fasta, chrom_sizes, terminal_exons, atlas, blacklist)
        }
        .set { ch_reference_bundle }

    PREPARE_REFERENCE_BUNDLE.out.reference_meta
        .map { meta, star_index, normalized_gtf, fasta, chrom_sizes, atlas, blacklist -> atlas }
        .set { ch_known_polya_reference }

    emit:
    reference_bundle  = ch_reference_bundle
    reference_manifest = PREPARE_REFERENCE_BUNDLE.out.manifest
    polya_reference_manifest = ch_polya_reference_manifest
    known_polya      = ch_known_polya_reference
    terminal_exons    = BUILD_TERMINAL_EXON_CATALOG.out.terminal_exons
}
