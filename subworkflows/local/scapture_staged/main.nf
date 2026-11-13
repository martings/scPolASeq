/*
 * Opt-in SCAPTURE path with a cache boundary at each expensive stage.
 * Input: tuple(meta, filtered_bam, bai, fasta, fai, gtf, chrom_sizes,
 *              cell_annotations, polya_db, has_polya_db).
 * Optional val_catalog_only (default false) runs discovery + catalog + manifest
 * without running PASquant (#28).
 */
include { SCAPTURE_DISCOVERY } from '../scapture_discovery/main'
include { SCAPTURE_CATALOG } from '../../../modules/local/apa/scapture_catalog/main'
include { SCAPTURE_PASQUANT } from '../../../modules/local/apa/scapture_pasquant/main'
include { SCAPTURE_MANIFEST } from '../../../modules/local/apa/scapture_manifest/main'

workflow SCAPTURE_STAGED {
    take:
    ch_input
    val_catalog_only

    main:
    def catalog_only = val_catalog_only ? true : false

    ch_discovery = ch_input.map { meta, bam, bai, fasta, fai, gtf, chrom_sizes, annotations, polya_db, has_polya_db ->
        tuple(meta, bam, bai, fasta, fai, gtf, chrom_sizes, polya_db, has_polya_db)
    }
    SCAPTURE_DISCOVERY(ch_discovery)

    ch_catalog = SCAPTURE_DISCOVERY.out.evaluated.map { meta, exonic, intronic, extended ->
        tuple(meta, exonic, intronic, extended,
              file("${moduleDir}/../../../bin/run_scapture_filter.py", checkIfExists: true))
    }
    SCAPTURE_CATALOG(ch_catalog)

    def ch_quant_matrix
    def ch_quant_manifest
    def ch_quant_metrics
    if (!catalog_only) {
        // The reference and BAM are keyed by library ID. A missing or duplicate
        // library fails the join instead of silently quantifying the wrong sample.
        ch_bams = ch_input.map { meta, bam, bai, fasta, fai, gtf, chrom_sizes, annotations, polya_db, has_polya_db ->
            tuple(meta.library_id, meta, bam, bai, annotations)
        }
        ch_quant = SCAPTURE_DISCOVERY.out.reference
            .map { meta, reference -> tuple(meta.library_id, reference) }
            .join(ch_bams, failOnDuplicate: true, failOnMismatch: true)
            .map { library_id, reference, meta, bam, bai, annotations ->
                tuple(meta, bam, bai, reference, annotations,
                      file("${moduleDir}/../../../bin/FilterBamByTag", checkIfExists: true))
            }
        SCAPTURE_PASQUANT(ch_quant)
        ch_quant_matrix   = SCAPTURE_PASQUANT.out.matrix
        ch_quant_manifest = SCAPTURE_PASQUANT.out.manifest
        ch_quant_metrics  = SCAPTURE_PASQUANT.out.metrics
    } else {
        ch_quant_matrix   = Channel.empty()
        ch_quant_manifest = SCAPTURE_CATALOG.out.metrics.map { meta, metrics ->
            tuple(meta, file("${moduleDir}/../../../assets/NO_FILE", checkIfExists: true))
        }
        ch_quant_metrics  = Channel.empty()
    }

    ch_manifest = ch_input
        .map { meta, bam, bai, fasta, fai, gtf, chrom_sizes, annotations, polya_db, has_polya_db ->
            tuple(meta.library_id, meta, bam, gtf, fasta, fai)
        }
        .join(SCAPTURE_DISCOVERY.out.prepared_gtf.map { meta, gtf -> tuple(meta.library_id, gtf) }, failOnDuplicate: true, failOnMismatch: true)
        .join(SCAPTURE_DISCOVERY.out.gene_name_map.map { meta, mapping -> tuple(meta.library_id, mapping) }, failOnDuplicate: true, failOnMismatch: true)
        .join(SCAPTURE_DISCOVERY.out.duplicate_gene_names.map { meta, report -> tuple(meta.library_id, report) }, failOnDuplicate: true, failOnMismatch: true)
        .join(SCAPTURE_DISCOVERY.out.pascall_metrics.map { meta, metrics -> tuple(meta.library_id, metrics) }, failOnDuplicate: true, failOnMismatch: true)
        .join(SCAPTURE_DISCOVERY.out.reference_metrics.map { meta, metrics -> tuple(meta.library_id, metrics) }, failOnDuplicate: true, failOnMismatch: true)
        .join(SCAPTURE_CATALOG.out.metrics.map { meta, metrics -> tuple(meta.library_id, metrics) }, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_quant_manifest.map { meta, manifest -> tuple(meta.library_id, manifest) }, failOnDuplicate: true, failOnMismatch: true)
        .map { row ->
            tuple(row[1], row[2], row[3], row[4], row[5], row[6], row[7],
                  row[8], row[9], row[10], row[11], row[12],
                  file("${moduleDir}/../../../bin/scapture_manifest.py", checkIfExists: true))
        }
    SCAPTURE_MANIFEST(ch_manifest)

    emit:
    site_catalog = SCAPTURE_CATALOG.out.catalog
    quant = ch_quant_matrix
    manifest = SCAPTURE_MANIFEST.out.manifest
    quant_manifest = catalog_only ? Channel.empty() : ch_quant_manifest
    catalog_metrics = SCAPTURE_CATALOG.out.metrics
    pasquant_metrics = ch_quant_metrics
    reference = SCAPTURE_DISCOVERY.out.reference
    reference_metrics = SCAPTURE_DISCOVERY.out.reference_metrics
    merge_status = SCAPTURE_DISCOVERY.out.merge_status
}
