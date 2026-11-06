/*
 * Opt-in SCAPTURE path with a cache boundary at each expensive stage.
 * Input: tuple(meta, filtered_bam, bai, fasta, fai, gtf, chrom_sizes,
 *              cell_annotations, polya_db, has_polya_db).
 */
include { SCAPTURE_DISCOVERY } from '../scapture_discovery/main'
include { SCAPTURE_CATALOG } from '../../../modules/local/apa/scapture_catalog/main'
include { SCAPTURE_PASQUANT } from '../../../modules/local/apa/scapture_pasquant/main'
include { SCAPTURE_MANIFEST } from '../../../modules/local/apa/scapture_manifest/main'

workflow SCAPTURE_STAGED {
    take:
    ch_input

    main:
    ch_discovery = ch_input.map { meta, bam, bai, fasta, fai, gtf, chrom_sizes, annotations, polya_db, has_polya_db ->
        tuple(meta, bam, bai, fasta, fai, gtf, chrom_sizes, polya_db, has_polya_db)
    }
    SCAPTURE_DISCOVERY(ch_discovery)

    ch_catalog = SCAPTURE_DISCOVERY.out.evaluated.map { meta, exonic, intronic, extended ->
        tuple(meta, exonic, intronic, extended,
              file("${moduleDir}/../../../bin/run_scapture_filter.py", checkIfExists: true))
    }
    SCAPTURE_CATALOG(ch_catalog)

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
        .join(SCAPTURE_PASQUANT.out.manifest.map { meta, manifest -> tuple(meta.library_id, manifest) }, failOnDuplicate: true, failOnMismatch: true)
        .map { row ->
            tuple(row[1], row[2], row[3], row[4], row[5], row[6], row[7],
                  row[8], row[9], row[10], row[11], row[12],
                  file("${moduleDir}/../../../bin/scapture_manifest.py", checkIfExists: true))
        }
    SCAPTURE_MANIFEST(ch_manifest)

    emit:
    site_catalog = SCAPTURE_CATALOG.out.catalog
    quant = SCAPTURE_PASQUANT.out.matrix
    manifest = SCAPTURE_MANIFEST.out.manifest
    quant_manifest = SCAPTURE_PASQUANT.out.manifest
    catalog_metrics = SCAPTURE_CATALOG.out.metrics
    pasquant_metrics = SCAPTURE_PASQUANT.out.metrics
    reference = SCAPTURE_DISCOVERY.out.reference
    reference_metrics = SCAPTURE_DISCOVERY.out.reference_metrics
    merge_status = SCAPTURE_DISCOVERY.out.merge_status
}
