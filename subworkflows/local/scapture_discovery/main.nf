/*
 * Incremental SCAPTURE discovery. This subworkflow stops after DeepPASS and
 * intentionally does not replace SCAPTURE_FILTER until PASmerge/PASquant and
 * real-output equivalence have been validated (issues #27-29).
 *
 * Input: tuple(meta, bam, bai, fasta, fai, gtf, chrom_sizes, polya_db,
 *              has_polya_db). library_id must be unique in this channel.
 */
include { SCAPTURE_PREPARE } from '../../../modules/local/apa/scapture_prepare/main'
include { SCAPTURE_ANNOTATION } from '../../../modules/local/apa/scapture_annotation/main'
include { SCAPTURE_PASCALL } from '../../../modules/local/apa/scapture_pascall/main'

workflow SCAPTURE_DISCOVERY {
    take:
    ch_input

    main:
    ch_by_library = ch_input.map { meta, bam, bai, fasta, fai, gtf, chrom_sizes, polya_db, has_polya_db ->
        tuple(meta.library_id, meta, bam, bai, fasta, fai, chrom_sizes, polya_db, has_polya_db)
    }
    SCAPTURE_PREPARE(ch_input.map { meta, bam, bai, fasta, fai, gtf, chrom_sizes, polya_db, has_polya_db ->
        tuple(meta, gtf, file("${moduleDir}/../../../bin/prepare_scapture_gtf.py", checkIfExists: true))
    })

    ch_prepared = SCAPTURE_PREPARE.out.gtf
        .map { meta, gtf -> tuple(meta.library_id, gtf) }
        .join(ch_by_library, failOnDuplicate: true, failOnMismatch: true)
        .map { id, scapture_gtf, meta, bam, bai, fasta, fai, chrom_sizes, polya_db, has_polya_db ->
            tuple(meta, scapture_gtf, fasta, fai, chrom_sizes)
        }
    SCAPTURE_ANNOTATION(ch_prepared)

    ch_annotated = SCAPTURE_ANNOTATION.out.annotation
        .map { meta, annotation_files -> tuple(meta.library_id, annotation_files) }
        .join(ch_by_library, failOnDuplicate: true, failOnMismatch: true)
        .map { id, annotation_files, meta, bam, bai, fasta, fai, chrom_sizes, polya_db, has_polya_db ->
            tuple(meta, bam, bai, fasta, fai, annotation_files, polya_db, has_polya_db,
                  file("${moduleDir}/../../../bin/DeepPASS_Predict.py", checkIfExists: true))
        }
    SCAPTURE_PASCALL(ch_annotated)

    emit:
    evaluated = SCAPTURE_PASCALL.out.evaluated
    pascall_metrics = SCAPTURE_PASCALL.out.metrics
    gene_name_map = SCAPTURE_PREPARE.out.gene_name_map
    duplicate_gene_names = SCAPTURE_PREPARE.out.duplicate_gene_names
    annotation_log = SCAPTURE_ANNOTATION.out.log
    pascall_log = SCAPTURE_PASCALL.out.log
    prepare_versions = SCAPTURE_PREPARE.out.versions
    annotation_versions = SCAPTURE_ANNOTATION.out.versions
    pascall_versions = SCAPTURE_PASCALL.out.versions
}
