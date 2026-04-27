include { STARSOLO_ALIGN_SC } from '../../../modules/local/alignment/starsolo_align_sc/main'
include { IMPORT_VALIDATE_SC } from '../../../modules/local/bam/import_validate_sc/main'

workflow ALIGN_OR_IMPORT_SC {
    take:
    samplesheet
    reference_bundle
    run_mode
    aligner
    protocol_mode

    main:
    if (aligner != 'starsolo') {
        log.warn("Aligner '${aligner}' is not yet fully implemented in the scaffold. Falling back to STARSOLO_ALIGN_SC for FASTQ inputs.")
    }

    def nofile = file("${projectDir}/assets/NO_FILE")

    samplesheet
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                sample_id   : row.sample_id,
                library_id  : row.library_id,
                protocol    : row.protocol ?: protocol_mode,
                chemistry   : row.chemistry,
                condition   : row.condition,
                replicate_id: row.replicate_id ?: row.sample_id,
            ]
            def reads = []
            if (row.fastq_r1) {
                reads << file(row.fastq_r1, checkIfExists: true)
            }
            if (row.fastq_r2) {
                reads << file(row.fastq_r2, checkIfExists: true)
            }
            def bam = row.bam ? file(row.bam, checkIfExists: true) : nofile
            def matrix = row.matrix_path ? file(row.matrix_path, checkIfExists: true) : nofile
            def whitelist = row.barcode_whitelist ? file(row.barcode_whitelist, checkIfExists: true) : nofile
            tuple(meta, reads, bam, matrix, whitelist)
        }
        .set { ch_rows }

    ch_rows
        .filter { meta, reads, bam, matrix, whitelist ->
            run_mode == 'full' && reads.size() >= 2
        }
        .combine(reference_bundle)
        .map { meta, reads, bam, matrix, whitelist, reference_meta, star_index, gtf, fasta, fasta_fai, chrom_sizes, terminal_exons, atlas, blacklist ->
            tuple(meta, reads, whitelist, star_index, gtf, meta.protocol)
        }
        .set { ch_fastq }

    ch_rows
        .filter { meta, reads, bam, matrix, whitelist ->
            bam.name != 'NO_FILE' && (run_mode in ['from_bam', 'feedback'] || reads.size() == 0)
        }
        .map { meta, reads, bam, matrix, whitelist ->
            tuple(meta, bam, matrix, whitelist)
        }
        .set { ch_bam }

    STARSOLO_ALIGN_SC(ch_fastq)
    IMPORT_VALIDATE_SC(ch_bam)

    STARSOLO_ALIGN_SC.out.bam_bundle.mix(IMPORT_VALIDATE_SC.out.bam_bundle).set { ch_bam_bundle }
    STARSOLO_ALIGN_SC.out.matrix_bundle.mix(IMPORT_VALIDATE_SC.out.matrix_bundle).set { ch_matrix_bundle }
    STARSOLO_ALIGN_SC.out.barcode_registry.mix(IMPORT_VALIDATE_SC.out.barcode_registry).set { ch_barcode_registry }
    STARSOLO_ALIGN_SC.out.alignment_manifest.mix(IMPORT_VALIDATE_SC.out.alignment_manifest).set { ch_alignment_manifest }

    emit:
    bam_bundle         = ch_bam_bundle
    matrix_bundle      = ch_matrix_bundle
    barcode_registry   = ch_barcode_registry
    alignment_manifest = ch_alignment_manifest
}
