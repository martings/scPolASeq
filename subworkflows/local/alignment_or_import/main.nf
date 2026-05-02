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
    def normalize = { value -> value ? value.toString().trim() : '' }
    def distinctOrThrow = { sampleId, entries, key, fallback = '' ->
        def values = entries.collect { entry -> normalize(entry[key]) }.findAll { it }.unique()
        if (values.size() > 1) {
            throw new IllegalArgumentException(
                "Sample '${sampleId}' has conflicting ${key} values across source libraries: ${values.join(', ')}"
            )
        }
        values ? values[0] : fallback
    }

    samplesheet
        .splitCsv(header: true)
        .map { row ->
            def sampleId = normalize(row.sample_id)
            def fastqR1 = normalize(row.fastq_r1)
            def fastqR2 = normalize(row.fastq_r2)
            def bamPath = normalize(row.bam)
            def matrixPath = normalize(row.matrix_path)
            def whitelistPath = normalize(row.barcode_whitelist)
            tuple(
                sampleId,
                [
                    sample_id        : sampleId,
                    source_library_id: normalize(row.library_id),
                    protocol         : normalize(row.protocol) ?: protocol_mode,
                    chemistry        : normalize(row.chemistry),
                    condition        : normalize(row.condition),
                    replicate_id     : normalize(row.replicate_id) ?: sampleId,
                    fastq_r1         : fastqR1 ? file(fastqR1, checkIfExists: true) : null,
                    fastq_r2         : fastqR2 ? file(fastqR2, checkIfExists: true) : null,
                    bam              : bamPath ? file(bamPath, checkIfExists: true) : nofile,
                    matrix           : matrixPath ? file(matrixPath, checkIfExists: true) : nofile,
                    whitelist        : whitelistPath ? file(whitelistPath, checkIfExists: true) : nofile,
                ]
            )
        }
        .groupTuple()
        .map { sampleId, entries ->
            def read1s = entries.collect { it.fastq_r1 }.findAll { it != null }
            def read2s = entries.collect { it.fastq_r2 }.findAll { it != null }
            if (read1s.size() != read2s.size()) {
                throw new IllegalArgumentException(
                    "Sample '${sampleId}' has mismatched FASTQ mates after sample-level grouping: R1=${read1s.size()} R2=${read2s.size()}"
                )
            }

            def bamInputs = entries.collect { it.bam }.findAll { it != null && it.name != 'NO_FILE' }
            def matrixInputs = entries.collect { it.matrix }.findAll { it != null && it.name != 'NO_FILE' }
            def whitelistCandidates = entries.collect { it.whitelist }.findAll { it != null && it.name != 'NO_FILE' }
            def whitelistPaths = whitelistCandidates.collect { it.toString() }.unique()
            if (whitelistPaths.size() > 1) {
                throw new IllegalArgumentException(
                    "Sample '${sampleId}' has multiple barcode whitelist files across source libraries: ${whitelistPaths.join(', ')}"
                )
            }

            def sourceLibraries = entries.collect { normalize(it.source_library_id) }.findAll { it }.unique().sort()
            def meta = [
                sample_id         : sampleId,
                library_id        : sampleId,
                source_library_ids: sourceLibraries.join(','),
                protocol          : distinctOrThrow(sampleId, entries, 'protocol', protocol_mode),
                chemistry         : distinctOrThrow(sampleId, entries, 'chemistry', ''),
                condition         : distinctOrThrow(sampleId, entries, 'condition', ''),
                replicate_id      : distinctOrThrow(sampleId, entries, 'replicate_id', sampleId),
            ]

            tuple(meta, read1s, read2s, bamInputs, matrixInputs, whitelistCandidates ? whitelistCandidates[0] : nofile)
        }
        .set { ch_rows }

    ch_rows
        .filter { meta, read1s, read2s, bamInputs, matrixInputs, whitelist ->
            run_mode == 'full' && read1s.size() >= 1 && read2s.size() >= 1
        }
        .combine(reference_bundle)
        .map { meta, read1s, read2s, bamInputs, matrixInputs, whitelist, reference_meta, star_index, gtf, fasta, chrom_sizes, terminal_exons, atlas, blacklist ->
            tuple(meta, read1s, read2s, whitelist, star_index, gtf, meta.protocol)
        }
        .set { ch_fastq }

    ch_rows
        .filter { meta, read1s, read2s, bamInputs, matrixInputs, whitelist ->
            bamInputs.size() >= 1 && (run_mode in ['from_bam', 'feedback'] || (read1s.size() == 0 && read2s.size() == 0))
        }
        .map { meta, read1s, read2s, bamInputs, matrixInputs, whitelist ->
            if (bamInputs.size() != 1) {
                throw new IllegalArgumentException(
                    "Sample-level BAM import currently expects exactly one BAM per sample; sample '${meta.sample_id}' received ${bamInputs.size()}"
                )
            }
            if (matrixInputs.size() > 1) {
                throw new IllegalArgumentException(
                    "Sample-level BAM import currently expects at most one matrix bundle per sample; sample '${meta.sample_id}' received ${matrixInputs.size()}"
                )
            }
            tuple(meta, bamInputs[0], matrixInputs ? matrixInputs[0] : nofile, whitelist)
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
