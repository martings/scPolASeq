/*
 * Standalone entrypoint for SCAPTURE staged subworkflow.
 *
 * Runs the opt-in, cache-separated SCAPTURE stages:
 *   SCAPTURE_PREPARE -> SCAPTURE_ANNOTATION -> SCAPTURE_PASCALL (GPU) ->
 *   SCAPTURE_MERGE (PASmerge + PASreference) -> SCAPTURE_CATALOG ->
 *   SCAPTURE_PASQUANT -> SCAPTURE_MANIFEST
 *
 * Input modes:
 *   1. Single sample:
 *      --bam <file.bam> --cell_annotations <file.tsv> [--library_id <id>] [--sample_id <id>] [--bai <file.bai>]
 *   2. Samplesheet CSV:
 *      --input <samplesheet.csv> (columns: sample_id,library_id,bam,bai,cell_annotations)
 */
nextflow.enable.dsl = 2

include { SCAPTURE_STAGED } from './subworkflows/local/scapture_staged/main'

params.input                       = null
params.library_id                  = null
params.sample_id                   = null
params.bam                         = null
params.bai                         = null
params.cell_annotations           = null
params.genome_fasta                = null
params.genome_fai                  = null
params.gtf                         = null
params.chrom_sizes                 = null
params.scapture_polya_db           = null
params.scapture_species            = 'human'
params.scapture_read_length        = 90
params.scapture_peak_width         = 400
params.scapture_normality_seed      = 1
params.scapture_deeppass_threshold  = 0.5
params.outdir                      = 'results'
params.publish_dir_mode            = 'copy'
params.apptainer_cache_dir         = null

workflow {
    if (!params.genome_fasta) {
        error "Missing required parameter: --genome_fasta"
    }
    if (!params.gtf) {
        error "Missing required parameter: --gtf"
    }

    def fasta = file(params.genome_fasta, checkIfExists: true)
    def fai = params.genome_fai ? file(params.genome_fai, checkIfExists: true) : file("${params.genome_fasta}.fai", checkIfExists: true)
    def gtf = file(params.gtf, checkIfExists: true)

    def chrom_sizes_path = params.chrom_sizes ?: "${fasta.parent}/chrom.sizes"
    def chrom_sizes = file(chrom_sizes_path, checkIfExists: true)

    def nofile = file("${projectDir}/assets/NO_FILE")
    def polya_db = params.scapture_polya_db ? file(params.scapture_polya_db, checkIfExists: true) : nofile
    def has_polya_db = params.scapture_polya_db ? true : false

    def ch_input

    if (params.input) {
        def sheet = file(params.input, checkIfExists: true)
        ch_input = Channel.fromPath(sheet)
            .splitCsv(header: true)
            .map { row ->
                def lib_id = (row.library_id ?: row.sample_id)?.trim()
                def smp_id = (row.sample_id ?: lib_id)?.trim()
                if (!lib_id) error "A row in ${params.input} is missing library_id/sample_id: ${row}"
                def meta = [library_id: lib_id, sample_id: smp_id]
                def bam_p = row.bam?.trim()
                if (!bam_p) error "Missing BAM path for ${lib_id} in ${params.input}"
                def bai_p = row.bai?.trim() ?: "${bam_p}.bai"
                def annot_p = row.cell_annotations?.trim()
                if (!annot_p) error "Missing cell_annotations path for ${lib_id} in ${params.input}"
                tuple(
                    meta,
                    file(bam_p, checkIfExists: true),
                    file(bai_p, checkIfExists: true),
                    fasta,
                    fai,
                    gtf,
                    chrom_sizes,
                    file(annot_p, checkIfExists: true),
                    polya_db,
                    has_polya_db
                )
            }
    } else if (params.bam && params.cell_annotations) {
        def bam_file = file(params.bam, checkIfExists: true)
        def bai_file = params.bai ? file(params.bai, checkIfExists: true) : file("${params.bam}.bai", checkIfExists: true)
        def annot_file = file(params.cell_annotations, checkIfExists: true)
        def lib_id = (params.library_id ?: params.sample_id ?: bam_file.name.replaceAll(/\.(filtered|Aligned).*/, ''))?.trim()
        def smp_id = (params.sample_id ?: lib_id)?.trim()
        def meta = [library_id: lib_id, sample_id: smp_id]
        ch_input = Channel.of(tuple(
            meta,
            bam_file,
            bai_file,
            fasta,
            fai,
            gtf,
            chrom_sizes,
            annot_file,
            polya_db,
            has_polya_db
        ))
    } else {
        error "Supply either --input <samplesheet.csv> or (--bam <file.bam> and --cell_annotations <file.tsv>)"
    }

    SCAPTURE_STAGED(ch_input)
}
