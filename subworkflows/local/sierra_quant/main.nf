include { SIERRA_QUANT_PREPARE } from '../../../modules/local/apa/sierra_quant_prepare/main'
include { SIERRA_COUNTPEAKS    } from '../../../modules/nf-core/sierra/countpeaks/main'
include { SIERRA_QUANT_TABLE   } from '../../../modules/local/apa/sierra_quant_table/main'

/*
 * Composite subworkflow adapting the atomic nf-core SIERRA_COUNTPEAKS module
 * to scPolASeq's group-scoped tabular PAS quantification contract.
 */
workflow SIERRA_QUANT {
    take:
    ch_sierra_input // tuple val(meta), val(group_level), val(group_id), path(grouped_bam), path(grouped_bai), path(pas_reference), path(cell_annotations), path(gtf)

    main:
    SIERRA_QUANT_PREPARE(ch_sierra_input)

    ch_prepared_keyed = SIERRA_QUANT_PREPARE.out.prepared.map { meta, group_level, group_id, whitelist, peak_sites, site_id_map, bam, bai, gtf, pas_ref, ann ->
        def safeGroupId    = group_id.toString().replaceAll(/[^A-Za-z0-9._-]+/, '_')
        def safeLibraryId  = meta.library_id.toString().replaceAll(/[^A-Za-z0-9._-]+/, '_')
        def safeGroupLevel = group_level.toString().replaceAll(/[^A-Za-z0-9._-]+/, '_')
        def prefix         = "${safeLibraryId}.${safeGroupLevel}.${safeGroupId}"
        def strand         = (meta.strandedness ?: (meta.chemistry == '10xv1' ? 'reverse' : 'forward')).toString().toLowerCase()
        def count_meta     = [
            id          : prefix,
            orig_meta   : meta,
            library_id  : meta.library_id,
            group_level : group_level,
            group_id    : group_id,
            strand      : strand,
            strandedness: meta.strandedness,
            chemistry   : meta.chemistry
        ]
        tuple(prefix, count_meta, whitelist, peak_sites, site_id_map, bam, bai, gtf, pas_ref, ann)
    }

    ch_count_inputs = ch_prepared_keyed.multiMap { prefix, count_meta, whitelist, peak_sites, site_id_map, bam, bai, gtf, pas_ref, ann ->
        bam_ch   : tuple(count_meta, bam, bai, whitelist)
        peaks_ch : peak_sites
        gtf_ch   : gtf
    }

    SIERRA_COUNTPEAKS(ch_count_inputs.bam_ch, ch_count_inputs.peaks_ch, ch_count_inputs.gtf_ch)

    ch_table_input = SIERRA_COUNTPEAKS.out.counts
        .map { count_meta, matrix_mtx, barcodes_tsv, sitenames_tsv ->
            tuple(count_meta.id, count_meta, matrix_mtx, barcodes_tsv, sitenames_tsv)
        }
        .join(
            SIERRA_COUNTPEAKS.out.log.map { count_meta, log_file -> tuple(count_meta.id, log_file) },
            failOnDuplicate: true,
            failOnMismatch: true
        )
        .join(
            ch_prepared_keyed.map { prefix, count_meta, whitelist, peak_sites, site_id_map, bam, bai, gtf, pas_ref, ann ->
                tuple(prefix, site_id_map, bam, bai, gtf, pas_ref, ann)
            },
            failOnDuplicate: true,
            failOnMismatch: true
        )
        .map { prefix, count_meta, matrix_mtx, barcodes_tsv, sitenames_tsv, log_file, site_id_map, bam, bai, gtf, pas_ref, ann ->
            tuple(count_meta.orig_meta, count_meta.group_level, count_meta.group_id, count_meta, matrix_mtx, barcodes_tsv, sitenames_tsv, log_file, site_id_map, bam, bai, gtf, pas_ref, ann)
        }

    SIERRA_QUANT_TABLE(ch_table_input)

    emit:
    quantified_pas      = SIERRA_QUANT_TABLE.out.quantified_pas
    manifest            = SIERRA_QUANT_TABLE.out.manifest
    quant_log           = SIERRA_QUANT_TABLE.out.log
    versions            = SIERRA_QUANT_TABLE.out.versions
    countpeaks_versions = SIERRA_COUNTPEAKS.out.versions
    prepare_versions    = SIERRA_QUANT_PREPARE.out.versions
}
