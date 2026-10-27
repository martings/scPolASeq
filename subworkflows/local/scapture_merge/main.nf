include { SCAPTURE_PASMERGE } from '../../../modules/local/apa/scapture_pasmerge/main'
include { SCAPTURE_PASREFERENCE } from '../../../modules/local/apa/scapture_pasreference/main'

workflow SCAPTURE_MERGE {
    take:
    ch_evaluated

    main:
    ch_categories = ch_evaluated.flatMap { meta, exonic, intronic, extended ->
        [tuple(meta, 'exonic', exonic),
         tuple(meta, 'intronic', intronic),
         tuple(meta, '3primeExtended', extended)]
    }
    SCAPTURE_PASMERGE(ch_categories)

    ch_reference_input = SCAPTURE_PASMERGE.out.integrated
        .map { meta, category, bed -> tuple(meta.library_id, meta, category, bed) }
        .groupTuple(by: 0, size: 3)
        .map { library_id, metas, categories, beds ->
            if (categories.toSet() != ['exonic', 'intronic', '3primeExtended'].toSet()) {
                error "Incomplete PASmerge categories for ${library_id}: ${categories}"
            }
            tuple(metas[0], beds)
        }
    SCAPTURE_PASREFERENCE(ch_reference_input)

    emit:
    reference = SCAPTURE_PASREFERENCE.out.reference
    reference_metrics = SCAPTURE_PASREFERENCE.out.metrics
    category_integrated = SCAPTURE_PASMERGE.out.integrated
    category_status = SCAPTURE_PASMERGE.out.status
    category_logs = SCAPTURE_PASMERGE.out.log
    merge_versions = SCAPTURE_PASMERGE.out.versions
    reference_versions = SCAPTURE_PASREFERENCE.out.versions
}
