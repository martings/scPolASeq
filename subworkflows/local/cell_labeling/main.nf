include { SCANPY_CLUSTER_SC } from '../../../modules/local/clustering/scanpy_cluster_sc/main'
include { BUILD_GROUP_MAP   } from '../../../modules/local/grouping/build_group_map/main'
include { WRITE_MANIFESTS as WRITE_LABEL_MANIFESTS } from '../../../modules/local/provenance/write_manifests/main'

workflow CELL_LABELING_SC {
    take:
    matrix_bundle
    cell_metadata
    enable_internal_clustering
    grouping_levels

    main:
    matrix_bundle
        .combine(cell_metadata)
        .map { meta, matrix_dir, metadata_file ->
            tuple(meta, matrix_dir, metadata_file, enable_internal_clustering)
        }
        .set { ch_cluster_input }

    SCANPY_CLUSTER_SC(ch_cluster_input)
    BUILD_GROUP_MAP(SCANPY_CLUSTER_SC.out.cell_annotations.map { meta, tsv -> tsv }.collect(), grouping_levels)
    WRITE_LABEL_MANIFESTS(BUILD_GROUP_MAP.out.group_map, Channel.value('labels'))

    emit:
    cell_annotations = BUILD_GROUP_MAP.out.cell_annotations
    group_map        = BUILD_GROUP_MAP.out.group_map
    group_summary    = BUILD_GROUP_MAP.out.group_summary
    label_manifest   = WRITE_LABEL_MANIFESTS.out.manifest
    h5ad             = SCANPY_CLUSTER_SC.out.h5ad
}
