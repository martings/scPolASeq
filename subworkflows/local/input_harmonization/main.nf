include { VALIDATE_INPUTS         } from '../../../modules/local/input/validate_inputs/main'
include { HARMONIZE_CELL_METADATA } from '../../../modules/local/metadata/harmonize_cell_metadata/main'

workflow INPUT_HARMONIZATION {
    take:
    input_sheet
    cell_metadata
    cluster_assignments
    cell_type_labels

    main:
    VALIDATE_INPUTS(input_sheet, file("${projectDir}/assets/schema_input.json"))
    HARMONIZE_CELL_METADATA(cell_metadata, cluster_assignments, cell_type_labels)

    emit:
    samplesheet       = VALIDATE_INPUTS.out.samplesheet
    input_manifest    = VALIDATE_INPUTS.out.manifest
    cell_metadata     = HARMONIZE_CELL_METADATA.out.cell_metadata
    metadata_manifest = HARMONIZE_CELL_METADATA.out.manifest
}
