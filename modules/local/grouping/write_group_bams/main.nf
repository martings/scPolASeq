process WRITE_GROUP_BAMS {
    tag "$meta.library_id"
    label 'process_single'
    label 'process_python'

    conda "${projectDir}/envs/python.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-python.sif" : null
    publishDir "${params.outdir}/coverage", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(group_map)

    output:
    tuple val(meta), path("${meta.library_id}.group_bams"), emit: group_bam_dir
    tuple val(meta), path("${meta.library_id}.group_bams_manifest.tsv"), emit: manifest

    script:
    """
    python3 ${projectDir}/bin/write_group_bams.py \\
        --group-map ${group_map} \\
        --sample-id ${meta.sample_id} \\
        --out-dir ${meta.library_id}.group_bams \\
        --out-manifest ${meta.library_id}.group_bams_manifest.tsv
    """

    stub:
    """
    mkdir -p ${meta.library_id}.group_bams
    touch ${meta.library_id}.group_bams/cluster__cluster_1.bam
    printf "group_level\\tgroup_id\\tplaceholder_bam\\ncluster\\tcluster_1\\tcluster__cluster_1.bam\\n" > ${meta.library_id}.group_bams_manifest.tsv
    """
}
