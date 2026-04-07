process GROUPED_3PRIME_COVERAGE {
    tag "$meta.library_id"
    label 'process_medium'
    label 'process_python'

    conda "${projectDir}/envs/python.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-python.sif" : null
    publishDir "${params.outdir}/coverage", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai), path(barcode_registry), path(group_map), path(chrom_sizes), val(emit_bigwigs)

    output:
    tuple val(meta), path("${meta.library_id}.grouped_3prime_counts.tsv"), emit: coverage_bundle
    tuple val(meta), path("${meta.library_id}.tracks"), emit: track_bundle
    tuple val(meta), path("${meta.library_id}.coverage_summary.tsv"), emit: coverage_summary

    script:
    """
    python ${projectDir}/bin/grouped_3prime_coverage.py \\
        --bam ${bam} \\
        --barcode-registry ${barcode_registry} \\
        --group-map ${group_map} \\
        --chrom-sizes ${chrom_sizes} \\
        --sample-id ${meta.sample_id} \\
        --library-id ${meta.library_id} \\
        --emit-bigwigs ${emit_bigwigs} \\
        --min-mapq 255 \\
        --out-counts ${meta.library_id}.grouped_3prime_counts.tsv \\
        --out-track-dir ${meta.library_id}.tracks \\
        --out-summary ${meta.library_id}.coverage_summary.tsv
    """

    stub:
    """
    printf "sample_id\\tlibrary_id\\tgroup_level\\tgroup_id\\tchrom\\tposition\\tstrand\\tread_count\\tumi_count\\n${meta.sample_id}\\t${meta.library_id}\\tcluster\\tcluster_1\\tchr1\\t100\\t+\\t1\\t1\\n" > ${meta.library_id}.grouped_3prime_counts.tsv
    mkdir -p ${meta.library_id}.tracks
    printf "chr1\\t99\\t100\\t1\\n" > ${meta.library_id}.tracks/cluster_1.bedGraph
    touch ${meta.library_id}.tracks/cluster_1.bigWig
    printf "sample_id\\tlibrary_id\\tmode\\tn_groups\\tn_loci\\n${meta.sample_id}\\t${meta.library_id}\\tstub\\t1\\t1\\n" > ${meta.library_id}.coverage_summary.tsv
    """
}
