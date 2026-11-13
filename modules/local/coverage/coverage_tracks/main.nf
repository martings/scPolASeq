// Stage 5 — Strand-aware coverage tracks (bedGraph + bigWig) per grouped BAM.
// Separates forward and reverse strand for 3' APA inference.
process COVERAGE_TRACKS {
    tag "${group_id}:${group_level}"
    label 'process_medium'

    conda "${projectDir}/envs/python.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-python.sif" : null
    publishDir "${params.outdir}/coverage/${group_level}", mode: params.publish_dir_mode, pattern: "*.bedGraph"

    input:
    tuple val(meta), val(group_level), val(group_id), path(bam), path(bai)
    path coverage_script

    output:
    tuple val(meta), val(group_level), val(group_id), path("${group_id}.fwd.bedGraph"), path("${group_id}.rev.bedGraph"), emit: bedgraphs
    path "${group_id}.coverage_stats.tsv", emit: coverage_stats

    script:
    """
    # Strand-aware coverage using pysam (bam_to_bedgraph.py) — single pass for both strands.
    test -s ${bam} || { echo 'Coverage input BAM is missing or empty' >&2; exit 1; }
    python3 ${coverage_script} ${bam} \\
        --fwd-output ${group_id}.fwd.bedGraph \\
        --rev-output ${group_id}.rev.bedGraph \\
        --cpus ${task.cpus} || { echo 'Coverage bedGraph generation failed' >&2; exit 1; }
    test -f ${group_id}.fwd.bedGraph && test -f ${group_id}.rev.bedGraph ||
        { echo 'Coverage bedGraph output is missing' >&2; exit 1; }

    # Coverage stats
    printf "group_level\tgroup_id\tstrand\ttotal_coverage\n" > ${group_id}.coverage_stats.tsv
    awk 'BEGIN{s=0}{s+=(\$4*(\$3-\$2))}END{printf "${group_level}\t${group_id}\t+\t%d\n",s}' \\
        ${group_id}.fwd.bedGraph >> ${group_id}.coverage_stats.tsv 2>/dev/null || true
    awk 'BEGIN{s=0}{s+=(\$4*(\$3-\$2))}END{printf "${group_level}\t${group_id}\t-\t%d\n",s}' \\
        ${group_id}.rev.bedGraph >> ${group_id}.coverage_stats.tsv 2>/dev/null || true
    """

    stub:
    """
    printf "chr1\t100\t200\t5\nchr1\t300\t400\t8\n" > ${group_id}.fwd.bedGraph
    printf "chr1\t100\t200\t3\nchr1\t300\t400\t4\n" > ${group_id}.rev.bedGraph
    printf "group_level\tgroup_id\tstrand\ttotal_coverage\n" > ${group_id}.coverage_stats.tsv
    printf "${group_level}\t${group_id}\t+\t1300\n${group_level}\t${group_id}\t-\t700\n" \
        >> ${group_id}.coverage_stats.tsv
    """
}
