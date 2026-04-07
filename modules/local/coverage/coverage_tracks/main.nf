// Stage 5 — Strand-aware coverage tracks (bedGraph + bigWig) per grouped BAM.
// Separates forward and reverse strand for 3' APA inference.
process COVERAGE_TRACKS {
    tag "${group_id}:${group_level}"
    label 'process_medium'

    conda "${projectDir}/envs/alignment.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-alignment.sif" : null
    publishDir "${params.outdir}/coverage/${group_level}", mode: params.publish_dir_mode, pattern: "*.bedGraph"
    publishDir "${params.outdir}/coverage/${group_level}", mode: params.publish_dir_mode, pattern: "*.bw"

    input:
    tuple val(meta), val(group_level), val(group_id), path(bam), path(bai)
    path chrom_sizes

    output:
    tuple val(meta), val(group_level), val(group_id), path("${group_id}.fwd.bedGraph"), path("${group_id}.rev.bedGraph"), emit: bedgraphs
    tuple val(meta), val(group_level), val(group_id), path("${group_id}.fwd.bw"),       path("${group_id}.rev.bw"),       emit: bigwigs, optional: true
    path "${group_id}.coverage_stats.tsv", emit: coverage_stats

    script:
    """
    # Strand-aware coverage using bedtools
    if command -v bedtools >/dev/null 2>&1 && [ -s ${bam} ]; then
        bedtools genomecov -ibam ${bam} -bga -strand + -split \\
            | sort -k1,1 -k2,2n > ${group_id}.fwd.bedGraph || touch ${group_id}.fwd.bedGraph
        bedtools genomecov -ibam ${bam} -bga -strand - -split \\
            | sort -k1,1 -k2,2n > ${group_id}.rev.bedGraph || touch ${group_id}.rev.bedGraph
    else
        touch ${group_id}.fwd.bedGraph ${group_id}.rev.bedGraph
    fi

    # Convert to bigWig if tool available and chrom_sizes is real
    if command -v bedGraphToBigWig >/dev/null 2>&1 && [ -s ${group_id}.fwd.bedGraph ] && [ -s ${chrom_sizes} ]; then
        bedGraphToBigWig ${group_id}.fwd.bedGraph ${chrom_sizes} ${group_id}.fwd.bw || touch ${group_id}.fwd.bw
        bedGraphToBigWig ${group_id}.rev.bedGraph ${chrom_sizes} ${group_id}.rev.bw || touch ${group_id}.rev.bw
    else
        touch ${group_id}.fwd.bw ${group_id}.rev.bw
    fi

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
    touch ${group_id}.fwd.bw ${group_id}.rev.bw
    printf "group_level\tgroup_id\tstrand\ttotal_coverage\n" > ${group_id}.coverage_stats.tsv
    printf "${group_level}\t${group_id}\t+\t1300\n${group_level}\t${group_id}\t-\t700\n" \
        >> ${group_id}.coverage_stats.tsv
    """
}
