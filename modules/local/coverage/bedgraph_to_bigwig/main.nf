// Stage 5b — Optional BigWig conversion from strand-aware bedGraph tracks.
// Fail-closed when enabled: never emits empty placeholder BigWig files.
process BEDGRAPH_TO_BIGWIG {
    tag "${group_id}:${group_level}"
    label 'process_single'

    conda "${projectDir}/envs/python.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-python.sif" : null
    publishDir "${params.outdir}/coverage/${group_level}", mode: params.publish_dir_mode, pattern: "*.bw"
    publishDir "${params.outdir}/coverage/${group_level}", mode: params.publish_dir_mode, pattern: "*.bigwig_status.tsv"

    input:
    tuple val(meta), val(group_level), val(group_id), path(fwd_bedgraph), path(rev_bedgraph)
    path chrom_sizes

    output:
    tuple val(meta), val(group_level), val(group_id), path("${group_id}.fwd.bw"), path("${group_id}.rev.bw"), emit: bigwigs
    path "${group_id}.bigwig_status.tsv", emit: status
    path "${group_id}.bigwig.versions.yml", emit: versions

    script:
    """
    command -v bedGraphToBigWig >/dev/null 2>&1 ||
        { echo 'bedGraphToBigWig executable unavailable; install ucsc-bedgraphtobigwig' >&2; exit 127; }
    test -s ${chrom_sizes} ||
        { echo 'Chromosome sizes file is missing or empty' >&2; exit 1; }
    test -s ${fwd_bedgraph} && test -s ${rev_bedgraph} ||
        { echo 'Cannot create BigWig from an empty bedGraph' >&2; exit 1; }

    bedGraphToBigWig ${fwd_bedgraph} ${chrom_sizes} ${group_id}.fwd.bw ||
        { echo 'Forward bedGraphToBigWig conversion failed' >&2; exit 1; }
    test -s ${group_id}.fwd.bw ||
        { echo 'Forward BigWig output is missing or empty' >&2; exit 1; }

    bedGraphToBigWig ${rev_bedgraph} ${chrom_sizes} ${group_id}.rev.bw ||
        { echo 'Reverse bedGraphToBigWig conversion failed' >&2; exit 1; }
    test -s ${group_id}.rev.bw ||
        { echo 'Reverse BigWig output is missing or empty' >&2; exit 1; }

    fwd_bytes=\$(wc -c < ${group_id}.fwd.bw | tr -d '[:space:]')
    rev_bytes=\$(wc -c < ${group_id}.rev.bw | tr -d '[:space:]')
    converter_ver=\$(bedGraphToBigWig 2>&1 | awk 'NF {print; exit}' || true)
    converter_ver=\${converter_ver:-unknown}

    printf "group_level\\tgroup_id\\tstatus\\tfwd_bw_bytes\\trev_bw_bytes\\tconverter\\n" > ${group_id}.bigwig_status.tsv
    printf "%s\\t%s\\tok\\t%s\\t%s\\t%s\\n" "${group_level}" "${group_id}" "\$fwd_bytes" "\$rev_bytes" "\$converter_ver" >> ${group_id}.bigwig_status.tsv

    cat <<-END_VERSIONS > ${group_id}.bigwig.versions.yml
    "${task.process}":
      bedGraphToBigWig: "\$converter_ver"
    END_VERSIONS
    """

    stub:
    """
    touch ${group_id}.fwd.bw ${group_id}.rev.bw
    printf "group_level\\tgroup_id\\tstatus\\tfwd_bw_bytes\\trev_bw_bytes\\tconverter\\n" > ${group_id}.bigwig_status.tsv
    printf "${group_level}\\t${group_id}\\tstub\\t0\\t0\\tstub\\n" >> ${group_id}.bigwig_status.tsv
    cat <<-END_VERSIONS > ${group_id}.bigwig.versions.yml
    "${task.process}":
      bedGraphToBigWig: "stub"
    END_VERSIONS
    """
}
