process SCAPTURE_PASMERGE {
    tag "${meta.library_id}:${category}"
    label 'process_medium'

    conda "${projectDir}/envs/scapture.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-scapture.sif" : null

    input:
    tuple val(meta), val(category), path(evaluated_bed)

    output:
    tuple val(meta), val(category), path("${meta.library_id}_merge.${category}.Integrated.bed"), emit: integrated
    tuple val(meta), val(category), path("${meta.library_id}_merge.${category}.status.tsv"), emit: status
    tuple val(meta), val(category), path("${meta.library_id}_merge.${category}.PASmerge.log"), emit: log
    path 'versions.yml', emit: versions

    script:
    if (!(category in ['exonic', 'intronic', '3primeExtended'])) {
        error "Invalid SCAPTURE PASmerge category: ${category}"
    }
    def prefix = meta.library_id
    """
    set -euo pipefail
    export TMPDIR=\$(pwd)/tmp
    mkdir -p "\$TMPDIR"
    test -e ${evaluated_bed}
    out=${prefix}_merge.${category}.Integrated.bed
    log=${prefix}_merge.${category}.PASmerge.log
    status=${prefix}_merge.${category}.status.tsv
    scored=\$(awk '\$2 ~ /^[0-9]+\$/ && \$3 ~ /^[0-9]+\$/ && NF >= 14 && \$14 != "NA" { n++ } END { print n+0 }' ${evaluated_bed})
    if [ "\$scored" -eq 0 ]; then
        : > "\$out"
        printf 'INFO: no DeepPASS-labeled rows in %s; category skipped\n' ${evaluated_bed} > "\$log"
        printf 'category\tstate\tscored_rows\tintegrated_rows\n%s\tempty_input\t0\t0\n' ${category} > "\$status"
    else
        printf '%s\t%s\n' ${prefix} ${evaluated_bed} > peaks_list.${category}.txt
        scapture -m PASmerge -o ${prefix}_merge.${category} --peak peaks_list.${category}.txt &> "\$log"
        if [ ! -s "\$out" ]; then
            echo 'ERROR: PASmerge produced no Integrated.bed for a scored category' >&2
            tail -80 "\$log" >&2 || true
            exit 1
        fi
        invalid=\$(awk 'NF < 12 || \$2 !~ /^[0-9]+\$/ || \$3 !~ /^[0-9]+\$/ { n++ } END { print n+0 }' "\$out")
        if [ "\$invalid" -ne 0 ]; then
            echo "ERROR: PASmerge produced \$invalid invalid BED12 rows" >&2
            exit 1
        fi
        printf 'category\tstate\tscored_rows\tintegrated_rows\n%s\tmerged\t%s\t%s\n' \
            ${category} "\$scored" "\$(wc -l < "\$out")" > "\$status"
    fi
    commit=\$(git -C /opt/scapture rev-parse HEAD 2>/dev/null || printf 'unknown')
    printf '"SCAPTURE_PASMERGE":\n  scapture_commit: "%s"\n  scapture_script_sha256: "%s"\n' \
        "\$commit" "\$(sha256sum "\$(command -v scapture)" | cut -d' ' -f1)" > versions.yml
    """
}
