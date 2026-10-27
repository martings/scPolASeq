process SCAPTURE_PASREFERENCE {
    tag "${meta.library_id}"
    label 'process_low'

    conda "${projectDir}/envs/scapture.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-scapture.sif" : null

    input:
    tuple val(meta), path(category_beds)

    output:
    tuple val(meta), path("${meta.library_id}_merge.Integrated.bed"), emit: reference
    tuple val(meta), path("${meta.library_id}_merge.reference.metrics.tsv"), emit: metrics
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.library_id
    """
    set -euo pipefail
    raw=${prefix}_merge.Integrated.raw.bed
    out=${prefix}_merge.Integrated.bed
    : > "\$raw"
    printf 'category\tintegrated_rows\n' > ${prefix}_merge.reference.metrics.tsv
    nonempty=0
    for cat in exonic intronic 3primeExtended; do
        bed=${prefix}_merge.\${cat}.Integrated.bed
        if [ ! -e "\$bed" ]; then
            echo "ERROR: missing PASmerge category \$cat" >&2
            exit 1
        fi
        rows=\$(wc -l < "\$bed")
        printf '%s\t%s\n' "\$cat" "\$rows" >> ${prefix}_merge.reference.metrics.tsv
        if [ "\$rows" -gt 0 ]; then
            cat "\$bed" >> "\$raw"
            nonempty=\$((nonempty+1))
        fi
    done
    if [ "\$nonempty" -eq 0 ]; then
        echo 'ERROR: all SCAPTURE PASmerge categories were empty' >&2
        exit 1
    fi
    invalid=\$(awk 'NF < 12 || \$2 !~ /^[0-9]+\$/ || \$3 !~ /^[0-9]+\$/ { n++ } END { print n+0 }' "\$raw")
    if [ "\$invalid" -ne 0 ]; then
        echo "ERROR: combined PAS reference contains \$invalid invalid BED12 rows" >&2
        exit 1
    fi
    sort -k1,1 -k2,2n "\$raw" \
        | awk 'BEGIN{OFS="\t"} {print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12}' > "\$out"
    test -s "\$out"
    printf 'total\t%s\n' "\$(wc -l < "\$out")" >> ${prefix}_merge.reference.metrics.tsv
    printf '"SCAPTURE_PASREFERENCE":\n  coreutils: "%s"\n' "\$(sort --version | head -1 | awk '{print \$NF}')" > versions.yml
    """
}
