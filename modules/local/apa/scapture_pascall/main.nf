process SCAPTURE_PASCALL {
    tag "${meta.library_id}"
    label 'process_high'

    conda "${projectDir}/envs/scapture.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-scapture.sif" : null
    // The patched predictor is a staged input, so changing its content changes
    // the PAScall task hash. The bind must point at this staged path.
    containerOptions { params.apptainer_cache_dir ? "--bind ${predict_script}:/opt/scapture/DeepPASS/Predict.py" : "" }

    input:
    tuple val(meta), path(bam), path(bai), path(genome_fasta), path(genome_fai), path(annotation_files), path(polya_db), val(has_polya_db), path(predict_script)

    output:
    tuple val(meta), path("${meta.library_id}.exonic.peaks.evaluated.bed"), path("${meta.library_id}.intronic.peaks.evaluated.bed"), path("${meta.library_id}.3primeExtended.peaks.evaluated.bed"), emit: evaluated
    tuple val(meta), path("${meta.library_id}.PAScall.log"), emit: log
    tuple val(meta), path("${meta.library_id}.scapture.pascall.metrics.tsv"), emit: metrics
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.library_id
    def read_len = params.scapture_read_length ?: 90
    def species = params.scapture_species ?: 'human'
    def peak_width = params.scapture_peak_width ?: 400
    def db_arg = has_polya_db ? "--polyaDB ${polya_db}" : ''
    """
    set -euo pipefail
    export TMPDIR=\$(pwd)/tmp
    export MPLCONFIGDIR=\$(pwd)/tmp/mpl
    mkdir -p "\$TMPDIR" "\$MPLCONFIGDIR"
    test -s ${bam}
    test -s ${bai}
    test -s ${genome_fasta}
    test -s ${genome_fai}
    test -s ${prefix}_annot.genetype.gtf
    test -s ${prefix}_annot.genetype.bed
    if ${has_polya_db}; then test -s ${polya_db}; fi

    ln -sfn "\$(readlink -f ${genome_fasta})" ${prefix}.scapture.genome.fa
    ln -sfn "\$(readlink -f ${genome_fai})" ${prefix}.scapture.genome.fa.fai
    scapture -m PAScall \
        -a ${prefix}_annot \
        -g ${prefix}.scapture.genome.fa \
        -b ${bam} \
        -l ${read_len} \
        -o ${prefix} \
        -p ${task.cpus ?: 1} \
        --species ${species} \
        -w ${peak_width} \
        ${db_arg} \
        &> ${prefix}.PAScall.log

    # Preserve the monolith's recovery of category files from per-gene output.
    per_gene=${prefix}_tmp/Per_Gene
    for cat in exonic intronic 3primeExtended; do
        bed=${prefix}.\${cat}.peaks.evaluated.bed
        if [ ! -s "\$bed" ]; then
            : > "\$bed"
            case "\$cat" in
                exonic) pattern='scapture.*.exonic.peak.isoform.normality.aggregate.bed' ;;
                intronic) pattern='scapture.*.intronic.peak.aggregate.bed' ;;
                *) pattern='' ;;
            esac
            if [ -n "\$pattern" ] && [ -d "\$per_gene" ]; then
                find "\$per_gene" -maxdepth 1 -type f -name "\$pattern" -exec cat {} + \
                    | sort -k1,1 -k2,2n > "\$bed"
            fi
        fi
    done

    # DeepPASS evaluated BED column 14 is a positive/negative prediction label.
    # It is not a numeric probability; NA denotes evaluation failure.
    scored_rows() {
        awk '\$2 ~ /^[0-9]+\$/ && \$3 ~ /^[0-9]+\$/ && NF >= 14 && \$14 != "NA" { n++ } END { print n+0 }' "\$1"
    }
    valid_rows() {
        awk '\$2 ~ /^[0-9]+\$/ && \$3 ~ /^[0-9]+\$/ { n++ } END { print n+0 }' "\$1"
    }
    valid=0
    scored=0
    printf 'metric\tvalue\n' > ${prefix}.scapture.pascall.metrics.tsv
    for cat in exonic intronic 3primeExtended; do
        bed=${prefix}.\${cat}.peaks.evaluated.bed
        nvalid=\$(valid_rows "\$bed")
        nscored=\$(scored_rows "\$bed")
        valid=\$((valid+nvalid))
        scored=\$((scored+nscored))
        printf '%s_valid\t%s\n%s_scored\t%s\n' "\$cat" "\$nvalid" "\$cat" "\$nscored" >> ${prefix}.scapture.pascall.metrics.tsv
    done
    printf 'valid_peaks\t%s\nscored_peaks\t%s\n' "\$valid" "\$scored" >> ${prefix}.scapture.pascall.metrics.tsv
    if [ "\$valid" -eq 0 ] || [ "\$scored" -eq 0 ]; then
        echo 'ERROR: PAScall/DeepPASS produced no usable scored peaks' >&2
        tail -80 ${prefix}.PAScall.log >&2 || true
        exit 1
    fi
    if [ ! -s ${prefix}.3primeExtended.peaks.evaluated.bed ]; then
        echo 'INFO: 3primeExtended category absent from PAScall output' >> ${prefix}.PAScall.log
    fi
    commit=\$(git -C /opt/scapture rev-parse HEAD 2>/dev/null || printf 'unknown')
    printf '"SCAPTURE_PASCALL":\n  scapture_commit: "%s"\n  scapture_script_sha256: "%s"\n' \
        "\$commit" "\$(sha256sum "\$(command -v scapture)" | cut -d' ' -f1)" > versions.yml
    """
}
