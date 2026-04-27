/*
 * SCAPTURE de novo PAS detection.
 *
 * Runs all three SCAPTURE modes (annotation -> PAScall -> PASmerge -> PASquant)
 * on the filtered library BAM and emits:
 *   - normalized site catalog TSV (converted from evaluated BED15 peaks)
 *   - per-cell UMI count matrix (PASquant, gzip-compressed)
 *   - manifest TSV
 *   - per-mode log files
 *
 * Reference: https://github.com/YangLab/SCAPTURE
 *
 * NOTE on PASmerge:
 *   SCAPTURE's PASmerge module is designed to merge ONE peak category at a time
 *   (exonic, intronic, or 3primeExtended). Each invocation produces a single
 *   <prefix>.Integrated.bed file. We therefore run PASmerge per-category and
 *   concatenate the resulting Integrated.bed files into a combined PAS reference
 *   for PASquant, mirroring the README's worked example.
 */
process SCAPTURE_FILTER {
    tag "${meta.library_id}"
    label 'process_high'

    conda "${projectDir}/envs/scapture.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-scapture.sif" : null
    // Bind-mount the patched Predict.py to fix TF 2.4 + pandas 1.1.5 incompatibility.
    // tf.convert_to_tensor cannot handle a pandas object-dtype Series of numpy arrays.
    containerOptions params.apptainer_cache_dir ? "--bind ${projectDir}/bin/DeepPASS_Predict.py:/opt/scapture/DeepPASS/Predict.py" : ""

    publishDir "${params.outdir}/scapture", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai), path(genome_fasta), path(genome_fai), path(gtf), path(chrom_sizes), path(cell_annotations), path(polyaDB), path(resume_pascall_dir)

    output:
    tuple val(meta), path("${meta.library_id}.scapture.site_catalog.tsv"), emit: site_catalog
    tuple val(meta), path("*.PASquant.KeepCell.UMIs.tsv.gz"),              emit: quant,    optional: true
    path "${meta.library_id}.scapture.manifest.tsv",                        emit: manifest
    path "${meta.library_id}.scapture.metrics.tsv",                         emit: metrics
    path "${meta.library_id}.scapture.gene_name_map.tsv",                   emit: gene_name_map
    path "${meta.library_id}.annotation.duplicate_gene_names.tsv",          emit: duplicate_gene_names, optional: true
    path "*.log",                                                            emit: logs

    script:
    def prefix     = "${meta.library_id}"
    def ncores     = task.cpus ?: 1
    def read_len   = params.scapture_read_length        ?: 90
    def species    = params.scapture_species             ?: 'human'
    def peak_width = params.scapture_peak_width          ?: 400
    def score_thr  = params.scapture_deeppass_threshold  ?: 0.5
    def db_arg     = (polyaDB.name != 'NO_FILE') ? "--polyaDB ${polyaDB}" : ""
    def has_resume = resume_pascall_dir.name != 'NO_FILE'
    """
    set -euo pipefail

    export TMPDIR=\$(pwd)/tmp
    export MPLCONFIGDIR=\$(pwd)/tmp/mpl
    mkdir -p "\$TMPDIR" "\$MPLCONFIGDIR"

    # FilterBamByTag is a Drop-seq tool missing from the container.
    # Copy our pysam replacement into the workdir (always bind-mounted) and
    # prepend it to PATH so scapture_quant's subprocess call resolves it.
    cp ${projectDir}/bin/FilterBamByTag ./FilterBamByTag
    chmod +x ./FilterBamByTag
    export PATH=\$(pwd):\$PATH

    if [ ! -s "${genome_fasta}" ]; then
        echo "ERROR: genome FASTA is not readable inside SCAPTURE task: ${genome_fasta}" >&2
        echo "If this is a symlink, make sure its target is mounted in the container." >&2
        exit 1
    fi

    if [ ! -s "${genome_fai}" ]; then
        echo "ERROR: genome FASTA index is missing or empty: ${genome_fai}" >&2
        exit 1
    fi

    if [ "${genome_fai}" != "${genome_fasta}.fai" ]; then
        cp -L "${genome_fai}" "${genome_fasta}.fai"
    fi

    valid_rows() {
        # \$1: BED file. Echoes rows whose start/end columns are integers.
        awk 'BEGIN{n=0} \$2 ~ /^[0-9]+\$/ && \$3 ~ /^[0-9]+\$/ {n++} END{print n+0}' "\$1"
    }

    scored_rows() {
        # SCAPTURE evaluated.bed col 14 is "positive" or "negative" (binary label).
        # Col 15 is the 200 bp genomic sequence — neither column is numeric.
        # Count rows where col 14 is a real prediction (not "NA").
        awk '
            BEGIN { n=0 }
            \$2 ~ /^[0-9]+\$/ && \$3 ~ /^[0-9]+\$/ && NF >= 14 && \$14 != "NA" { n++ }
            END { print n+0 }
        ' "\$1"
    }

    catalog_data_rows() {
        awk 'END{print NR > 0 ? NR - 1 : 0}' "\$1"
    }

    # ---------------------------------------------------------------------
    # 1. Extract cell barcodes for this library/sample
    # ---------------------------------------------------------------------
    awk -v sid="${meta.sample_id ?: meta.library_id}" '
        NR==1 {
            for(i=1;i<=NF;i++) {
                if(\$i=="barcode_corrected") bc=i
                else if(\$i=="barcode_raw" && !bc) bc=i
                if(\$i=="sample_id") si=i
            }
            if (!bc) {
                print "ERROR: no barcode_corrected or barcode_raw column" > "/dev/stderr"
                exit 1
            }
            next
        }
        si && \$si != sid { next }
        { print \$bc }
    ' ${cell_annotations} | sort -u > barcodes.txt

    if [ ! -s barcodes.txt ]; then
        echo "ERROR: no barcodes extracted for ${meta.sample_id ?: meta.library_id}" >&2
        exit 1
    fi

    if ${has_resume}; then
        # -----------------------------------------------------------------
        # 2R. Rescue mode: reuse PAScall peak annotation and rerun DeepPASS.
        # -----------------------------------------------------------------
        RESUME_DIR="${resume_pascall_dir}"
        echo "SCAPTURE rescue mode: reusing PAScall annotated peaks from \$RESUME_DIR" \\
            > ${prefix}.annotation.log
        : > ${prefix}.PAScall.log

        if [ -f "\$RESUME_DIR/${prefix}.PAScall.log" ]; then
            {
                echo "----- original PAScall log from rescue directory -----"
                cat "\$RESUME_DIR/${prefix}.PAScall.log"
                echo
                echo "----- rerun DeepPASS evaluation only -----"
            } >> ${prefix}.PAScall.log
        fi

        if [ -f "\$RESUME_DIR/${prefix}.scapture.gene_name_map.tsv" ]; then
            cp -L "\$RESUME_DIR/${prefix}.scapture.gene_name_map.tsv" \\
                ${prefix}.scapture.gene_name_map.tsv
        else
            printf "gene_id\\tscapture_gene_name\\toriginal_gene_name\\n" \\
                > ${prefix}.scapture.gene_name_map.tsv
        fi

        if [ -f "\$RESUME_DIR/${prefix}.annotation.duplicate_gene_names.tsv" ]; then
            cp -L "\$RESUME_DIR/${prefix}.annotation.duplicate_gene_names.tsv" \\
                ${prefix}.annotation.duplicate_gene_names.tsv
        else
            printf "gene_count\\toriginal_gene_name\\n" \\
                > ${prefix}.annotation.duplicate_gene_names.tsv
        fi

        SCAPTUREPATH="\$(dirname "\$(command -v scapture)")/"
        for cat in exonic intronic 3primeExtended; do
            src="\$RESUME_DIR/${prefix}.\${cat}.peaks.annotated.bed"
            inbed="${prefix}.\${cat}.peaks.annotated.bed"
            outprefix="${prefix}.\${cat}.peaks"
            if [ -e "\$src" ]; then
                cp -L "\$src" "\$inbed"
            else
                : > "\$inbed"
            fi

            if [ "\$(valid_rows "\$inbed")" -gt 0 ]; then
                scapture_evaluate \\
                    --peak "\$inbed" \\
                    -o "\$outprefix" \\
                    --species ${species} \\
                    -g ${genome_fasta} \\
                    ${db_arg} \\
                    --path "\$SCAPTUREPATH" \\
                    &>> ${prefix}.PAScall.log
            else
                echo "INFO: skipping DeepPASS for \${cat}; no valid annotated peak rows" \\
                    >> ${prefix}.PAScall.log
                : > "\$outprefix.evaluated.bed"
            fi
        done
    else
        # -----------------------------------------------------------------
        # 2. Prepare a SCAPTURE-specific GTF with collision-free gene tokens.
        # -----------------------------------------------------------------
        prepare_scapture_gtf.py \\
            --gtf ${gtf} \\
            --out-gtf ${prefix}.scapture.annotation.gtf \\
            --out-map ${prefix}.scapture.gene_name_map.tsv \\
            --out-duplicate-report ${prefix}.annotation.duplicate_gene_names.tsv

        # -----------------------------------------------------------------
        # 3. SCAPTURE annotation
        # -----------------------------------------------------------------
        scapture -m annotation \\
            -o ${prefix}_annot \\
            -g ${genome_fasta} \\
            --gtf ${prefix}.scapture.annotation.gtf \\
            --cs ${chrom_sizes} \\
            --extend 2000 \\
            &> ${prefix}.annotation.log

        if [ -s ${prefix}.annotation.duplicate_gene_names.tsv ]; then
            {
                echo
                echo "SCAPTURE duplicate original gene_name values detected and normalized."
                echo "Per-gene SCAPTURE jobs use gene_id-backed names from"
                echo "${prefix}.scapture.annotation.gtf to avoid filename collisions."
                echo
                head -20 ${prefix}.annotation.duplicate_gene_names.tsv
            } >> ${prefix}.annotation.log
        fi

        # -----------------------------------------------------------------
        # 4. SCAPTURE PAScall (peak calling + annotation + DeepPASS evaluation)
        # -----------------------------------------------------------------
        scapture -m PAScall \\
            -a ${prefix}_annot \\
            -g ${genome_fasta} \\
            -b ${bam} \\
            -l ${read_len} \\
            -o ${prefix} \\
            -p ${ncores} \\
            --species ${species} \\
            -w ${peak_width} \\
            ${db_arg} \\
            &> ${prefix}.PAScall.log
    fi

    # ---------------------------------------------------------------------
    # 5. Reconstruct evaluated BEDs if PAScall produced only per-gene outputs.
    #
    #   Global outputs that SCAPTURE should produce:
    #     ${prefix}.exonic.peaks.evaluated.bed
    #     ${prefix}.intronic.peaks.evaluated.bed
    #     ${prefix}.3primeExtended.peaks.evaluated.bed
    #
    #   If absent, exonic and intronic can be rebuilt from per-gene aggregates
    #   under ${prefix}_tmp/Per_Gene/. The 3primeExtended category is derived
    #   globally inside PAScall (no per-gene file) and cannot be reconstructed
    #   here — we fall back to whatever PAScall left in place, or an empty file.
    # ---------------------------------------------------------------------
    PER_GENE="${prefix}_tmp/Per_Gene"

    # exonic
    if [ ! -s ${prefix}.exonic.peaks.evaluated.bed ]; then
        : > ${prefix}.exonic.peaks.evaluated.bed
        if find "\$PER_GENE" -maxdepth 1 -type f \\
            -name "scapture.*.exonic.peak.isoform.normality.aggregate.bed" | grep -q .; then
            find "\$PER_GENE" -maxdepth 1 -type f \\
                -name "scapture.*.exonic.peak.isoform.normality.aggregate.bed" \\
                -exec cat {} + \\
              | sort -k1,1 -k2,2n > ${prefix}.exonic.peaks.evaluated.bed
        fi
    fi

    # intronic
    if [ ! -s ${prefix}.intronic.peaks.evaluated.bed ]; then
        : > ${prefix}.intronic.peaks.evaluated.bed
        if find "\$PER_GENE" -maxdepth 1 -type f \\
            -name "scapture.*.intronic.peak.aggregate.bed" | grep -q .; then
            find "\$PER_GENE" -maxdepth 1 -type f \\
                -name "scapture.*.intronic.peak.aggregate.bed" \\
                -exec cat {} + \\
              | sort -k1,1 -k2,2n > ${prefix}.intronic.peaks.evaluated.bed
        fi
    fi

    # 3primeExtended: produced globally, no per-gene fallback available
    if [ ! -e ${prefix}.3primeExtended.peaks.evaluated.bed ]; then
        : > ${prefix}.3primeExtended.peaks.evaluated.bed
    fi

    # ---------------------------------------------------------------------
    # 6. Metrics for debugging
    # ---------------------------------------------------------------------
    {
        echo -e "metric\\tvalue"
        echo -e "barcodes\\t\$(wc -l < barcodes.txt)"
        echo -e "exonic_evaluated_rows\\t\$(wc -l < ${prefix}.exonic.peaks.evaluated.bed)"
        echo -e "intronic_evaluated_rows\\t\$(wc -l < ${prefix}.intronic.peaks.evaluated.bed)"
        echo -e "extended_evaluated_rows\\t\$(wc -l < ${prefix}.3primeExtended.peaks.evaluated.bed)"
        echo -e "exonic_scored_rows\\t\$(scored_rows ${prefix}.exonic.peaks.evaluated.bed)"
        echo -e "intronic_scored_rows\\t\$(scored_rows ${prefix}.intronic.peaks.evaluated.bed)"
        echo -e "extended_scored_rows\\t\$(scored_rows ${prefix}.3primeExtended.peaks.evaluated.bed)"
        if [ -f ${prefix}.CallPasPerGene.sh ]; then
            echo -e "callpas_jobs\\t\$(wc -l < ${prefix}.CallPasPerGene.sh)"
        else
            echo -e "callpas_jobs\\tNA"
        fi
        echo -e "duplicate_original_gene_names\\t\$(awk 'END{print NR>0 ? NR-1 : 0}' ${prefix}.annotation.duplicate_gene_names.tsv)"
    } > ${prefix}.scapture.metrics.tsv

    # ---------------------------------------------------------------------
    # 7. Validate that at least one category produced usable BED rows.
    # ---------------------------------------------------------------------
    EX_VALID=\$(valid_rows ${prefix}.exonic.peaks.evaluated.bed)
    IN_VALID=\$(valid_rows ${prefix}.intronic.peaks.evaluated.bed)
    EXT_VALID=\$(valid_rows ${prefix}.3primeExtended.peaks.evaluated.bed)
    VALID_PEAKS=\$(( EX_VALID + IN_VALID + EXT_VALID ))
    EX_SCORED=\$(scored_rows ${prefix}.exonic.peaks.evaluated.bed)
    IN_SCORED=\$(scored_rows ${prefix}.intronic.peaks.evaluated.bed)
    EXT_SCORED=\$(scored_rows ${prefix}.3primeExtended.peaks.evaluated.bed)
    SCORED_PEAKS=\$(( EX_SCORED + IN_SCORED + EXT_SCORED ))

    if [ "\$VALID_PEAKS" -eq 0 ]; then
        echo "ERROR: SCAPTURE PAScall produced no valid evaluated peaks" >&2
        echo "Last PAScall log lines:" >&2
        tail -80 ${prefix}.PAScall.log >&2 || true
        exit 1
    fi

    if [ "\$SCORED_PEAKS" -eq 0 ]; then
        echo "ERROR: SCAPTURE DeepPASS produced no numeric prediction scores" >&2
        echo "This usually means FASTA extraction failed before DeepPASS prediction." >&2
        echo "Last PAScall log lines:" >&2
        tail -80 ${prefix}.PAScall.log >&2 || true
        for log in ${prefix}.*.DeepPASS.predict.log; do
            [ -f "\$log" ] || continue
            echo "Last DeepPASS log lines from \$log:" >&2
            tail -40 "\$log" >&2 || true
        done
        exit 1
    fi

    # ---------------------------------------------------------------------
    # 8. SCAPTURE PASmerge — run ONCE PER CATEGORY.
    #
    #    The --peak input is a 2-column TSV (sample_name, peak_file_path).
    #    Each run of PASmerge writes \${out_prefix}.Integrated.bed.
    #    We then concatenate the per-category Integrated.bed files into a
    #    single PAS reference for quantification, following the README's
    #    "Run scapture PASquant module" recipe.
    # ---------------------------------------------------------------------
    : > ${prefix}_merge.Integrated.bed
    MERGED_CATS=0

    for cat in exonic intronic 3primeExtended; do
        bed="${prefix}.\${cat}.peaks.evaluated.bed"
        nrows=\$(scored_rows "\$bed")
        if [ "\$nrows" -eq 0 ]; then
            echo "INFO: skipping PASmerge for \${cat} (no DeepPASS-scored rows in \$bed)" \\
                >> ${prefix}.PASmerge.log
            continue
        fi

        printf "${prefix}\\t\$bed\\n" > peaks_list.\${cat}.txt

        scapture -m PASmerge \\
            -o ${prefix}_merge.\${cat} \\
            --peak peaks_list.\${cat}.txt \\
            &>> ${prefix}.PASmerge.log

        if [ -s ${prefix}_merge.\${cat}.Integrated.bed ] \\
        && [ "\$(valid_rows ${prefix}_merge.\${cat}.Integrated.bed)" -gt 0 ]; then
            cat ${prefix}_merge.\${cat}.Integrated.bed >> ${prefix}_merge.Integrated.bed
            MERGED_CATS=\$(( MERGED_CATS + 1 ))
        else
            echo "WARN: PASmerge produced no valid Integrated.bed for \${cat}" \\
                >> ${prefix}.PASmerge.log
        fi
    done

    if [ "\$MERGED_CATS" -eq 0 ] || [ ! -s ${prefix}_merge.Integrated.bed ]; then
        echo "ERROR: SCAPTURE PASmerge produced no valid Integrated.bed across categories" >&2
        tail -80 ${prefix}.PASmerge.log >&2 || true
        exit 1
    fi

    # Sort the combined PAS reference for PASquant
    sort -k1,1 -k2,2n ${prefix}_merge.Integrated.bed -o ${prefix}_merge.Integrated.bed

    # Normalize to BED12: featureCounts rejects files where line column counts
    # differ (PASmerge exonic=15 cols, intronic=14 cols → mixed after concat).
    awk 'BEGIN{OFS="\\t"} NF>=12{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12}' \\
        ${prefix}_merge.Integrated.bed > ${prefix}_merge.Integrated.bed.tmp && \\
        mv ${prefix}_merge.Integrated.bed.tmp ${prefix}_merge.Integrated.bed

    # ---------------------------------------------------------------------
    # 9. SCAPTURE PASquant
    #
    #    Note: scapture forwards -p only on the --celltype branch, not the
    #    --celllist branch, so threading is governed by scapture_quant's
    #    internals here. Kept -p for forward compatibility.
    # ---------------------------------------------------------------------
    scapture -m PASquant \\
        -b ${bam} \\
        --celllist barcodes.txt \\
        --pas ${prefix}_merge.Integrated.bed \\
        -o ${prefix} \\
        -p ${ncores} \\
        &> ${prefix}.PASquant.log || true

    # ---------------------------------------------------------------------
    # 10. Build normalized site catalog
    # ---------------------------------------------------------------------
    run_scapture_filter.py \\
        --exonic    ${prefix}.exonic.peaks.evaluated.bed \\
        --intronic  ${prefix}.intronic.peaks.evaluated.bed \\
        --extended  ${prefix}.3primeExtended.peaks.evaluated.bed \\
        --threshold ${score_thr} \\
        --out-tsv   ${prefix}.scapture.site_catalog.tsv

    SITE_ROWS=\$(catalog_data_rows ${prefix}.scapture.site_catalog.tsv)
    if [ "\$SITE_ROWS" -eq 0 ]; then
        echo "ERROR: empty SCAPTURE site catalog after DeepPASS threshold ${score_thr}" >&2
        exit 1
    fi

    # ---------------------------------------------------------------------
    # 11. Manifest
    # ---------------------------------------------------------------------
    {
        printf "field\\tvalue\\n"
        printf "library_id\\t${prefix}\\n"
        printf "bam\\t${bam.name}\\n"
        printf "gtf\\t${gtf.name}\\n"
        printf "scapture_gtf\\t${prefix}.scapture.annotation.gtf\\n"
        printf "gene_name_map\\t${prefix}.scapture.gene_name_map.tsv\\n"
        printf "genome_fasta\\t${genome_fasta.name}\\n"
        printf "genome_fai\\t${genome_fai.name}\\n"
        printf "rescue_pascall_dir\\t${has_resume ? resume_pascall_dir.name : 'none'}\\n"
        printf "species\\t${species}\\n"
        printf "read_length\\t${read_len}\\n"
        printf "peak_width\\t${peak_width}\\n"
        printf "deeppass_threshold\\t${score_thr}\\n"
        printf "exonic_valid\\t\$EX_VALID\\n"
        printf "intronic_valid\\t\$IN_VALID\\n"
        printf "extended_valid\\t\$EXT_VALID\\n"
        printf "valid_peaks\\t\$VALID_PEAKS\\n"
        printf "exonic_scored\\t\$EX_SCORED\\n"
        printf "intronic_scored\\t\$IN_SCORED\\n"
        printf "extended_scored\\t\$EXT_SCORED\\n"
        printf "scored_peaks\\t\$SCORED_PEAKS\\n"
        printf "site_catalog_rows\\t\$SITE_ROWS\\n"
        printf "merged_categories\\t\$MERGED_CATS\\n"
        printf "duplicate_original_gene_names\\t\$(awk 'END{print NR>0 ? NR-1 : 0}' ${prefix}.annotation.duplicate_gene_names.tsv)\\n"
        printf "source\\tSCAPTURE\\n"
    } > ${prefix}.scapture.manifest.tsv
    """

    stub:
    def prefix = "${meta.library_id}"
    """
    printf "site_id\\tgene_id\\tchrom\\tstart\\tend\\tstrand\\tsite_class\\tsite_source\\tpriming_flag\\n" \\
        > ${prefix}.scapture.site_catalog.tsv
    printf "field\\tvalue\\nmode\\tstub\\n" > ${prefix}.scapture.manifest.tsv
    printf "metric\\tvalue\\nmode\\tstub\\n" > ${prefix}.scapture.metrics.tsv
    printf "gene_id\\tscapture_gene_name\\toriginal_gene_name\\n" > ${prefix}.scapture.gene_name_map.tsv
    printf "gene_count\\toriginal_gene_name\\n" > ${prefix}.annotation.duplicate_gene_names.tsv
    touch ${prefix}.annotation.log ${prefix}.PAScall.log ${prefix}.PASmerge.log ${prefix}.PASquant.log
    """
}
