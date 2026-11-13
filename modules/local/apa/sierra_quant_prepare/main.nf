process SIERRA_QUANT_PREPARE {
    tag "${meta.library_id}:${group_level}:${group_id}"
    label 'process_single'

    conda "${projectDir}/envs/r.yml"
    container "ghcr.io/martings/sierra:84b71c4"

    input:
    tuple val(meta), val(group_level), val(group_id), path(grouped_bam), path(grouped_bai), path(pas_reference), path(cell_annotations), path(gtf)

    output:
    tuple val(meta), val(group_level), val(group_id), path("${prefix}.whitelist.txt"), path("${prefix}.peak_sites.tsv"), path("${prefix}.site_id_map.tsv"), path(grouped_bam), path(grouped_bai), path(gtf), path(pas_reference), path(cell_annotations), emit: prepared
    path "${prefix}.prepare.versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def safeGroupId    = group_id.toString().replaceAll(/[^A-Za-z0-9._-]+/, '_')
    def safeLibraryId  = meta.library_id.toString().replaceAll(/[^A-Za-z0-9._-]+/, '_')
    def safeGroupLevel = group_level.toString().replaceAll(/[^A-Za-z0-9._-]+/, '_')
    prefix             = task.ext.prefix ?: "${safeLibraryId}.${safeGroupLevel}.${safeGroupId}"
    def b64LibraryId   = meta.library_id.toString().bytes.encodeBase64().toString()
    def b64GroupLevel  = group_level.toString().bytes.encodeBase64().toString()
    def b64GroupId     = group_id.toString().bytes.encodeBase64().toString()
    """
    test -s "${pas_reference}" && test -s "${cell_annotations}" || { echo 'Sierra reference/annotations missing or empty' >&2; exit 1; }

    Rscript - <<'RSCRIPT' "${cell_annotations}" "${pas_reference}" "${prefix}.whitelist.txt" "${prefix}.peak_sites.tsv" "${prefix}.site_id_map.tsv" "${b64LibraryId}" "${b64GroupLevel}" "${b64GroupId}"
    suppressPackageStartupMessages({
        library(data.table)
        library(jsonlite)
    })
    args <- commandArgs(trailingOnly = TRUE)
    ann_path  <- args[[1]]
    pas_path  <- args[[2]]
    wl_out    <- args[[3]]
    peak_out  <- args[[4]]
    map_out   <- args[[5]]
    lib_id    <- rawToChar(base64_dec(args[[6]]))
    grp_level <- rawToChar(base64_dec(args[[7]]))
    grp_id    <- rawToChar(base64_dec(args[[8]]))

    ann <- fread(ann_path)
    missing_cols <- setdiff(c("library_id", "barcode_corrected"), colnames(ann))
    if (length(missing_cols) > 0) {
        stop(paste("cell_annotations.tsv missing required columns:", paste(missing_cols, collapse = ", ")))
    }
    ann <- ann[ann[["library_id"]] == lib_id, ]
    if (grp_level == "cluster") {
        if (!("cluster_id" %in% colnames(ann))) stop("cell_annotations.tsv missing cluster_id")
        ann <- ann[cluster_id == grp_id]
    } else if (grp_level == "cell_type") {
        if (!("cell_type" %in% colnames(ann))) stop("cell_annotations.tsv missing cell_type")
        ann <- ann[cell_type == grp_id]
    } else if (grp_level == "cell") {
        if (!("cell_id" %in% colnames(ann))) stop("cell_annotations.tsv missing cell_id")
        ann <- ann[cell_id == grp_id]
    }
    whitelist_bcs <- unique(ann[["barcode_corrected"]])
    whitelist_bcs <- whitelist_bcs[!is.na(whitelist_bcs) & nzchar(whitelist_bcs)]
    writeLines(whitelist_bcs, wl_out)

    pas_ref <- fread(pas_path)
    s_start <- as.integer(pas_ref[["start"]])
    s_end   <- pmax(as.integer(pas_ref[["end"]]), s_start + 1L)
    s_str   <- ifelse(pas_ref[["strand"]] == "+", 1L, -1L)

    # First column must be row names matching Sierra::CountPeaks read.table(..., row.names=1)
    peak_dt <- data.table(
        PeakID    = pas_ref[["pas_reference_id"]],
        Chr       = pas_ref[["chrom"]],
        Start     = s_start,
        End       = as.integer(pas_ref[["end"]]),
        Strand    = s_str,
        Gene      = pas_ref[["gene_id"]],
        Fit.start = s_start,
        Fit.end   = s_end
    )
    con <- file(peak_out, "wt")
    writeLines(paste(c("Chr", "Start", "End", "Strand", "Gene", "Fit.start", "Fit.end"), collapse = "\\t"), con)
    close(con)
    fwrite(peak_dt, peak_out, sep = "\\t", append = TRUE, col.names = FALSE, quote = FALSE)

    sierra_names <- paste(pas_ref[["gene_id"]], pas_ref[["chrom"]], paste(s_start, s_end, sep = "-"), as.character(s_str), sep = ":")
    map_dt <- data.table(sierra_name = sierra_names, site_id = pas_ref[["site_id"]])
    fwrite(map_dt, map_out, sep = "\\t", quote = FALSE)
    RSCRIPT

    cat <<-END_VERSIONS > ${prefix}.prepare.versions.yml
    "${task.process}":
        r-data.table: \$(Rscript -e "cat(as.character(packageVersion('data.table')))")
    END_VERSIONS
    """

    stub:
    def safeGroupId    = group_id.toString().replaceAll(/[^A-Za-z0-9._-]+/, '_')
    def safeLibraryId  = meta.library_id.toString().replaceAll(/[^A-Za-z0-9._-]+/, '_')
    def safeGroupLevel = group_level.toString().replaceAll(/[^A-Za-z0-9._-]+/, '_')
    prefix             = task.ext.prefix ?: "${safeLibraryId}.${safeGroupLevel}.${safeGroupId}"
    """
    printf 'AAACCCAAGAAACACT-1\n' > "${prefix}.whitelist.txt"
    printf 'Chr\tStart\tEnd\tStrand\tGene\tFit.start\tFit.end\npas_1\tchr1\t100\t120\t1\tG1\t100\t120\n' > "${prefix}.peak_sites.tsv"
    printf 'sierra_name\tsite_id\nG1:chr1:100-120:1\tsite_1\n' > "${prefix}.site_id_map.tsv"
    printf '"%s":\n  r-data.table: "1.15.4"\n' "${task.process}" > "${prefix}.prepare.versions.yml"
    """
}
