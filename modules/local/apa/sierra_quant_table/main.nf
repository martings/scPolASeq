process SIERRA_QUANT_TABLE {
    tag "${meta.library_id}:${meta.group_level}:${meta.group_id}"
    label 'process_single'

    conda "${projectDir}/envs/r.yml"
    container "ghcr.io/martings/sierra:84b71c4"

    input:
    tuple val(orig_meta), val(group_level), val(group_id), val(meta), path(matrix_mtx), path(barcodes_tsv), path(sitenames_tsv), path(countpeaks_log), path(site_id_map), path(grouped_bam), path(grouped_bai), path(gtf), path(pas_reference), path(cell_annotations)

    output:
    tuple val(orig_meta), val(group_level), val(group_id), path("${prefix}.sierra_quant.tsv"), emit: quantified_pas
    path "${prefix}.sierra_quant.manifest.tsv", emit: manifest
    path "${prefix}.sierra_quant.log", emit: log
    path "${prefix}.versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix            = task.ext.prefix ?: "${meta.id}"
    def strand        = (meta.strand ?: 'forward').toString().toLowerCase()
    def b64LibraryId  = meta.library_id.toString().bytes.encodeBase64().toString()
    def b64GroupLevel = meta.group_level.toString().bytes.encodeBase64().toString()
    def b64GroupId    = meta.group_id.toString().bytes.encodeBase64().toString()
    def b64Stranded   = strand.bytes.encodeBase64().toString()
    """
    test -s "${matrix_mtx}" && test -s "${barcodes_tsv}" && test -e "${sitenames_tsv}" && test -s "${site_id_map}" || { echo 'Sierra MEX inputs missing' >&2; exit 1; }

    Rscript - <<'RSCRIPT' "${matrix_mtx}" "${barcodes_tsv}" "${sitenames_tsv}" "${site_id_map}" "${prefix}.sierra_quant.tsv" "${b64LibraryId}" "${b64GroupLevel}" "${b64GroupId}"
    suppressPackageStartupMessages({
        library(data.table)
        library(Matrix)
        library(jsonlite)
    })
    args      <- commandArgs(trailingOnly = TRUE)
    mat_file  <- args[[1]]
    bc_file   <- args[[2]]
    feat_file <- args[[3]]
    map_file  <- args[[4]]
    out_file  <- args[[5]]
    lib_id    <- rawToChar(base64_dec(args[[6]]))
    grp_level <- rawToChar(base64_dec(args[[7]]))
    grp_id    <- rawToChar(base64_dec(args[[8]]))

    empty_dt <- data.table(
        library_id   = character(),
        group_level  = character(),
        group_id     = character(),
        site_id      = character(),
        cell_barcode = character(),
        umi_count    = integer()
    )

    sm       <- Matrix::readMM(mat_file)
    barcodes <- if (file.info(bc_file)[["size"]] > 0) fread(bc_file, sep = "\\n", header = FALSE)[["V1"]] else character(0)
    features <- if (file.info(feat_file)[["size"]] > 0) fread(feat_file, sep = "\\n", header = FALSE)[["V1"]] else character(0)

    if (nrow(sm) == 0 || ncol(sm) == 0 || length(features) == 0 || length(barcodes) == 0) {
        fwrite(empty_dt, out_file, sep = "\\t")
        quit(status = 0)
    }

    rownames(sm) <- features
    colnames(sm) <- barcodes
    idx <- which(sm != 0, arr.ind = TRUE)
    if (length(idx) == 0) {
        fwrite(empty_dt, out_file, sep = "\\t")
        quit(status = 0)
    }

    map_dt <- fread(map_file)
    sierra_to_site_id <- setNames(map_dt[["site_id"]], map_dt[["sierra_name"]])
    raw_pas_ids       <- rownames(sm)[idx[, 1]]
    canonical_ids     <- sierra_to_site_id[raw_pas_ids]
    canonical_ids[is.na(canonical_ids)] <- raw_pas_ids[is.na(canonical_ids)]

    dt <- data.table(
        library_id   = lib_id,
        group_level  = grp_level,
        group_id     = grp_id,
        site_id      = canonical_ids,
        cell_barcode = colnames(sm)[idx[, 2]],
        umi_count    = as.integer(sm[idx])
    )
    fwrite(dt, out_file, sep = "\\t")
    RSCRIPT

    cp "${countpeaks_log}" "${prefix}.sierra_quant.log"

    Rscript - <<'RSCRIPT' "${prefix}.sierra_quant.manifest.tsv" "${grouped_bam.name}" "${grouped_bai.name}" "${gtf.name}" "${pas_reference.name}" "${cell_annotations.name}" "${b64Stranded}" "${b64GroupId}"
    suppressPackageStartupMessages(library(jsonlite))
    args <- commandArgs(trailingOnly = TRUE)
    stranded <- rawToChar(base64_dec(args[[7]]))
    grp_id   <- rawToChar(base64_dec(args[[8]]))
    lines <- c(
        "field\\tvalue",
        paste0("grouped_bam\\t", args[[2]]),
        paste0("grouped_bai\\t", args[[3]]),
        paste0("gtf\\t", args[[4]]),
        paste0("pas_reference\\t", args[[5]]),
        paste0("cell_annotations\\t", args[[6]]),
        paste0("strandedness\\t", stranded),
        "quant_source\\tSierra::CountPeaks",
        paste0("group_id\\t", grp_id)
    )
    writeLines(lines, args[[1]])
    RSCRIPT

    cat <<-END_VERSIONS > "${prefix}.versions.yml"
    SIERRA_QUANT:
        Sierra: '\$(Rscript -e "cat(as.character(packageVersion('Sierra')))")'
        R: '\$(Rscript -e "cat(as.character(getRversion()))")'
    END_VERSIONS

    test -s "${prefix}.sierra_quant.tsv" && test -s "${prefix}.versions.yml" || { echo 'Sierra table outputs missing or empty' >&2; exit 1; }
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf "library_id\\tgroup_level\\tgroup_id\\tsite_id\\tcell_barcode\\tumi_count\\n" > "${prefix}.sierra_quant.tsv"
    printf "field\\tvalue\\nquant_source\\tSierra::CountPeaks\\nmode\\tstub\\n" > "${prefix}.sierra_quant.manifest.tsv"
    printf "stub Sierra::CountPeaks\\n" > "${prefix}.sierra_quant.log"
    printf 'SIERRA_QUANT:\n  Sierra: "0.99.27.9001"\n  R: "4.3.3"\n' > "${prefix}.versions.yml"
    """
}
