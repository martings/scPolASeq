process SIERRA_COUNTPEAKS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/martings/sierra:84b71c4"

    input:
    tuple val(meta), path(bam), path(bai), path(whitelist)
    path peak_sites
    path gtf

    output:
    tuple val(meta), path("${prefix}/matrix.mtx.gz"), path("${prefix}/barcodes.tsv.gz"), path("${prefix}/sitenames.tsv.gz"), emit: counts
    tuple val(meta), path("${prefix}.sierra_countpeaks.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    def ncores = task.cpus ?: 1
    """
    test -e "${bam}" && test -e "${bai}" || { echo 'Sierra BAM/BAI missing' >&2; exit 1; }
    test -s "${peak_sites}" && test -s "${gtf}" || { echo 'Sierra peak_sites/GTF missing or empty' >&2; exit 1; }
    mkdir -p "${prefix}"

    Rscript - <<'RSCRIPT' ${args}
    suppressPackageStartupMessages({
        library(Sierra)
        library(Matrix)
        library(Rsamtools)
        library(R.utils)
    })

    cli <- commandArgs(trailingOnly = TRUE)
    get_opt <- function(flag, default = NULL) {
        idx <- which(cli == flag)
        if (length(idx) == 0) return(default)
        pos <- tail(idx, 1)
        if (pos >= length(cli)) stop("Missing value for ", flag)
        cli[[pos + 1]]
    }

    strand_mode <- get_opt("--strand-mode", NULL)
    strandedness <- get_opt("--strandedness", "forward")
    if (is.null(strand_mode)) {
        strand_mode <- switch(tolower(strandedness),
            "reverse"    = 2L,
            "unstranded" = 0L,
            "forward"    = 1L,
            stop("Unsupported strandedness: ", strandedness)
        )
    } else {
        strand_mode <- as.integer(strand_mode)
    }
    count_umi       <- as.logical(get_opt("--count-umi", "TRUE"))
    filter_chr      <- as.logical(get_opt("--filter-chr", "FALSE"))
    cb_tag          <- get_opt("--cb-tag", "CB")
    umi_tag         <- get_opt("--umi-tag", "UB")
    gene_symbol_ref <- get_opt("--gene-symbol-ref", "gene_name")

    bam_file       <- "${bam}"
    whitelist_file <- "${whitelist}"
    peak_file      <- "${peak_sites}"
    gtf_file       <- "${gtf}"
    out_dir        <- "${prefix}"
    log_file       <- "${prefix}.sierra_countpeaks.log"
    ncores         <- as.integer(${ncores})

    log_con <- file(log_file, open = "wt")
    log_msg <- function(...) {
        msg <- paste0(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
        message(msg)
        writeLines(msg, log_con)
    }

    write_empty_mex <- function(dir_path, barcodes_vec) {
        dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)
        empty_mat <- Matrix::Matrix(0, nrow = 0, ncol = length(barcodes_vec), sparse = TRUE)
        Matrix::writeMM(empty_mat, file = file.path(dir_path, "matrix.mtx"))
        R.utils::gzip(file.path(dir_path, "matrix.mtx"), overwrite = TRUE)
        writeLines(barcodes_vec, con = file.path(dir_path, "barcodes.tsv"))
        R.utils::gzip(file.path(dir_path, "barcodes.tsv"), overwrite = TRUE)
        writeLines(character(0), con = file.path(dir_path, "sitenames.tsv"))
        R.utils::gzip(file.path(dir_path, "sitenames.tsv"), overwrite = TRUE)
    }

    bam_size <- file.info(bam_file)[["size"]]
    if (is.na(bam_size) || bam_size == 0) {
        log_msg("WARNING: BAM file is empty — writing empty MEX output")
        write_empty_mex(out_dir, character(0))
        close(log_con)
        quit(status = 0)
    }

    # Always validate BAM readability first so corrupt BAMs fail closed even with empty whitelists.
    diag_aln <- tryCatch(
        Rsamtools::scanBam(bam_file, param = Rsamtools::ScanBamParam(tag = cb_tag, what = character(0)))[[1]],
        error = function(e) {
            log_msg("ERROR: BAM diagnostic scan failed: ", conditionMessage(e))
            close(log_con)
            stop(e)
        }
    )

    whitelist_bcs <- if (file.exists(whitelist_file) && file.info(whitelist_file)[["size"]] > 0) {
        unique(trimws(readLines(whitelist_file, warn = FALSE)))
    } else {
        character(0)
    }
    whitelist_bcs <- whitelist_bcs[nzchar(whitelist_bcs)]
    log_msg("Whitelist barcodes: ", length(whitelist_bcs))

    if (length(whitelist_bcs) == 0) {
        log_msg("WARNING: whitelist is empty — writing empty MEX output")
        write_empty_mex(out_dir, character(0))
        close(log_con)
        quit(status = 0)
    }

    bam_cbs   <- unique(na.omit(diag_aln[["tag"]][[cb_tag]]))
    n_overlap <- sum(bam_cbs %in% whitelist_bcs)
    log_msg("BAM ", cb_tag, " tags sampled: ", length(bam_cbs), "; overlap: ", n_overlap)
    if (n_overlap == 0) {
        log_msg("WARNING: zero BAM/whitelist barcode overlap — writing empty MEX output")
        write_empty_mex(out_dir, whitelist_bcs)
        close(log_con)
        quit(status = 0)
    }

    log_msg("Running Sierra::CountPeaks (strandMode = ", strand_mode, ", ncores = ", ncores, ")")
    tryCatch(
        Sierra::CountPeaks(
            peak.sites.file = peak_file,
            gtf.file        = gtf_file,
            bamfile         = bam_file,
            whitelist.file  = whitelist_file,
            output.dir      = out_dir,
            countUMI        = count_umi,
            ncores          = ncores,
            filter.chr      = filter_chr,
            gene.symbol.ref = gene_symbol_ref,
            CBtag           = cb_tag,
            UMItag          = umi_tag,
            strandMode      = strand_mode
        ),
        error = function(e) {
            log_msg("ERROR: Sierra::CountPeaks failed: ", conditionMessage(e))
            close(log_con)
            stop(e)
        }
    )

    for (fname in c("matrix.mtx.gz", "barcodes.tsv.gz", "sitenames.tsv.gz")) {
        fpath <- file.path(out_dir, fname)
        if (!file.exists(fpath) || file.info(fpath)[["size"]] == 0) {
            log_msg("ERROR: Missing MEX output: ", fpath)
            close(log_con)
            stop("Sierra::CountPeaks did not produce ", fname)
        }
    }
    log_msg("Sierra::CountPeaks complete.")
    close(log_con)
    RSCRIPT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sierra: \$(Rscript -e "cat(as.character(packageVersion('Sierra')))")
        r-base: \$(Rscript -e "cat(as.character(getRversion()))")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}"
    printf '%%%%MatrixMarket matrix coordinate integer general\n%%\n0 0 0\n' | gzip -c > "${prefix}/matrix.mtx.gz"
    printf 'AAACCCAAGAAACACT-1\n' | gzip -c > "${prefix}/barcodes.tsv.gz"
    : | gzip -c > "${prefix}/sitenames.tsv.gz"
    printf 'stub Sierra::CountPeaks\n' > "${prefix}.sierra_countpeaks.log"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sierra: 0.99.27.9001
        r-base: 4.3.3
    END_VERSIONS
    """
}
