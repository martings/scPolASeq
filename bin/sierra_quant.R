#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
    library(Sierra)
    library(Matrix)
    library(Rsamtools)
    library(foreach)
})

option_list <- list(
    make_option("--bam",               type = "character", help = "Grouped BAM file"),
    make_option("--gtf",               type = "character", help = "GTF annotation file"),
    make_option("--pas-reference",     type = "character", help = "pas_reference.tsv"),
    make_option("--cell-annotations",  type = "character", help = "cell_annotations.tsv"),
    make_option("--library-id",        type = "character", help = "Library ID"),
    make_option("--group-level",       type = "character", help = "Grouping level"),
    make_option("--group-id",          type = "character", help = "Group ID"),
    make_option("--output",            type = "character", help = "Output TSV path"),
    make_option("--log",               type = "character", help = "Log file path"),
    make_option("--versions",          type = "character", help = "YAML runtime version output"),
    make_option("--ncores",            type = "integer",   default = 1L),
    make_option("--strandedness",      type = "character", default = "forward", help = "Library strandedness: forward, reverse, or unstranded (default: forward)"),
    make_option("--strand-mode",       type = "integer",   default = NULL,      help = "Bioconductor strandMode integer: 1 (forward), 2 (reverse), 0 (unstranded)")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!is.null(opt$versions)) {
    writeLines(c("SIERRA_QUANT:",
                 paste0("  Sierra: '", as.character(packageVersion("Sierra")), "'"),
                 paste0("  R: '", getRversion(), "'")), opt$versions)
}

log_con <- file(opt$log, open = "wt")
log <- function(...) {
    msg <- paste0(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
    message(msg)
    writeLines(msg, log_con)
}

# ── Bail out cleanly on empty / invalid BAM ───────────────────────────────────
empty_dt <- data.table(
    library_id   = character(),
    group_level  = character(),
    group_id     = character(),
    site_id      = character(),
    cell_barcode = character(),
    umi_count    = integer()
)
if (file.info(opt$bam)$size == 0) {
    log("WARNING: BAM file is empty — writing empty TSV and skipping")
    fwrite(empty_dt, opt$output, sep = "\t")
    close(log_con)
    quit(status = 0)
}

required_ann_cols <- c("library_id", "barcode_corrected")
subset_annotations <- function(ann, library_id, group_level, group_id) {
    missing_cols <- setdiff(required_ann_cols, colnames(ann))
    if (length(missing_cols) > 0) {
        stop(paste("cell_annotations.tsv missing required columns:", paste(missing_cols, collapse = ", ")))
    }

    library_mask <- ann[["library_id"]] == library_id
    ann <- ann[library_mask, ]

    if (group_level == "cluster") {
        if (!("cluster_id" %in% colnames(ann))) {
            stop("cell_annotations.tsv missing cluster_id needed for cluster-level Sierra quant")
        }
        ann <- ann[cluster_id == group_id]
    } else if (group_level == "cell_type") {
        if (!("cell_type" %in% colnames(ann))) {
            stop("cell_annotations.tsv missing cell_type needed for cell_type-level Sierra quant")
        }
        ann <- ann[cell_type == group_id]
    } else if (group_level == "cell") {
        if (!("cell_id" %in% colnames(ann))) {
            stop("cell_annotations.tsv missing cell_id needed for cell-level Sierra quant")
        }
        ann <- ann[cell_id == group_id]
    }

    ann
}

# ── Diagnostics: peek at BAM CB tags vs whitelist ─────────────────────────────
log("Reading cell annotations: ", opt$`cell-annotations`)
ann <- fread(opt$`cell-annotations`)
ann <- subset_annotations(ann, opt$`library-id`, opt$`group-level`, opt$`group-id`)
whitelist_barcodes <- unique(ann$barcode_corrected)
log("Scoped annotations rows: ", nrow(ann))
log("Whitelist barcodes after scoping: ", length(whitelist_barcodes),
    if (length(whitelist_barcodes) > 0) paste0(" (first: ", whitelist_barcodes[1], ")") else "")
if (length(whitelist_barcodes) == 0) {
    log("WARNING: scoped whitelist is empty — writing empty TSV")
    empty_dt <- data.table(
        library_id   = character(),
        group_level  = character(),
        group_id     = character(),
        site_id      = character(),
        cell_barcode = character(),
        umi_count    = integer()
    )
    fwrite(empty_dt, opt$output, sep = "\t")
    close(log_con)
    quit(status = 0)
}

param_diag <- ScanBamParam(tag = "CB", what = character(0))
diag_aln   <- tryCatch(
    scanBam(opt$bam, param = param_diag)[[1]],
    error = function(e) {
        log("ERROR: BAM diagnostic scan failed: ", conditionMessage(e))
        close(log_con)
        stop(e)
    }
)
bam_cbs    <- unique(na.omit(diag_aln$tag$CB))
n_overlap   <- sum(bam_cbs %in% whitelist_barcodes)
log("BAM CB tags sampled: ", length(bam_cbs), " unique",
    if (length(bam_cbs) > 0) paste0(" (first: ", bam_cbs[1], ")") else "")
log("Barcode overlap with whitelist: ", n_overlap, " / ", length(bam_cbs))
if (n_overlap == 0) {
    log("No BAM barcodes overlap the scoped whitelist — writing an empty TSV")
    fwrite(empty_dt, opt$output, sep = "\t")
    close(log_con)
    quit(status = 0)
}

whitelist_file <- tempfile(fileext = ".txt")
writeLines(whitelist_barcodes, whitelist_file)

# ── Build peak sites file ──────────────────────────────────────────────────────
log("Reading PAS reference: ", opt$`pas-reference`)
pas_ref <- fread(opt$`pas-reference`)

# Sierra internally filters Fit.start < Fit.end — ensure at least 1 bp window
peak_df <- data.frame(
    Chr       = pas_ref$chrom,
    Start     = as.integer(pas_ref$start),
    End       = as.integer(pas_ref$end),
    Strand    = ifelse(pas_ref$strand == "+", 1L, -1L),
    Gene      = pas_ref$gene_id,
    Fit.start = as.integer(pas_ref$start),
    Fit.end   = pmax(as.integer(pas_ref$end), as.integer(pas_ref$start) + 1L),
    stringsAsFactors = FALSE,
    row.names = pas_ref$pas_reference_id
)
log("Peak sites: ", nrow(peak_df), " (Fit.start < Fit.end for ",
    sum(peak_df$Fit.start < peak_df$Fit.end), ")")

# Build Sierra-name → canonical site_id lookup.
# Sierra constructs peak names as gene_id:chrom:start-sierra_end:strand_int
# (chr-prefix chrom, integer strand, end=pmax(end,start+1)).
# This map translates back to the site catalog's site_id format.
sierra_end_vec   <- pmax(as.integer(pas_ref$end), as.integer(pas_ref$start) + 1L)
sierra_strand_vec <- ifelse(pas_ref$strand == "+", "1", "-1")
sierra_names     <- paste(pas_ref$gene_id, pas_ref$chrom,
                          paste(as.integer(pas_ref$start), sierra_end_vec, sep="-"),
                          sierra_strand_vec, sep=":")
sierra_to_site_id <- setNames(pas_ref$site_id, sierra_names)

peak_file <- tempfile(fileext = ".tsv")
write.table(peak_df, peak_file, sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

# ── Resolve strandedness / strandMode ─────────────────────────────────────────
strand_mode <- if (!is.null(opt$`strand-mode`)) {
    as.integer(opt$`strand-mode`)
} else {
    switch(tolower(opt$strandedness),
        "reverse"    = 2L,
        "unstranded" = 0L,
        1L # default forward
    )
}
log("Strandedness configured: ", opt$strandedness, " (Bioconductor strandMode: ", strand_mode, ")")

# ── Extensible Sierra::CountPeaks implementation with strandMode support ──────
CountPeaks_with_strandMode <- function (peak.sites.file, gtf.file, bamfile, whitelist.file,
    output.dir, countUMI = TRUE, ncores = 1, chr.names = NULL,
    filter.chr = FALSE, gene.symbol.ref = "gene_name", CBtag = "CB",
    UMItag = "UB", strandMode = 1L)
{
    lock <- tempfile()
    whitelist.bc <- read.table(whitelist.file, stringsAsFactors = FALSE)
    whitelist.bc <- whitelist.bc[, 1]
    n.bcs <- length(whitelist.bc)
    message("There are ", n.bcs, " whitelist barcodes.")
    n.columns <- n.bcs + 1
    genes.ref <- Sierra:::make_reference(gtf_file = gtf.file, chr.names = chr.names,
        filter.chr = filter.chr, gene.symbol.ref = gene.symbol.ref)
    chr.names <- as.character(unique(genes.ref$chr))
    n.genes <- nrow(genes.ref)
    peak.sites <- read.table(peak.sites.file, header = T, sep = "\t",
        quote = "", stringsAsFactors = FALSE)
    n.total.sites <- nrow(peak.sites)
    message("There are ", n.total.sites, " sites")
    message("Doing counting for each site (strandMode = ", strandMode, ")...")
    system.name <- Sys.info()["sysname"]
    new_cl <- FALSE
    if (system.name == "Windows") {
        new_cl <- TRUE
        cluster <- parallel::makePSOCKcluster(rep("localhost", ncores))
        doParallel::registerDoParallel(cluster)
    } else {
        doParallel::registerDoParallel(cores = ncores)
    }
    mat.to.write <- foreach::foreach(each.chr = chr.names, .combine = "rbind",
        .packages = c("magrittr")) %dopar% {
        mat.per.chr <- c()
        message("Processing chr: ", each.chr)
        for (strand in c(1, -1)) {
            message(" and strand ", strand)
            if (strandMode == 2L) {
                # Reverse stranded (10x v1, 10x 5' SE): cDNA read is sense
                isMinusStrand <- if (strand == 1) TRUE else FALSE
                flag_param <- Rsamtools::scanBamFlag(isMinusStrand = isMinusStrand)
            } else if (strandMode == 0L) {
                # Unstranded: accept reads from either strand
                flag_param <- Rsamtools::scanBamFlag()
            } else {
                # Forward stranded (default): cDNA read is antisense
                isMinusStrand <- if (strand == 1) FALSE else TRUE
                flag_param <- Rsamtools::scanBamFlag(isMinusStrand = isMinusStrand)
            }
            peak.sites.chr <- dplyr::filter(peak.sites, Chr ==
                each.chr & Strand == strand) %>% dplyr::select(Gene,
                Chr, Fit.start, Fit.end, Strand)
            peak.sites.chr$Fit.start <- as.integer(peak.sites.chr$Fit.start)
            peak.sites.chr$Fit.end <- as.integer(peak.sites.chr$Fit.end)
            peak.sites.chr <- dplyr::filter(peak.sites.chr, Fit.start <
                Fit.end)
            if (nrow(peak.sites.chr) == 0) {
                next
            }
            which <- GenomicRanges::GRanges(seqnames = each.chr,
                ranges = IRanges::IRanges(1, max(peak.sites.chr$Fit.end)))
            param <- Rsamtools::ScanBamParam(tag = c(CBtag, UMItag),
                which = which, flag = flag_param)
            aln <- GenomicAlignments::readGAlignments(bamfile,
                param = param)
            nobarcodes <- which(unlist(is.na(GenomicRanges::mcols(aln)[CBtag])))
            noUMI <- which(unlist(is.na(GenomicRanges::mcols(aln)[UMItag])))
            to.remove <- dplyr::union(nobarcodes, noUMI)
            if (length(to.remove) > 0) {
                aln <- aln[-to.remove]
            }
            whitelist.pos <- which(unlist(GenomicRanges::mcols(aln)[CBtag]) %in%
                whitelist.bc)
            aln <- aln[whitelist.pos]
            if (length(aln) == 0) {
                next
            }
            if (countUMI) {
                GenomicRanges::mcols(aln)$CB_UB <- paste0(unlist(GenomicRanges::mcols(aln)[CBtag]),
                  "_", unlist(GenomicRanges::mcols(aln)[UMItag]))
                uniqUMIs <- which(!duplicated(GenomicRanges::mcols(aln)$CB_UB))
                aln <- aln[uniqUMIs]
            }
            aln <- GenomicRanges::split(aln, unlist(GenomicRanges::mcols(aln)[CBtag]))
            polyA.GR <- GenomicRanges::GRanges(seqnames = peak.sites.chr$Chr,
                IRanges::IRanges(start = peak.sites.chr$Fit.start,
                  end = as.integer(peak.sites.chr$Fit.end)))
            n.polyA <- length(polyA.GR)
            barcodes.gene <- names(aln)
            res <- sapply(barcodes.gene, function(x) GenomicRanges::countOverlaps(polyA.GR,
                aln[[x]]))
            res.mat <- matrix(0L, nrow = n.polyA, ncol = n.bcs)
            res.mat[, match(barcodes.gene, whitelist.bc)] <- res
            mat.per.strand <- Matrix::Matrix(res.mat, sparse = TRUE)
            polyA.ids <- paste0(peak.sites.chr$Gene, ":", peak.sites.chr$Chr,
                ":", peak.sites.chr$Fit.start, "-", peak.sites.chr$Fit.end,
                ":", peak.sites.chr$Strand)
            rownames(mat.per.strand) <- polyA.ids
            if (is.null(mat.per.chr)) {
                mat.per.chr <- mat.per.strand
            } else {
                mat.per.chr <- rbind(mat.per.chr, mat.per.strand)
            }
        }
        return(mat.per.chr)
    }
    if (new_cl) {
        parallel::stopCluster(cluster)
    }
    if (!dir.exists(output.dir)) {
        dir.create(output.dir)
    }
    Matrix::writeMM(mat.to.write, file = paste0(output.dir, "/matrix.mtx"))
    write.table(whitelist.bc, file = paste0(output.dir, "/barcodes.tsv"),
        quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(rownames(mat.to.write), file = paste0(output.dir,
        "/sitenames.tsv"), quote = FALSE, row.names = FALSE,
        col.names = FALSE)
    R.utils::gzip(paste0(output.dir, "/matrix.mtx"), overwrite = TRUE)
    R.utils::gzip(paste0(output.dir, "/barcodes.tsv"), overwrite = TRUE)
    R.utils::gzip(paste0(output.dir, "/sitenames.tsv"), overwrite = TRUE)
}

# ── Run Sierra::CountPeaks ────────────────────────────────────────────────────
out_dir <- paste0(opt$output, "_sierra_mex")
dir.create(out_dir, showWarnings = FALSE)

log("Running Sierra::CountPeaks (strandMode = ", strand_mode, ") on: ", opt$bam)
tryCatch(
    CountPeaks_with_strandMode(
        peak.sites.file = peak_file,
        gtf.file        = opt$gtf,
        bamfile         = opt$bam,
        whitelist.file  = whitelist_file,
        output.dir      = out_dir,
        countUMI        = TRUE,
        ncores          = opt$ncores,
        strandMode      = strand_mode
    ),
    error = function(e) {
        log("ERROR: CountPeaks failed: ", conditionMessage(e))
        close(log_con)
        stop(e)
    }
)

# ── Parse Sierra output; a missing matrix after valid input is a task failure ──
mat_file     <- file.path(out_dir, "matrix.mtx.gz")
barcode_file <- file.path(out_dir, "barcodes.tsv.gz")
feature_file <- file.path(out_dir, "sitenames.tsv.gz")

if (!file.exists(mat_file)) {
    log("ERROR: Sierra returned without a matrix for a BAM with matching barcodes")
    close(log_con)
    stop("Sierra matrix.mtx.gz is missing")
}

log("CountPeaks complete, reading MEX output")
sm       <- readMM(mat_file)
barcodes <- fread(barcode_file, sep = "\n", header = FALSE)$V1
features <- fread(feature_file, sep = "\n", header = FALSE)$V1
rownames(sm) <- features
colnames(sm) <- barcodes
log("Matrix loaded: ", nrow(sm), " peaks x ", ncol(sm), " cells")

idx <- which(sm != 0, arr.ind = TRUE)
if (length(idx) == 0) {
    log("WARNING: count matrix is all zeros — writing empty TSV")
    fwrite(empty_dt, opt$output, sep = "\t")
} else {
    raw_pas_ids   <- rownames(sm)[idx[, 1]]
    canonical_ids <- sierra_to_site_id[raw_pas_ids]
    n_total       <- length(unique(raw_pas_ids))
    n_unmatched   <- sum(is.na(canonical_ids))
    log("site_id translation: ", n_total - n_unmatched, "/", n_total, " unique peaks matched to canonical site_id")
    if (n_unmatched > 0)
        log("WARNING: ", n_unmatched, " peaks had no site_id match — Sierra name format may have changed; kept as-is")
    canonical_ids[is.na(canonical_ids)] <- raw_pas_ids[is.na(canonical_ids)]

    dt <- data.table(
        library_id   = opt$`library-id`,
        group_level  = opt$`group-level`,
        group_id     = opt$`group-id`,
        site_id      = canonical_ids,
        cell_barcode = colnames(sm)[idx[, 2]],
        umi_count    = sm[idx]
    )
    log("Non-zero peak-cell pairs: ", nrow(dt))
    fwrite(dt, opt$output, sep = "\t")
}

log("Output written: ", opt$output)
close(log_con)
