#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
    library(Sierra)
    library(Matrix)
    library(Rsamtools)
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
    make_option("--ncores",            type = "integer",   default = 1L)
)
opt <- parse_args(OptionParser(option_list = option_list))

log_con <- file(opt$log, open = "wt")
log <- function(...) {
    msg <- paste0(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
    message(msg)
    writeLines(msg, log_con)
}

# ── Diagnostics: peek at BAM CB tags vs whitelist ─────────────────────────────
log("Reading cell annotations: ", opt$`cell-annotations`)
ann <- fread(opt$`cell-annotations`)
whitelist_barcodes <- ann$barcode_corrected
log("Whitelist barcodes: ", length(whitelist_barcodes), " (first: ", whitelist_barcodes[1], ")")

param_diag <- ScanBamParam(tag = "CB", what = character(0))
diag_aln   <- scanBam(opt$bam, param = param_diag)[[1]]
bam_cbs    <- unique(na.omit(diag_aln$tag$CB))
n_overlap   <- sum(bam_cbs %in% whitelist_barcodes)
log("BAM CB tags sampled: ", length(bam_cbs), " unique (first: ", bam_cbs[1], ")")
log("Barcode overlap with whitelist: ", n_overlap, " / ", length(bam_cbs))

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

# ── Run Sierra::CountPeaks ────────────────────────────────────────────────────
out_dir <- paste0(opt$output, "_sierra_mex")
dir.create(out_dir, showWarnings = FALSE)

log("Running Sierra::CountPeaks on: ", opt$bam)
result <- tryCatch(
    Sierra::CountPeaks(
        peak.sites.file = peak_file,
        gtf.file        = opt$gtf,
        bamfile         = opt$bam,
        whitelist.file  = whitelist_file,
        output.dir      = out_dir,
        countUMI        = TRUE,
        ncores          = opt$ncores
    ),
    error = function(e) {
        log("CountPeaks error: ", conditionMessage(e))
        NULL
    }
)

# ── Parse output or write empty TSV if no reads were found ───────────────────
empty_dt <- data.table(
    library_id   = character(),
    group_level  = character(),
    group_id     = character(),
    site_id      = character(),
    cell_barcode = character(),
    umi_count    = integer()
)

mat_file     <- file.path(out_dir, "matrix.mtx.gz")
barcode_file <- file.path(out_dir, "barcodes.tsv.gz")
feature_file <- file.path(out_dir, "sitenames.tsv.gz")

if (!file.exists(mat_file)) {
    log("WARNING: no matrix output from Sierra — writing empty TSV")
    fwrite(empty_dt, opt$output, sep = "\t")
    close(log_con)
    quit(status = 0)
}

log("CountPeaks complete, reading MEX output")
sm       <- readMM(mat_file)
barcodes <- fread(barcode_file, header = FALSE)$V1
features <- fread(feature_file, header = FALSE)$V1
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
