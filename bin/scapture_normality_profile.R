# Loaded by Rscript only during SCAPTURE PAScall. Seed Normality.R independently
# for each transcript, so GNU Parallel scheduling cannot change its samples.
argv <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", argv, value = TRUE)
if (length(script_arg) == 1L &&
    identical(basename(sub("^--file=", "", script_arg)), "Normality.R")) {
    base_text <- Sys.getenv("SCPOLASEQ_SCAPTURE_NORMALITY_SEED", unset = "1")
    if (!grepl("^[0-9]+$", base_text)) {
        stop("SCPOLASEQ_SCAPTURE_NORMALITY_SEED must be an integer from 0 to 2147483646")
    }
    base_seed <- suppressWarnings(as.numeric(base_text))
    if (!is.finite(base_seed) || base_seed >= 2147483647) {
        stop("SCPOLASEQ_SCAPTURE_NORMALITY_SEED must be an integer from 0 to 2147483646")
    }
    normality_args <- commandArgs(trailingOnly = TRUE)
    if (length(normality_args) < 1L) {
        stop("Normality.R requires the peak file as its first argument")
    }
    seed <- base_seed
    for (byte in utf8ToInt(basename(normality_args[[1L]]))) {
        seed <- (seed * 131 + byte) %% 2147483647
    }
    set.seed(as.integer(seed))
}
