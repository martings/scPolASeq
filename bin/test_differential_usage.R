#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  result <- list()
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    value <- args[[i + 1]]
    result[[key]] <- value
    i <- i + 2
  }
  result
}

parsed <- parse_args(args)

usage <- read.delim(parsed[["apa-usage"]], sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

if (nrow(usage) == 0) {
  stats <- data.frame(
    contrast_id = character(),
    gene_id = character(),
    site_id = character(),
    group_level = character(),
    logFC = numeric(),
    delta_usage = numeric(),
    pvalue = numeric(),
    fdr = numeric(),
    test_class = character(),
    stringsAsFactors = FALSE
  )
  write.table(stats, file = parsed[["out-stats"]], sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(metric = "n_tests", value = 0), file = parsed[["out-summary"]], sep = "\t", quote = FALSE, row.names = FALSE)
  quit(save = "no", status = 0)
}

usage$usage_fraction <- as.numeric(usage$usage_fraction)
usage$group_level <- as.character(usage$group_level)
usage$group_id <- as.character(usage$group_id)

split_group <- function(group_id) {
  tokens <- strsplit(group_id, "__", fixed = TRUE)[[1]]
  c(
    sample_id = ifelse(length(tokens) >= 1, tokens[[1]], NA_character_),
    condition = ifelse(length(tokens) >= 2, tokens[[2]], NA_character_),
    cell_type = ifelse(length(tokens) >= 3, tokens[[3]], NA_character_)
  )
}

rows <- list()
for (gene_id in unique(usage$gene_id)) {
  gene_usage <- usage[usage$gene_id == gene_id, , drop = FALSE]
  for (site_id in unique(gene_usage$site_id)) {
    site_usage <- gene_usage[gene_usage$site_id == site_id, , drop = FALSE]
    sample_rows <- site_usage[site_usage$group_level == "sample_condition_celltype", , drop = FALSE]
    if (nrow(sample_rows) > 0) {
      group_parts <- t(vapply(sample_rows$group_id, split_group, FUN.VALUE = c(sample_id = "", condition = "", cell_type = "")))
      sample_rows$sample_id <- group_parts[, "sample_id"]
      sample_rows$condition <- group_parts[, "condition"]
      sample_rows$cell_type <- group_parts[, "cell_type"]

      for (cell_type in unique(sample_rows$cell_type)) {
        cell_rows <- sample_rows[sample_rows$cell_type == cell_type, , drop = FALSE]
        conditions <- unique(cell_rows$condition)
        if (length(conditions) < 2) {
          next
        }
        cond_a <- conditions[[1]]
        cond_b <- conditions[[2]]
        values_a <- cell_rows$usage_fraction[cell_rows$condition == cond_a]
        values_b <- cell_rows$usage_fraction[cell_rows$condition == cond_b]
        if (length(values_a) == 0 || length(values_b) == 0) {
          next
        }
        delta_usage <- mean(values_b) - mean(values_a)
        logFC <- log2((mean(values_b) + 1e-6) / (mean(values_a) + 1e-6))
        if (length(values_a) >= 2 && length(values_b) >= 2) {
          test <- t.test(values_b, values_a)
          pvalue <- test$p.value
          test_class <- "replicate_aware_t_test"
        } else {
          pvalue <- NA_real_
          test_class <- "descriptive_only"
        }
        rows[[length(rows) + 1]] <- data.frame(
          contrast_id = paste(cond_b, "vs", cond_a, cell_type, sep = "__"),
          gene_id = gene_id,
          site_id = site_id,
          group_level = "sample_condition_celltype",
          logFC = logFC,
          delta_usage = delta_usage,
          pvalue = pvalue,
          fdr = NA_real_,
          test_class = test_class,
          stringsAsFactors = FALSE
        )
      }
    }

    for (level in setdiff(unique(site_usage$group_level), "sample_condition_celltype")) {
      level_rows <- site_usage[site_usage$group_level == level, , drop = FALSE]
      rows[[length(rows) + 1]] <- data.frame(
        contrast_id = paste("descriptive", level, sep = "__"),
        gene_id = gene_id,
        site_id = site_id,
        group_level = level,
        logFC = 0,
        delta_usage = ifelse(nrow(level_rows) > 0, max(level_rows$usage_fraction) - min(level_rows$usage_fraction), 0),
        pvalue = NA_real_,
        fdr = NA_real_,
        test_class = "descriptive_only",
        stringsAsFactors = FALSE
      )
    }
  }
}

stats <- if (length(rows) > 0) do.call(rbind, rows) else data.frame(
  contrast_id = character(),
  gene_id = character(),
  site_id = character(),
  group_level = character(),
  logFC = numeric(),
  delta_usage = numeric(),
  pvalue = numeric(),
  fdr = numeric(),
  test_class = character(),
  stringsAsFactors = FALSE
)

if (nrow(stats) > 0) {
  stats$fdr <- p.adjust(stats$pvalue, method = "BH")
}

write.table(stats, file = parsed[["out-stats"]], sep = "\t", quote = FALSE, row.names = FALSE)

summary_df <- data.frame(
  metric = c("n_usage_rows", "n_stats_rows"),
  value = c(nrow(usage), nrow(stats)),
  stringsAsFactors = FALSE
)
write.table(summary_df, file = parsed[["out-summary"]], sep = "\t", quote = FALSE, row.names = FALSE)
