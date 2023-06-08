#!/usr/bin/env Rscript

library(dplyr)
library(purrr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: get_counts.R <feature_counts_files>", call. = FALSE)
}

files <- args[1]

map(files, ~ {
  feature_counts <- read.table(.x, header = TRUE) %>%
    select(1, length(feature_counts))
}) %>%
  reduce(inner_join, by = "Geneid") -> count_table

write.table(count_table, file = "count_table.txt", row.names = FALSE, quote = FALSE)

colnames(count_table) <- sapply(strsplit(colnames(count_table), split = "_"), "[[", 1)
rownames(count_table) <- count_table$Geneid
count_table$Geneid <- NULL

count_table <- as.matrix(count_table)

saveRDS(count_table, file = "count_table.rds")
