#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
    stop(
        "Usage: report_utr_amplicon.R ",
        "<sample_id> <amp_id> <amplicon_tsv> ",
        "<readnames_tsv> <output_tsv>"
    )
}

sample_id     <- args[[1]]
amp_id        <- args[[2]]
amplicon_file <- args[[3]]
readnames_file <- args[[4]]
output_file   <- args[[5]]

amplicon_reads <- read_tsv(
    amplicon_file,
    show_col_types = FALSE
)

cell_readnames <- read_tsv(
    readnames_file,
    show_col_types = FALSE
)

utr_reads <- amplicon_reads %>%
    left_join(
        cell_readnames,
        by = "read_name"
    ) %>%
    mutate(
        sample_id = sample_id,
        amp_id = amp_id
    ) %>%
    select(
        sample_id,
        amp_id,
        everything()
    )

write_tsv(
    utr_reads,
    output_file
)