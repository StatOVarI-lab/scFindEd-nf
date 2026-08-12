#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
    stop(
        "Usage: report_utr_amplicon.R ",
        "<sample_id> <amp_id> <amplicon_tsv> ",
        "<readnames_tsv> <output_tsv>",
        call. = FALSE
    )
}

sample_id      <- args[[1]]
amp_id         <- args[[2]]
amplicon_file  <- args[[3]]
readnames_file <- args[[4]]
output_file    <- args[[5]]

# Read barcode/read-name table.
# This file has no header, so assign column names explicitly.
cell_readnames <- fread(
    readnames_file,
    header = FALSE,
    col.names = c("read_name", "cell_barcode")
)

# Read UTR alignment data.
amplicon_reads <- fread(
    amplicon_file
)

# Check expected join column exists.
if (!"read_name" %in% names(amplicon_reads)) {
    stop(
        "Column 'read_name' not found in amplicon input.",
        call. = FALSE
    )
}

# Key tables for efficient matching.
setkey(cell_readnames, read_name)
setkey(amplicon_reads, read_name)

# Keep all UTR reads and append the matching cell barcode.
utr_reads <- cell_readnames[
    amplicon_reads,
    on = "read_name"
]

# Add sample and amplicon identifiers by reference.
utr_reads[, sample_id := sample_id]
utr_reads[, amp_id := amp_id]

# Put identifiers first.
setcolorder(
    utr_reads,
    c(
        "sample_id",
        "amp_id",
        setdiff(names(utr_reads), c("sample_id", "amp_id"))
    )
)

fwrite(
    utr_reads,
    output_file,
    sep = "\t",
    quote = FALSE
)