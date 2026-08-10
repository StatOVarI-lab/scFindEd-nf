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

# Read only required columns from the barcode/read-name table
cell_readnames <- fread(
    readnames_file,
    select = c("read_name", "cell_barcode")
)

# Read UTR alignment data
amplicon_reads <- fread(
    amplicon_file
)

# Key both tables for an efficient join
setkey(cell_readnames, read_name)
setkey(amplicon_reads, read_name)

# Add cell barcode information.
# Keep all UTR reads, even if no barcode is matched.
utr_reads <- cell_readnames[
    amplicon_reads,
    on = "read_name"
]

# Add sample/amplicon identifiers by reference
utr_reads[, sample_id := sample_id]
utr_reads[, amp_id := amp_id]

# Put identifiers first
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