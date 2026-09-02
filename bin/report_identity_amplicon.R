#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
    stop(
        "Usage: report_identity_amplicon.R ",
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

# Read the full alignment table for the requested identity amplicon.
amplicon_reads <- fread(
    amplicon_file
)

# Basic sanity check.
if (!"read_name" %in% names(amplicon_reads)) {
    stop(
        "Column 'read_name' not found in identity amplicon TSV.",
        call. = FALSE
    )
}

# Read-name / cell-barcode table has no header.
cell_readnames <- fread(
    readnames_file,
    header = FALSE,
    col.names = c(
        "read_name",
        "cell_barcode"
    )
)

# Remove reads without a valid cell barcode.
cell_readnames <- cell_readnames[
    !is.na(cell_barcode) &
    cell_barcode != ""
]

# Keep only identity-amplicon alignments associated with a cell barcode.
identity_reads <- merge(
    amplicon_reads,
    cell_readnames,
    by = "read_name",
    all = FALSE,
    sort = FALSE
)

# Add sample and amplicon identifiers.
identity_reads[, sample_id := sample_id]
identity_reads[, amp_id := amp_id]

# Put the identifying columns first.
first_cols <- c(
    "sample_id",
    "cell_barcode",
    "amp_id",
    "read_name"
)

other_cols <- setdiff(
    names(identity_reads),
    first_cols
)

setcolorder(
    identity_reads,
    c(
        first_cols,
        other_cols
    )
)

fwrite(
    identity_reads,
    output_file,
    sep = "\t",
    quote = FALSE
)