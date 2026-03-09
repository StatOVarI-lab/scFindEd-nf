#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: reference_from_designSummary.R <design_summary_csv>", call. = FALSE)
}

amp_data <- read.csv(
  args[1],
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fill = TRUE,
  na.strings = character(0)   # keep literal "NA" as text
)

if (!("AmpID" %in% names(amp_data))) {
  stop("Input must contain column: AmpID", call. = FALSE)
}

has_amplicon_seq <- "amplicon_seq" %in% names(amp_data)
has_amplicon_sequence <- "Amplicon sequence" %in% names(amp_data)

if (!has_amplicon_seq && !has_amplicon_sequence) {
  stop("Input must contain 'amplicon_seq' and/or 'Amplicon sequence' column", call. = FALSE)
}

seq1 <- if (has_amplicon_seq) amp_data[["amplicon_seq"]] else rep("", nrow(amp_data))
seq2 <- if (has_amplicon_sequence) amp_data[["Amplicon sequence"]] else rep("", nrow(amp_data))

# Use amplicon_seq if it is not empty and not the literal string "NA" (case-insensitive),
# otherwise fall back to Amplicon sequence.
bad_seq1 <- is.na(seq1) | !nzchar(seq1) | toupper(seq1) == "NA"
seqs <- seq1
seqs[bad_seq1] <- seq2[bad_seq1]

# Drop rows still missing usable sequence
keep <- !is.na(seqs) & nzchar(seqs) & toupper(seqs) != "NA"
amp_data <- amp_data[keep, , drop = FALSE]
seqs <- seqs[keep]

output_odd  <- paste0(">", amp_data$AmpID)
output_even <- seqs

output <- rep("", length(output_odd) * 2)
output[seq(1, length(output) - 1, by = 2)] <- output_odd
output[seq(2, length(output),     by = 2)] <- output_even

writeLines(output, con = "ref_amplicons.fa")