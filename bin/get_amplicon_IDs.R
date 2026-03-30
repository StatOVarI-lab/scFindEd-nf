#!/usr/bin/env Rscript

library(Biostrings)
library(cli)
library(readr)
library(readxl)

get_amplicon_IDs <- function(
  SNP_info,
  design_summary,
  verbose = FALSE
) {
  if (verbose) {
    cli::cli_process_start(
      msg = "Assigning correct AmpID for each SNP in SNP_info...",
      msg_done = "Correct AmpID assigned for each SNP in SNP_info."
    )
  }

  # Remove UTR barcode amplicon row from design summary file
  temp_design_summary <- design_summary[-nrow(design_summary), ]

  # Basic input checks
  required_snp_cols <- c("HDR_template")
  missing_snp_cols <- setdiff(required_snp_cols, colnames(SNP_info))
  if (length(missing_snp_cols) > 0) {
    cli::cli_abort(c(
      "Missing required column(s) in SNP_info:",
      "x" = paste(missing_snp_cols, collapse = ", ")
    ))
  }

  required_design_cols <- c("amplicon_seq", "AmpID")
  missing_design_cols <- setdiff(required_design_cols, colnames(design_summary))
  if (length(missing_design_cols) > 0) {
    cli::cli_abort(c(
      "Missing required column(s) in design_summary:",
      "x" = paste(missing_design_cols, collapse = ", ")
    ))
  }

  # Extract important information as vectors
  HDR_templates <- SNP_info$HDR_template
  amp_sequences <- Biostrings::DNAStringSet(temp_design_summary$amplicon_seq)
  ref_AmpIDs <- temp_design_summary$AmpID

  AmpID <- c()

  for (i in seq_along(HDR_templates)) {
    HDR_template <- Biostrings::DNAString(HDR_templates[i])

    # Check for matches with HDR_template
    matches <- Biostrings::vcountPattern(
      HDR_template,
      amp_sequences,
      max.mismatch = 5,
      fixed = FALSE
    )

    # Abort process if multiple matches found
    if (sum(matches) > 1) {
      cli::cli_abort("HDR_template {i} matches multiple amplicon sequences.")
    }

    # If no matches found, try instead with the reverse complement
    if (sum(matches) == 0) {
      matches <- Biostrings::vcountPattern(
        Biostrings::reverseComplement(HDR_template),
        amp_sequences,
        max.mismatch = 5,
        fixed = FALSE
      )

      if (sum(matches) > 1) {
        cli::cli_abort("HDR_template {i} matches multiple amplicon sequences.")
      }
    }

    if (sum(matches) == 0) {
      AmpID <- c(AmpID, NA)
    } else if (sum(matches) == 1) {
      match_idx <- which(matches > 0)[1]
      AmpID <- c(AmpID, ref_AmpIDs[match_idx])
    }
  }

  # Add AmpID as a column of SNP_info
  SNP_info <- cbind(SNP_info[, 1:2], AmpID, SNP_info[, 3:ncol(SNP_info)])

  if (verbose) {
    cli::cli_process_done()
  }

  return(SNP_info)
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    "Usage: get_amplicon_ids.R <snp_info_xlsx> <design_summary_csv> <output_csv> [--verbose]",
    call. = FALSE
  )
}

snp_info_file <- args[1]
design_summary_file <- args[2]
output_csv <- args[3]
verbose <- "--verbose" %in% args

design_summary <- readr::read_csv(design_summary_file, show_col_types = FALSE)
SNP_info <- readxl::read_excel(snp_info_file)

SNP_info_with_AmpIDs <- get_amplicon_IDs(
  SNP_info = SNP_info,
  design_summary = design_summary,
  verbose = verbose
)

readr::write_csv(SNP_info_with_AmpIDs, output_csv)