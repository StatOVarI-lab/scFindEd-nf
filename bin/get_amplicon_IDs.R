#!/usr/bin/env Rscript

library(Biostrings)
library(cli)
library(readr)
library(readxl)
library(dplyr)

infer_edit_blocks_df <- function(hdr_template, amplicon_seq, max_mismatch = 5) {
  hdr <- Biostrings::DNAString(hdr_template)
  amp <- Biostrings::DNAString(amplicon_seq)

  hits <- Biostrings::matchPattern(
    pattern = hdr,
    subject = amp,
    max.mismatch = max_mismatch,
    fixed = FALSE
  )
  orientation <- "forward"
  hdr_oriented <- hdr

  if (length(hits) == 0) {
    hdr_rc <- Biostrings::reverseComplement(hdr)
    hits <- Biostrings::matchPattern(
      pattern = hdr_rc,
      subject = amp,
      max.mismatch = max_mismatch,
      fixed = FALSE
    )
    orientation <- "reverse_complement"
    hdr_oriented <- hdr_rc
  }

  if (length(hits) == 0) {
    return(tibble(
      hdr_orientation = NA_character_,
      block_index = NA_integer_,
      start = NA_integer_,
      length = NA_integer_,
      ref = NA_character_,
      change = NA_character_
    ))
  }

  if (length(hits) > 1) {
    cli::cli_abort("HDR template matches multiple locations in amplicon.")
  }

  hit_start <- start(hits)[1]
  hit_end <- end(hits)[1]

  amp_subseq <- Biostrings::subseq(amp, start = hit_start, end = hit_end)

  amp_chars <- strsplit(as.character(amp_subseq), "")[[1]]
  hdr_chars <- strsplit(as.character(hdr_oriented), "")[[1]]

  mismatch_idx <- which(amp_chars != hdr_chars)

  if (length(mismatch_idx) == 0) {
    return(tibble(
      hdr_orientation = orientation,
      block_index = NA_integer_,
      start = NA_integer_,
      length = NA_integer_,
      ref = NA_character_,
      change = NA_character_
    ))
  }

  block_id <- cumsum(c(1, diff(mismatch_idx) != 1))
  blocks <- split(mismatch_idx, block_id)

  tibble(
  hdr_orientation = orientation,
  block_index = seq_along(blocks),
  start = vapply(blocks, function(idx) as.integer(hit_start + min(idx) - 1), integer(1)),
  length = vapply(blocks, function(idx) as.integer(length(idx)), integer(1)),
  ref = vapply(blocks, function(idx) paste0(amp_chars[idx], collapse = ""), character(1)),
  change = vapply(blocks, function(idx) paste0(hdr_chars[idx], collapse = ""), character(1))
)
}

get_amplicon_IDs_and_blocks <- function(
  SNP_info,
  design_summary,
  verbose = FALSE
) {
  if (verbose) {
    cli::cli_process_start(
      msg = "Assigning AmpIDs and mutation blocks...",
      msg_done = "AmpIDs and mutation blocks assigned."
    )
  }

  temp_design_summary <- design_summary[-nrow(design_summary), ]

  required_snp_cols <- c("HDR_template")
  missing_snp_cols <- setdiff(required_snp_cols, colnames(SNP_info))
  if (length(missing_snp_cols) > 0) {
    cli::cli_abort(c(
      "Missing required column(s) in SNP_info:",
      "x" = paste(missing_snp_cols, collapse = ", ")
    ))
  }

  required_design_cols <- c("amplicon_seq", "AmpID")
  missing_design_cols <- setdiff(required_design_cols, colnames(temp_design_summary))
  if (length(missing_design_cols) > 0) {
    cli::cli_abort(c(
      "Missing required column(s) in design_summary:",
      "x" = paste(missing_design_cols, collapse = ", ")
    ))
  }

  amp_sequences <- Biostrings::DNAStringSet(temp_design_summary$amplicon_seq)
  ref_AmpIDs <- temp_design_summary$AmpID

  out_list <- vector("list", length = nrow(SNP_info))

  for (i in seq_len(nrow(SNP_info))) {
    hdr_template <- SNP_info$HDR_template[i]
    hdr_dna <- Biostrings::DNAString(hdr_template)

    matches <- Biostrings::vcountPattern(
      hdr_dna,
      amp_sequences,
      max.mismatch = 5,
      fixed = FALSE
    )

    if (sum(matches > 0) == 0) {
      matches <- Biostrings::vcountPattern(
        Biostrings::reverseComplement(hdr_dna),
        amp_sequences,
        max.mismatch = 5,
        fixed = FALSE
      )
    }

    if (sum(matches > 0) > 1) {
      cli::cli_abort("HDR_template {i} matches multiple amplicon sequences.")
    }

    if (sum(matches > 0) == 0) {
      block_df <- tibble(
        hdr_orientation = NA_character_,
        block_index = NA_integer_,
        start = NA_integer_,
        length = NA_integer_,
        ref = NA_character_,
        change = NA_character_
      )
      AmpID <- NA_character_
    } else {
      match_idx <- which(matches > 0)[1]
      AmpID <- ref_AmpIDs[match_idx]

      block_df <- infer_edit_blocks_df(
        hdr_template = hdr_template,
        amplicon_seq = temp_design_summary$amplicon_seq[match_idx],
        max_mismatch = 5
      )
    }

    snp_row <- SNP_info[i, , drop = FALSE]
    snp_row$AmpID <- AmpID

    out_list[[i]] <- bind_cols(
      snp_row[rep(1, nrow(block_df)), , drop = FALSE],
      block_df
    )
  }

  out <- bind_rows(out_list)

  if (verbose) {
    cli::cli_process_done()
  }

  out
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

result <- get_amplicon_IDs_and_blocks(
  SNP_info = SNP_info,
  design_summary = design_summary,
  verbose = verbose
)

readr::write_csv(result, output_csv)