#!/usr/bin/env Rscript

library(Biostrings)
library(readr)
library(readxl)
library(dplyr)


infer_edit_blocks_df <- function(
  hdr_template,
  amplicon_seq,
  max_mismatch = 5
) {

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

    return(
      tibble(
        hdr_orientation = NA_character_,
        block_index = NA_integer_,
        start = NA_integer_,
        length = NA_integer_,
        ref = NA_character_,
        change = NA_character_
      )
    )
  }

  if (length(hits) > 1) {
    stop(
      "HDR template matches multiple locations in amplicon.",
      call. = FALSE
    )
  }

  hit_start <- start(hits)[1]
  hit_end <- end(hits)[1]

  amp_subseq <- Biostrings::subseq(
    amp,
    start = hit_start,
    end = hit_end
  )

  amp_chars <- strsplit(
    as.character(amp_subseq),
    ""
  )[[1]]

  hdr_chars <- strsplit(
    as.character(hdr_oriented),
    ""
  )[[1]]

  mismatch_idx <- which(
    amp_chars != hdr_chars
  )

  if (length(mismatch_idx) == 0) {

    return(
      tibble(
        hdr_orientation = orientation,
        block_index = NA_integer_,
        start = NA_integer_,
        length = NA_integer_,
        ref = NA_character_,
        change = NA_character_
      )
    )
  }

  block_id <- cumsum(
    c(
      1,
      diff(mismatch_idx) != 1
    )
  )

  blocks <- split(
    mismatch_idx,
    block_id
  )

  tibble(
    hdr_orientation = orientation,

    block_index = seq_along(blocks),

    start = vapply(
      blocks,
      function(idx) {
        as.integer(
          hit_start + min(idx) - 1
        )
      },
      integer(1)
    ),

    length = vapply(
      blocks,
      function(idx) {
        as.integer(length(idx))
      },
      integer(1)
    ),

    ref = vapply(
      blocks,
      function(idx) {
        paste0(
          amp_chars[idx],
          collapse = ""
        )
      },
      character(1)
    ),

    change = vapply(
      blocks,
      function(idx) {
        paste0(
          hdr_chars[idx],
          collapse = ""
        )
      },
      character(1)
    )
  )
}


get_amplicon_IDs_and_blocks <- function(
  SNP_info,
  design_summary,
  identity_amplicons = character(0),
  verbose = FALSE
) {

  if (verbose) {
    message(
      "Assigning AmpIDs and mutation blocks..."
    )
  }


  # ---------------------------------------------------------
  # Validate AmpID column
  # ---------------------------------------------------------
  if (!"AmpID" %in% colnames(design_summary)) {
    stop(
      "Missing required column in design_summary: AmpID",
      call. = FALSE
    )
  }


  # ---------------------------------------------------------
  # Check that requested identity amplicons exist
  # ---------------------------------------------------------
  if (length(identity_amplicons) > 0) {

    missing_identity_amplicons <- setdiff(
      identity_amplicons,
      design_summary$AmpID
    )

    if (length(missing_identity_amplicons) > 0) {

      stop(
        "Identity amplicon(s) not found in design_summary: ",
        paste(
          missing_identity_amplicons,
          collapse = ", "
        ),
        call. = FALSE
      )
    }
  }


  # ---------------------------------------------------------
  # Exclude identity amplicons BY AmpID.
  #
  # Their position in the design CSV does not matter.
  # ---------------------------------------------------------
  temp_design_summary <- design_summary[
    !design_summary$AmpID %in% identity_amplicons,
    ,
    drop = FALSE
  ]


  if (verbose && length(identity_amplicons) > 0) {

    message(
      "Excluding identity amplicon(s) from HDR matching: ",
      paste(
        identity_amplicons,
        collapse = ", "
      )
    )
  }


  # ---------------------------------------------------------
  # Required SNP columns
  # ---------------------------------------------------------
  required_snp_cols <- c(
    "HDR_template"
  )

  missing_snp_cols <- setdiff(
    required_snp_cols,
    colnames(SNP_info)
  )

  if (length(missing_snp_cols) > 0) {

    stop(
      paste(
        "Missing required column(s) in SNP_info:",
        paste(
          missing_snp_cols,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }


  # ---------------------------------------------------------
  # Required design columns
  # ---------------------------------------------------------
  required_design_cols <- c(
    "amplicon_seq",
    "AmpID"
  )

  missing_design_cols <- setdiff(
    required_design_cols,
    colnames(temp_design_summary)
  )

  if (length(missing_design_cols) > 0) {

    stop(
      paste(
        "Missing required column(s) in design_summary:",
        paste(
          missing_design_cols,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }


  if (nrow(temp_design_summary) == 0) {

    stop(
      "No non-identity amplicons remain in design_summary.",
      call. = FALSE
    )
  }


  # ---------------------------------------------------------
  # Build searchable amplicon sequence set
  # ---------------------------------------------------------
  amp_sequences <- Biostrings::DNAStringSet(
    temp_design_summary$amplicon_seq
  )

  ref_AmpIDs <- temp_design_summary$AmpID


  # ---------------------------------------------------------
  # Match every HDR template against the available
  # non-identity amplicons
  # ---------------------------------------------------------
  out_list <- vector(
    "list",
    length = nrow(SNP_info)
  )


  for (i in seq_len(nrow(SNP_info))) {

    hdr_template <- SNP_info$HDR_template[i]

    hdr_dna <- Biostrings::DNAString(
      hdr_template
    )


    matches <- Biostrings::vcountPattern(
      hdr_dna,
      amp_sequences,
      max.mismatch = 5,
      fixed = FALSE
    )


    # Try reverse complement if no forward match
    if (sum(matches > 0) == 0) {

      matches <- Biostrings::vcountPattern(
        Biostrings::reverseComplement(
          hdr_dna
        ),
        amp_sequences,
        max.mismatch = 5,
        fixed = FALSE
      )
    }


    # More than one amplicon matched
    if (sum(matches > 0) > 1) {

      stop(
        sprintf(
          "HDR_template %s matches multiple amplicon sequences.",
          i
        ),
        call. = FALSE
      )
    }


    # No matching amplicon
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

      match_idx <- which(
        matches > 0
      )[1]

      AmpID <- ref_AmpIDs[
        match_idx
      ]


      block_df <- infer_edit_blocks_df(
        hdr_template = hdr_template,
        amplicon_seq =
          temp_design_summary$amplicon_seq[
            match_idx
          ],
        max_mismatch = 5
      )
    }


    snp_row <- SNP_info[
      i,
      ,
      drop = FALSE
    ]

    snp_row$AmpID <- AmpID


    out_list[[i]] <- bind_cols(
      snp_row[
        rep(
          1,
          nrow(block_df)
        ),
        ,
        drop = FALSE
      ],
      block_df
    )
  }


  out <- bind_rows(
    out_list
  )


  if (verbose) {
    message(
      "AmpIDs and mutation blocks assigned."
    )
  }


  out
}


# ============================================================
# Command-line arguments
# ============================================================

args <- commandArgs(
  trailingOnly = TRUE
)


if (length(args) < 3) {

  stop(
    paste(
      "Usage:",
      "get_amplicon_ids.R",
      "<snp_info_csv>",
      "<design_summary_csv>",
      "<output_csv>",
      "[identity_amplicon1,identity_amplicon2,...]",
      "[--verbose]"
    ),
    call. = FALSE
  )
}


snp_info_file <- args[1]
design_summary_file <- args[2]
output_csv <- args[3]


# ---------------------------------------------------------
# Optional comma-separated identity amplicons
#
# Examples:
#
# AMPL902654
#
# AMPL902654,AMPL902655
# ---------------------------------------------------------
identity_amplicon_arg <- ""

if (
  length(args) >= 4 &&
  args[4] != "--verbose"
) {

  identity_amplicon_arg <- args[4]
}


if (nzchar(identity_amplicon_arg)) {

  identity_amplicons <- strsplit(
    identity_amplicon_arg,
    ",",
    fixed = TRUE
  )[[1]]

  identity_amplicons <- trimws(
    identity_amplicons
  )

  identity_amplicons <- identity_amplicons[
    nzchar(identity_amplicons)
  ]

} else {

  identity_amplicons <- character(0)
}


verbose <- "--verbose" %in% args


# ============================================================
# Read inputs
# ============================================================

design_summary <- readr::read_csv(
  design_summary_file,
  show_col_types = FALSE
)

SNP_info <- readr::read_csv(
  snp_info_file,
  show_col_types = FALSE
)


# ============================================================
# Run
# ============================================================

result <- get_amplicon_IDs_and_blocks(
  SNP_info = SNP_info,
  design_summary = design_summary,
  identity_amplicons = identity_amplicons,
  verbose = verbose
)


# ============================================================
# Write output
# ============================================================

readr::write_csv(
  result,
  output_csv
)