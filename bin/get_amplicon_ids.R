#!/usr/bin/env Rscript

library(Biostrings)
library(readr)
library(readxl)
library(dplyr)


infer_edit_blocks_df <- function(
  hdr_template,
  amplicon_seq,
  max_mismatch = 5,
  snp_row = NA_integer_,
  amp_id = NA_character_,
  gene = NA_character_,
  snp_name = NA_character_
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


  # ---------------------------------------------------------
  # Try reverse complement if there is no forward match
  # ---------------------------------------------------------
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


  # ---------------------------------------------------------
  # No matching region found
  # ---------------------------------------------------------
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


  # ---------------------------------------------------------
  # More than one match within the same amplicon is ambiguous
  # ---------------------------------------------------------
  if (length(hits) > 1) {

    stop(
      "HDR template matches multiple locations in amplicon.",
      call. = FALSE
    )
  }


  # ---------------------------------------------------------
  # Coordinates of HDR match within the amplicon
  # ---------------------------------------------------------
  hit_start <- start(hits)[1]
  hit_end <- end(hits)[1]

  amp_length <- length(amp)
  hdr_length <- length(hdr_oriented)


  # ---------------------------------------------------------
  # Handle HDR templates that extend beyond the amplicon
  # boundary.
  #
  # Example:
  #
  #   amplicon length = 232 bp
  #   HDR length      = 100 bp
  #   HDR match       = 135-234
  #
  # The amplicon ends at 232, so the final 2 bp of the HDR
  # cannot be compared.
  #
  # In this case, use only the part of the HDR that overlaps
  # the amplicon and issue a visible warning.
  # ---------------------------------------------------------

  left_bases_outside <- max(
    0L,
    1L - hit_start
  )

  right_bases_outside <- max(
    0L,
    hit_end - amp_length
  )


  if (
    left_bases_outside > 0 ||
    right_bases_outside > 0
  ) {

    original_hit_start <- hit_start
    original_hit_end <- hit_end

    overlap_start <- max(
      1L,
      hit_start
    )

    overlap_end <- min(
      amp_length,
      hit_end
    )

    overlap_length <-
      overlap_end -
      overlap_start +
      1L


    # -------------------------------------------------------
    # Safety check: there must be some overlap
    # -------------------------------------------------------
    if (overlap_length <= 0) {

      stop(
        paste0(
          "HDR template has no valid overlap with amplicon ",
          amp_id,
          "."
        ),
        call. = FALSE
      )
    }


    # -------------------------------------------------------
    # Work out which part of the oriented HDR corresponds
    # to the sequence that is actually present in the
    # amplicon.
    # -------------------------------------------------------
    hdr_overlap_start <-
      1L +
      left_bases_outside

    hdr_overlap_end <-
      hdr_length -
      right_bases_outside


    if (
      hdr_overlap_start > hdr_overlap_end ||
      hdr_overlap_start < 1 ||
      hdr_overlap_end > hdr_length
    ) {

      stop(
        paste0(
          "Unable to determine valid HDR overlap for ",
          "SNP row ", snp_row,
          ", AmpID ", amp_id,
          "."
        ),
        call. = FALSE
      )
    }


    warning_text <- paste0(
      "\n",
      "============================================================\n",
      "WARNING: HDR TEMPLATE EXTENDS BEYOND AMPLICON BOUNDARY\n",
      "------------------------------------------------------------\n",
      "SNP row: ", snp_row, "\n",
      "Gene: ", gene, "\n",
      "SNP: ", snp_name, "\n",
      "AmpID: ", amp_id, "\n",
      "HDR orientation: ", orientation, "\n",
      "HDR length: ", hdr_length, " bp\n",
      "Amplicon length: ", amp_length, " bp\n",
      "Original HDR match: ",
      original_hit_start, "-", original_hit_end, "\n",
      "Bases outside left boundary: ",
      left_bases_outside, " bp\n",
      "Bases outside right boundary: ",
      right_bases_outside, " bp\n",
      "Overlap used for mutation inference: ",
      overlap_start, "-", overlap_end,
      " (", overlap_length, " bp)\n",
      "The part of the HDR outside the amplicon cannot be ",
      "evaluated and will be ignored.\n",
      "Please check whether the amplicon boundary is correct.\n",
      "============================================================\n"
    )


    # -------------------------------------------------------
    # Print to stdout.
    #
    # With `debug true` in the Nextflow process this will
    # also be visible in the main LSF .out file.
    # -------------------------------------------------------
    cat(
      warning_text,
      file = stdout()
    )


    # -------------------------------------------------------
    # Also retain the warning in task stderr (.command.err)
    # -------------------------------------------------------
    warning(
      warning_text,
      call. = FALSE
    )


    # -------------------------------------------------------
    # Restrict both sequences to their overlapping region
    # -------------------------------------------------------
    hit_start <- overlap_start
    hit_end <- overlap_end

    hdr_oriented <- Biostrings::subseq(
      hdr_oriented,
      start = hdr_overlap_start,
      end = hdr_overlap_end
    )
  }


  # ---------------------------------------------------------
  # Extract the matching region from the amplicon
  # ---------------------------------------------------------
  amp_subseq <- Biostrings::subseq(
    amp,
    start = hit_start,
    end = hit_end
  )


  # ---------------------------------------------------------
  # Safety check:
  # the amplicon and HDR sections must now have the same
  # length before base-by-base comparison
  # ---------------------------------------------------------
  if (
    length(amp_subseq) !=
    length(hdr_oriented)
  ) {

    stop(
      paste0(
        "HDR/amplicon overlap length mismatch for ",
        "SNP row ", snp_row,
        ", AmpID ", amp_id,
        ": amplicon overlap = ",
        length(amp_subseq),
        " bp; HDR overlap = ",
        length(hdr_oriented),
        " bp."
      ),
      call. = FALSE
    )
  }


  amp_chars <- strsplit(
    as.character(amp_subseq),
    ""
  )[[1]]

  hdr_chars <- strsplit(
    as.character(hdr_oriented),
    ""
  )[[1]]


  # ---------------------------------------------------------
  # Identify bases that differ between the amplicon and HDR
  # ---------------------------------------------------------
  mismatch_idx <- which(
    amp_chars != hdr_chars
  )


  # ---------------------------------------------------------
  # HDR matches the available amplicon sequence exactly
  # ---------------------------------------------------------
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


  # ---------------------------------------------------------
  # Group consecutive mismatches into edit blocks
  # ---------------------------------------------------------
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

    block_index = seq_along(
      blocks
    ),

    start = vapply(
      blocks,
      function(idx) {

        as.integer(
          hit_start +
          min(idx) -
          1
        )
      },
      integer(1)
    ),

    length = vapply(
      blocks,
      function(idx) {

        as.integer(
          length(idx)
        )
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


  if (
    verbose &&
    length(identity_amplicons) > 0
  ) {

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

    hdr_template <-
      SNP_info$HDR_template[i]

    hdr_dna <- Biostrings::DNAString(
      hdr_template
    )


    matches <- Biostrings::vcountPattern(
      hdr_dna,
      amp_sequences,
      max.mismatch = 5,
      fixed = FALSE
    )


    # -------------------------------------------------------
    # Try reverse complement if no forward match
    # -------------------------------------------------------
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


    # -------------------------------------------------------
    # More than one amplicon matched
    # -------------------------------------------------------
    if (sum(matches > 0) > 1) {

      stop(
        sprintf(
          "HDR_template %s matches multiple amplicon sequences.",
          i
        ),
        call. = FALSE
      )
    }


    # -------------------------------------------------------
    # No matching amplicon
    # -------------------------------------------------------
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


      # -----------------------------------------------------
      # Infer mutation blocks.
      #
      # Pass identifying information so that any boundary
      # warning clearly identifies the affected SNP/AmpID.
      # -----------------------------------------------------
      block_df <- infer_edit_blocks_df(
        hdr_template = hdr_template,

        amplicon_seq =
          temp_design_summary$amplicon_seq[
            match_idx
          ],

        max_mismatch = 5,

        snp_row = i,

        amp_id = AmpID,

        gene =
          if ("Gene" %in% colnames(SNP_info)) {
            SNP_info$Gene[i]
          } else {
            NA_character_
          },

        snp_name =
          if ("SNP" %in% colnames(SNP_info)) {
            SNP_info$SNP[i]
          } else {
            NA_character_
          }
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
# TAMPL194941
#
# TAMPL194941,TAMPLxxxxxx
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