library(Biostrings)
library(cli)
library(readr)
library(readxl)


#' @keywords internal
get_amplicon_IDs <- function(
  SNP_info,
  design_summary,
  verbose
) {
  if (verbose) cli::cli_process_start(
    msg = "Assigning correct AmpID for each SNP in SNP_info...",
    msg_done = "Correct AmpID assigned for each SNP in SNP_info.")
  
  # Remove UTR barcode amplicon row from design summary file
  temp_design_summary <- design_summary[-nrow(design_summary), ]
  
  # Extract out vector of important information as vectors
  HDR_templates <- SNP_info$HDR_template
  amp_sequences <- Biostrings::DNAStringSet(temp_design_summary$amplicon_seq)
  ref_AmpIDs <- temp_design_summary$AmpID
  
  AmpID <- c()
  
  for (i in 1:length(HDR_templates)) {
    HDR_template <- Biostrings::DNAString(HDR_templates[i])
    
    # Check for matches with both HDR_template and its reverse complement
    matches <- Biostrings::vcountPattern(HDR_template, amp_sequences, 
                                         max.mismatch = 5, fixed = FALSE)
    # Abort process if multiple matches found
    if (sum(matches) > 1) {
      cli::cli_abort("HDR_template {i} matches multiple amplicon sequences.")
    }
    
    # If no matches found, try instead with the reverse complement
    if (sum(matches) == 0) {
      matches <- Biostrings::vcountPattern(Biostrings::reverseComplement(HDR_template), 
                                           amp_sequences, max.mismatch = 5, fixed = FALSE)
      # Abort process if multiple matches found
      if (sum(matches) > 1) {
        cli::cli_abort("HDR_template {i} matches multiple amplicon sequences.")
      }
    } 
    
    if (sum(matches) == 0) {
      # Assign no AmpID if no matches found for either strand
      AmpID <- c(AmpID, NA)
    } else if (sum(matches) == 1) {
      # Assign associated AmpID if exactly one match found
      match_idx <- which(matches > 0)[1]
      AmpID <- c(AmpID, ref_AmpIDs[match_idx])
    }
  }
  
  # Add AmpID as a column of SNP_info
  SNP_info <- cbind(SNP_info[,1:2], AmpID, SNP_info[,3:ncol(SNP_info)])
  
  if (verbose) cli::cli_process_done()
  
  return(SNP_info)
}


# Example use case
# Needs copy of 9824-design-summary-with-amp-seqs.csv
design_summary <- readr::read_csv("9824-design-summary-with-amp-seqs.csv",
                                  show_col_types = FALSE)
# Needs copy of "20260310_Pilot1_SNP_final.xlsx"
SNP_info <- readxl::read_excel("20260310_Pilot1_SNP_final.xlsx")

SNP_info_with_AmpIDs <- get_amplicon_IDs(SNP_info, design_summary, verbose = TRUE)


