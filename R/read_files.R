#' Read Samplesheet
#'
#' Read the Samplesheet file that is provided by users.
#' @details This function is to read the samplesheet file that is provided by users.
#' Spreadsheet has 4 columns with headers: Sample_ID, Condition, Rep, Fraction.
#' Sample_IDs need to be unique, Fractions should be Total, Supernatant, or Pellet.
#' @param file The input Samplesheet file name
#' @return A single data frame
#' @keywords samplesheet
#' @importFrom readr read_tsv
read_samplesheet <- function(file)
{
  load_samplesheet <- read_tsv(file, comment = "#", show_col_types = FALSE)
  # Check for correct column headers
  required_columns <- c("Sample_ID", "Condition", "Rep", "Fraction")
  if (!all(required_columns %in% names(load_samplesheet))) {
    stop("The input file must contain the following columns: Sample_ID, Condition, Rep, Fraction")
  }

  # Check for unique Sample_IDs
  if (any(duplicated(load_samplesheet$Sample_ID))) {
    stop("Sample_IDs need to be unique")
  }

  # Check for appropriately labeled Fractions
  allowed_fractions <- c("Total", "Supernatant", "Pellet")
  if (any(!load_samplesheet$Fraction %in% allowed_fractions)) {
    stop("Fractions should be one of: Total, Supernatant, or Pellet")
  }
  return(load_samplesheet)
}

#' Read count files
#'
#' Read a matrix of counts for each sample
#' @details This function is to read a count matrix containing a column of transcript_IDs
#' then a column of counts for every Sample_ID. It is highly recommended to first filter this for verified ORFs.
#' @param file The input count file
#' @keywords counts
#' @importFrom readr read_tsv
read_count_files <- function(file)
{
  count_data <- read_tsv(file, comment = "#", show_col_types = FALSE)

  # Check for unique Sample_IDs in column headers
  sample_ids <- names(count_data)[-1] # Exclude the "transcript_ID" column
  if (any(duplicated(sample_ids))) {
    stop("Sample_IDs in column headers need to be unique")
  }

  return(count_data)
}

#' Load the count matrix
#'
#' Load the count matrix and combine with the samplesheet
#' @details In this function, function read_count_files from this package is called.
#' Output a wider data frame for representing the count value of each transcript
#' for each fraction.
#' @param count_file the count matrix file
#' @param samplesheet sample description spreadsheet, either a tibble or a file
#' @keywords counts each fraction
#' @import dplyr tidyr
#' @export
load_counts <- function(count_file, samplesheet) {
  input_counts <- read_count_files(count_file) %>%
    pivot_longer(cols = c(-transcript_ID), names_to = "Sample_ID", values_to = "Count")

  if (is.character(samplesheet)) {
    # If samplesheet is a file path (character), read the samplesheet from the file
    samples <- read_samplesheet(samplesheet)
  } else if (is.data.frame(samplesheet)) {
    # If samplesheet is a tibble, use it directly
    samples <- samplesheet
  } else {
    stop("samplesheet must be either a file path or a tibble")
  }

  count_data <- samples %>%
    left_join(input_counts, by = "Sample_ID")
  return(count_data)
}

