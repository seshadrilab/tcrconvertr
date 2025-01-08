# Load required packages
library(tools)  # For file path manipulation
library(utils)  # For reading and writing data files

# Set up logging
log_message <- function(level, message) {
  cat(paste0("[", level, "] ", message, "\n"))
}

# Standard column names for different sources of TCR data
col_ref <- list(
  adaptive = c("v_resolved", "d_resolved", "j_resolved"),
  adaptivev2 = c("vMaxResolved", "dMaxResolved", "jMaxResolved"),
  imgt = c("v_gene", "d_gene", "j_gene", "c_gene"),
  tenx = c("v_gene", "d_gene", "j_gene", "c_gene")
)

#' Determine which lookup table to use and get filepath.
#'
#' @param frm A string, the input format of TCR data. Must be one of "tenx", "adaptive", "adaptivev2", or "imgt".
#' @param to A string, the output format of TCR data. Must be one of "tenx", "adaptive", "adaptivev2", or "imgt".
#' @param species A string, the species folder name under `tcrconvertr/inst/extdata/`. Optional; defaults to "human".
#'
#' @return A string, the path to correct lookup table.
choose_lookup <- function(frm, to, species = "human") {
  # Determine the lookup table path
  data_path <- system.file("data", species, package = "tcrconvert")
  
  if (frm == "tenx") {
    lookup_f <- file.path(data_path, "lookup_from_tenx.csv")
    log_message("WARNING", "Converting from 10X which lacks allele info. Choosing *01 as allele for all genes.")
  } else if (frm %in% c("adaptive", "adaptivev2")) {
    lookup_f <- file.path(data_path, "lookup_from_adaptive.csv")
    if (to == "imgt") {
      log_message("WARNING", "Converting from Adaptive to IMGT. If a gene lacks allele, will choose *01 as allele.")
    }
  } else {
    lookup_f <- file.path(data_path, "lookup.csv")
  }
  
  # Check if file exists
  if (file.exists(lookup_f)) {
    return(lookup_f)
  } else {
    log_message("ERROR", "Lookup table not found, please ensure reference files are available.")
    stop("FileNotFoundError")
  }
}

#' Determine input columns to use for converting gene names.
#'
#' @param df Dataframe containing TCR gene names.
#' @param frm A string, the input format of TCR data. Must be one of "tenx", "adaptive", "adaptivev2", or "imgt".
#' @param frm_cols A character vector, the custom column names to use.
#'
#' @return A character vector, column names to use.
which_frm_cols <- function(df, frm, frm_cols = NULL) {
  # Determine input columns for conversion
  if (frm == "imgt" && is.null(frm_cols)) {
    cols_from <- col_ref$tenx
    log_message("INFO", paste("No column names provided for IMGT data. Using 10X column names:", paste(cols_from, collapse = ", ")))
  } else if (!is.null(frm_cols)) {
    missing_cols <- setdiff(frm_cols, colnames(df))
    if (length(missing_cols) > 0) {
      log_message("ERROR", paste("These columns are not in the input dataframe:", paste(missing_cols, collapse = ", ")))
      stop("ValueError")
    } else {
      cols_from <- frm_cols
      log_message("INFO", paste("Using custom column names:", paste(cols_from, collapse = ", ")))
    }
  } else {
    cols_from <- col_ref[[frm]]
  }
  
  return(cols_from)
}

#' Convert T-cell receptor V, D, J, and/or C gene names from one naming convention to another.
#'
#' @param df Dataframe containing TCR gene names.
#' @param frm A string, the input format of TCR data. Must be one of "tenx", "adaptive", "adaptivev2", or "imgt".
#' @param to A string, the output format of TCR data. Must be one of "tenx", "adaptive", "adaptivev2", or "imgt".
#' @param species A string, the species folder name under `tcrconvertr/inst/extdata/`. Optional; defaults to "human".
#' @param frm_cols A character vector of custom V/D/J/C gene column names. Optional; defaults to NULL.
#' @param quiet A boolean, whether to suppress warning messages. Optional; defaults to FALSE.
#'
#' @return A dataframe of converted TCR data.
#' @export
convert_gene <- function(df, frm, to, species = "human", frm_cols = NULL, quiet = FALSE) {
  # Check input and output formats
  if (frm == to) {
    log_message("ERROR", '"frm" and "to" formats should be different.')
    stop("ValueError")
  }
  if (nrow(df) == 0) {
    log_message("ERROR", "Input data is empty.")
    stop("ValueError")
  }
  
  if (to %in% c("adaptive", "adaptivev2")) {
    log_message("WARNING", "Adaptive only captures VDJ genes, any C genes will become NA.")
  }
  
  # Load lookup table
  lookup_f <- choose_lookup(frm, to, species)
  lookup <- read.csv(lookup_f, stringsAsFactors = FALSE)
  
  # Determine columns to use
  cols_from <- which_frm_cols(df, frm, frm_cols)
  
  # Perform conversion
  new_genes <- list()
  bad_genes <- c()
  
  for (col in cols_from) {
    if (col %in% colnames(df)) {
      merged <- merge(df[, col, drop = FALSE], lookup, by.x = col, by.y = frm, all.x = TRUE)
      new_genes[[col]] <- merged[, to]
      bad_genes <- c(bad_genes, merged[is.na(merged[, to]), col])
    }
  }
  
  if (length(bad_genes) > 0) {
    log_message("WARNING", paste("These genes are not in IMGT for this species and will be replaced with NA:", paste(unique(bad_genes), collapse = ", ")))
  }
  
  # Replace columns in original dataframe
  df_out <- df
  for (col in names(new_genes)) {
    df_out[, col] <- new_genes[[col]]
  }
  
  return(df_out)
}

# Example usage of convert_gene function
# Uncomment to test:
# input_df <- read.csv("input.csv", stringsAsFactors = FALSE)
# converted_df <- convert_gene(input_df, frm = "tenx", to = "adaptive", species = "human")
# write.csv(converted_df, "output.csv", row.names = FALSE)
