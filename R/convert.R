# Standard column names for different sources of TCR data
col_ref <- list(
  adaptive = c("v_resolved", "d_resolved", "j_resolved"),
  adaptivev2 = c("vMaxResolved", "dMaxResolved", "jMaxResolved"),
  imgt = c("v_gene", "d_gene", "j_gene", "c_gene"),
  tenx = c("v_gene", "d_gene", "j_gene", "c_gene")
)

#' Choose lookup table
#'
#' `choose_lookup()` determines which CSV lookup table to use based on the the
#' input format (`frm`) and returns the path to that file.
#'
#' @param frm A string, the input format of TCR data. Must be one of "tenx", "adaptive", "adaptivev2", or "imgt".
#' @param to A string, the output format of TCR data. Must be one of "tenx", "adaptive", "adaptivev2", or "imgt".
#' @param species A string, the species. Optional; defaults to "human".
#' @param verbose A boolean, whether to show messages. Optional; defaults to TRUE
#'
#' @return A string, the path to correct lookup table.
#' @importFrom rappdirs user_data_dir
#' @export
#' @keywords internal
#' @examples
#' choose_lookup("imgt", "adaptive")
choose_lookup <- function(frm, to, species = "human", verbose = TRUE) {
  # Determine where to find lookup tables
  if (species %in% c("human", "mouse", "rhesus")) {
    lookup_dir <- system.file("extdata", package = "TCRconvertR")
  } else {
    lookup_dir <- rappdirs::user_data_dir("TCRconvertR", "Emmma Bishop")
  }

  data_path <- file.path(lookup_dir, species)

  if (frm == "tenx") {
    lookup_f <- file.path(data_path, "lookup_from_tenx.csv")
    if (verbose) {
      message("Converting from 10X. Using *01 as allele for all genes.")
    }
  } else if (frm %in% c("adaptive", "adaptivev2")) {
    lookup_f <- file.path(data_path, "lookup_from_adaptive.csv")
    if (to == "imgt" && verbose) {
      message("Converting from Adaptive to IMGT. Using *01 for genes lacking alleles.")
    }
  } else {
    lookup_f <- file.path(data_path, "lookup.csv")
  }

  # Check if file exists
  if (file.exists(lookup_f)) {
    return(lookup_f)
  } else {
    stop("Lookup table not found, please ensure reference files are available.")
  }
}

#' Determine input columns to use
#'
#' `which_frm_cols()` determines the columns that are expected to hold gene
#' name information in the input file based on the input format (`frm`). It
#' returns a vector of those column names.
#'
#' @param df Dataframe containing TCR gene names.
#' @param frm A string, the input format of TCR data. Must be one of "tenx", "adaptive", "adaptivev2", or "imgt".
#' @param frm_cols A character vector, the custom column names to use.
#' @param verbose A boolean, whether to show messages. Optional; defaults to TRUE
#'
#' @return A character vector, column names to use.
#' @export
#' @keywords internal
#' @examples
#' tcr_file <- get_example_path("tenx.csv")
#' df <- read.csv(tcr_file)
#' which_frm_cols(df, "tenx")
which_frm_cols <- function(df, frm, frm_cols = NULL, verbose = TRUE) {
  if (frm == "imgt" && is.null(frm_cols)) {
    cols_from <- col_ref$tenx
    warning(paste("No column names for IMGT data. Using 10X columns:", paste(cols_from, collapse = ", ")))
  } else if (!is.null(frm_cols)) {
    missing_cols <- setdiff(frm_cols, colnames(df))
    if (length(missing_cols) > 0) {
      stop(paste("These columns are not in the input dataframe:", paste(missing_cols, collapse = ", ")))
    } else {
      cols_from <- frm_cols
      if (verbose) {
        message(paste("Using custom column names:", paste(cols_from, collapse = ", ")))
      }
    }
  } else {
    cols_from <- col_ref[[frm]]
  }

  cols_from
}

#' Convert gene names
#'
#' @description
#' `convert_gene()` converts T-cell receptor (TCR) gene names between the IMGT,
#' 10X, and Adaptive formats. It determines the columns to convert based
#' on the input format (`frm`) unless specified by the user (`frm_cols`). It
#' returns a modified version of the input data frame with converted gene names
#' while preserving row order.
#'
#' @details
#' Gene names are converted by performing a `merge` between the relevant
#' input columns and a species-specific lookup table containing IMGT reference
#' genes in all three formats.
#'
#' **Behavioral Notes**
#' - If a gene name cannot be mapped, it is replaced with `NA` and a warning is
#' raised.
#' - If the input format is IMGT and `frm_cols` is not provided, 10X column
#' names are assumed.
#' - Constant (C) genes are set to `NA` when converting to Adaptive formats,
#' as Adaptive does not capture constant regions.
#' - The input does not need to include all gene types. Partial inputs
#' (e.g., only V genes) are supported.
#' - If no values in a custom column can be mapped (e.g. a CDR3 column) it is
#' skipped and a warning is raised.
#'
#' **Standard Column Names**
#'
#' If `frm_cols` is not provided, these column names will be used if present:
#' - **IMGT**: `"v_gene"`, `"d_gene"`, `"j_gene"`, `"c_gene"`
#' - **10X**: `"v_gene"`, `"d_gene"`, `"j_gene"`, `"c_gene"`
#' - **Adaptive**: `"v_resolved"`, `"d_resolved"`, `"j_resolved"`
#' - **Adaptive v2**: `"vMaxResolved"`, `"dMaxResolved"`, `"jMaxResolved"`
#'
#' @param df A dataframe containing TCR gene names.
#' @param frm A string, the input format of TCR data. Must be one of `"imgt"`, `"tenx"`, `"adaptive"`, or `"adaptivev2"`.
#' @param to A string, the output format of TCR data. Must be one of `"imgt"`, `"tenx"`, `"adaptive"`, or `"adaptivev2"`.
#' @param species A string, the species folder name under `tcrconvertr/inst/extdata/`. Optional; defaults to `"human"`.
#' @param frm_cols A character vector of customgene column names. Optional; defaults to `NULL` (auto-detect).
#' @param verbose A boolean, whether to display messages. Optional; defaults to `TRUE`.
#'
#' @return A dataframe with converted TCR gene names.
#' @autoglobal
#' @export
#' @examples
#' tcr_file <- get_example_path("tenx.csv")
#' df <- read.csv(tcr_file)[c("barcode", "v_gene", "j_gene", "cdr3")]
#' convert_gene(df, "tenx", "adaptive", verbose = FALSE)
convert_gene <- function(df, frm, to, species = "human", frm_cols = NULL, verbose = TRUE) {
  if (frm == to) {
    stop('"frm" and "to" formats should be different.')
  }
  if (!is.data.frame(df)) {
    stop("Input is not a data.frame.")
  }
  if (nrow(df) == 0) {
    stop("Input data is empty.")
  }
  if (to %in% c("adaptive", "adaptivev2")) {
    warning("Adaptive captures only VDJ genes; C genes will be NA.")
  }

  # Load lookup table and determine input columns
  lookup_f <- choose_lookup(frm, to, species, verbose)
  lookup <- utils::read.csv(lookup_f, stringsAsFactors = FALSE)
  cols_from <- which_frm_cols(df, frm, frm_cols, verbose)

  # Add column of row numbers so we can keep order straight
  df$id <- seq_len(nrow(df))

  # Loop over gene columns, doing a merge to get converted gene names
  new_genes <- list()
  bad_genes_all <- c()

  for (col in cols_from) {
    if (col %in% colnames(df)) {
      merged <- merge(df[, c(col, "id"), drop = FALSE], lookup, by.x = col, by.y = frm, all.x = TRUE)
      merged <- merged[order(merged$id), ]
      good_genes <- merged[, to]
      # Note genes where the merge produced an NA on the 'to' format side
      new_bad_genes <- merged[is.na(merged[, to]), col]
      # We don't expect the entire column of genes to be empty.
      if (length(new_bad_genes) < length(good_genes)) {
        new_genes[[col]] <- good_genes
        bad_genes_all <- c(bad_genes_all, new_bad_genes)
      } else {
        warning(paste("Column", col, "has no valid genes and was skipped."))
      }
    }
  }

  # Display genes we couldn't convert
  if (!all(is.na(bad_genes_all))) {
    bad_genes <- unique(bad_genes_all)
    warning(paste(
      "These genes are not in IMGT for this species and will be replaced with NA:\n",
      paste(unique(bad_genes), collapse = ", ")
    ))
  }

  # Build output
  df_out <- df
  for (col in names(new_genes)) {
    df_out[, col] <- new_genes[[col]]
  }
  df_out <- subset(df_out, select = -c(id)) # Remove row number column
  df_out[df_out == "NoData"] <- NA_character_

  df_out
}
