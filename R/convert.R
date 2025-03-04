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
      message("Converting from 10X which lacks allele info. Choosing *01 as allele for all genes.")
    }
  } else if (frm %in% c("adaptive", "adaptivev2")) {
    lookup_f <- file.path(data_path, "lookup_from_adaptive.csv")
    if (to == "imgt" && verbose) {
      message("Converting from Adaptive to IMGT. If a gene lacks allele, will choose *01 as allele.")
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

#' Determine input columns to use for converting gene names.
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
  # Determine input columns for conversion
  if (frm == "imgt" && is.null(frm_cols)) {
    cols_from <- col_ref$tenx
    warning(paste("No column names provided for IMGT data. Using 10X column names:", paste(cols_from, collapse = ", ")))
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

  return(cols_from)
}

#' Convert T-cell receptor V, D, J, and/or C gene names from one naming convention to another.
#'
#' @param df Dataframe containing TCR gene names.
#' @param frm A string, the input format of TCR data. Must be one of "tenx", "adaptive", "adaptivev2", or "imgt".
#' @param to A string, the output format of TCR data. Must be one of "tenx", "adaptive", "adaptivev2", or "imgt".
#' @param species A string, the species. Optional; defaults to "human".
#' @param frm_cols A character vector of custom V/D/J/C gene column names. Optional; defaults to NULL.
#' @param verbose A boolean, whether to show messages. Optional; defaults to TRUE
#'
#' @return A dataframe of converted TCR data.
#' @autoglobal
#' @export
#' @examples
#' tcr_file <- get_example_path("tenx.csv")
#' df <- read.csv(tcr_file)
#' convert_gene(df, "tenx", "adaptive", verbose = FALSE)
convert_gene <- function(df, frm, to, species = "human", frm_cols = NULL, verbose = TRUE) {
  # Check that input is ok
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
    warning("Adaptive only captures VDJ genes, any C genes will become NA.")
  }

  # Load lookup table
  lookup_f <- choose_lookup(frm, to, species, verbose)
  lookup <- utils::read.csv(lookup_f, stringsAsFactors = FALSE)

  # Determine columns to use
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
        warning(paste("The input column", col, "doesn't contain any valid genes and was skipped."))
      }
    }
  }

  # Display genes we couldn't convert
  if (!all(is.na(bad_genes_all))) {
    bad_genes <- unique(bad_genes_all)
    warning(paste("These genes are not in IMGT for this species and will be replaced with NA:\n", paste(unique(bad_genes), collapse = ", ")))
  }

  # Replace columns in original dataframe
  df_out <- df
  for (col in names(new_genes)) {
    df_out[, col] <- new_genes[[col]]
  }

  # Remove row numbers column
  df_out <- subset(df_out, select = -c(id))

  # Replace NoData with NA
  df_out[df_out == "NoData"] <- NA_character_

  return(df_out)
}
