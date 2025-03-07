# Special cases in gene nomenclature
adaptive_replacements <- c(
  "TRAV14/DV4" = "TRAV14-1",
  "TRAV23/DV6" = "TRAV23-1",
  "TRAV29/DV5" = "TRAV29-1",
  "TRAV36/DV7" = "TRAV36-1",
  "TRAV38-2/DV8" = "TRAV38-2",
  "TRAV4-4/DV10" = "TRAV4-4/",
  "TRAV6-7/DV9" = "TRAV6-7",
  "TRAV13-4/DV7" = "TRAV13-4",
  "TRAV14D-3/DV8" = "TRAV14D-3",
  "TRAV15D-1/DV6D-1" = "TRAV15D-1",
  "TRAV15-1/DV6-1" = "TRAV15-1",
  "TRAV16D/DV11" = "TRAV16D-1",
  "TRAV21/DV12" = "TRAV21-1",
  "TRAV15-2/DV6-2" = "TRAV15-2",
  "TRAV15D-2/DV6D-2" = "TRAV15D-2",
  "TR" = "TCR", "-" = "-0", "/OR9-02" = "-or09_02"
)

#' Extract gene names from a reference FASTA
#'
#' `parse_imgt_fasta()` extracts the second element from a "|"-delimited FASTA
#' header, which will be the gene name for IMGT reference FASTAs.
#'
#' @param infile A string, the path to FASTA file.
#'
#' @return A character vector of gene names.
#' @export
#' @keywords internal
#' @examples
#' # Given a FASTA file containing this header:
#' #   >SomeText|TRBV29/OR9-2*01|MoreText|
#' #   >SomeText|TRBVA/OR9-2*01|MoreText|
#'
#' fasta <- get_example_path("fasta_dir/test_trbv.fa")
#' parse_imgt_fasta(fasta)
parse_imgt_fasta <- function(infile) {
  lines <- readLines(infile)

  imgt_list <- vapply(lines[grep("^>", lines)], function(line) {
    strsplit(line, "\\|")[[1]][2]
  }, FUN.VALUE = character(1), USE.NAMES = FALSE)

  imgt_list
}

#' Extract all gene names from a folder of FASTAs
#'
#' `extract_imgt_genes()` first runs `parse_imgt_fasta()` on all FASTA files in
#' a given folder to pull out the gene names. Then it returns those names in an
#' alphabetically sorted dataframe.
#'
#' @param data_dir A string, the path to directory containing FASTA files.
#'
#' @return A dataframe of gene names.
#' @export
#' @keywords internal
#' @examples
#' # Given a folder with FASTA files containing these headers:
#' # >SomeText|TRAC*01|MoreText|
#' # >SomeText|TRAV1-1*01|MoreText|
#' # >SomeText|TRAV1-1*02|MoreText|
#' # >SomeText|TRAV14/DV4*01|MoreText|
#' # >SomeText|TRAV38-2/DV8*01|MoreText|
#' # >SomeText|TRBV29/OR9-2*01|MoreText|
#' # >SomeText|TRBVA/OR9-2*01|MoreText|
#'
#' fastadir <- get_example_path("fasta_dir/")
#' extract_imgt_genes(fastadir)
extract_imgt_genes <- function(data_dir) {
  fasta_files <- list.files(data_dir, pattern = "\\.(fa|fasta)$", full.names = TRUE)
  imgt <- unlist(lapply(fasta_files, parse_imgt_fasta))

  # Create and sort output data frame
  lookup <- data.frame(imgt = imgt, stringsAsFactors = FALSE)
  lookup_sorted <- lookup[order(lookup[["imgt"]]), , drop = FALSE]
  rownames(lookup_sorted) <- NULL

  lookup_sorted
}

#' Add `-01` to gene names lacking gene-level info
#'
#' Some genes just have the IMGT subgroup (e.g. TRBV2) and allele (e.g. *01)
#' designation. The Adaptive format always includes an IMGT gene (e.g. -01)
#' designation, with "-01" as the apparent default. `add_dash_one()` adds a
#' default gene-level designation if it's missing.
#'
#' @param gene_str A string, the gene name.
#'
#' @return A string, the updated gene name.
#' @export
#' @keywords internal
#' @examples
#' add_dash_one("TRBV2*01")
add_dash_one <- function(gene_str) {
  if (!length(gene_str) == 1) {
    stop("Attempting to add '-01' to more than one gene name at a time.")
  }
  if (!grepl("-", gene_str)) {
    return(sub("\\*", "-01*", gene_str))
  }

  gene_str
}

#' Add a `0` to single-digit gene-level designation
#'
#' `pad_single_digit()` takes a gene name and ensures that any single-digit
#' number following a sequence of letters is padded with a leading zero.
#' This is to match the Adaptive format.
#'
#' @param gene_str A string, the gene name.
#'
#' @return A string, the updated gene name.
#' @export
#' @keywords internal
#' @examples
#' pad_single_digit("TCRBV1-2")
pad_single_digit <- function(gene_str) {
  gsub("([A-Za-z]+)(\\d)([-\\*])", "\\10\\2\\3", gene_str)
}

#' Save a lookup table to a CSV file
#'
#' `save_lookup()` saves a data frame as a CSV file (without row names) in the
#' specified directory.
#'
#' @param df A data frame containing the lookup table data.
#' @param savedir A string, the path to the save directory.
#' @param name A string, the file name (should end in `.csv`).
#'
#' @return Nothing
#' @export
#' @keywords internal
#' @examples
#' # Create a temp save directory and load an example
#' save_dir <- file.path(tempdir(), "TCRconvertR_tmp")
#' dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
#' dat <- read.csv(get_example_path("fasta_dir/lookup.csv"))
#'
#' save_lookup(dat, save_dir, "newlookup.csv")
#'
#' # Clean up temporary folder
#' unlink(save_dir, recursive = TRUE)
save_lookup <- function(df, savedir, name) {
  # Ensure valid inputs
  if (!dir.exists(savedir)) {
    dir.create(savedir, showWarnings = FALSE, recursive = TRUE)
  }
  if (!is.data.frame(df)) {
    stop("'df' must be a data frame")
  }

  file_path <- file.path(savedir, name)
  utils::write.csv(df, file_path, row.names = FALSE)
}

#' Create lookup tables
#'
#' @description
#' `build_lookup_from_fastas()` processes IMGT reference FASTA files in a given
#' folder to generate lookup tables used for making gene name conversions. It
#' extracts all gene names and transforms them into 10X and Adaptive formats
#' following predefined conversion rules. The resulting files are created:
#'
#' - `lookup.csv`: IMGT gene names and their 10X and Adaptive equivalents.
#' - `lookup_from_tenx.csv`: Gene names aggregated by their 10X identifiers, with one representative allele (`*01`) for each.
#' - `lookup_from_adaptive.csv`: Adaptive gene names, with or without alleles, and their IMGT and 10X equivalents.
#'
#' The files are stored in a given subfolder (`species`) within the appropriate
#' application folder via `rappdirs`. For example:
#'    - MacOS: ``~/Library/Application Support/<AppName>``
#'    - Windows: ``C:\Documents and Settings\<User>\Application Data\Local Settings\<AppAuthor>\<AppName>``
#'    - Linux: ``~/.local/share/<AppName>``
#'
#' If a folder named `species` already exists in that location, it will be replaced.
#'
#' @details
#' Key transformations from IMGT:
#' - **10X:**
#'     - Remove allele information (e.g., `*01`) and modify `/DV` occurrences.
#' - **Adaptive:**
#'     - Apply renaming rules, such as adding gene-level designations and zero-padding single-digit numbers.
#'     - Convert constant genes to `"NoData"` (Adaptive only captures VDJ) which become `NA` after the merge in `convert_gene()`.
#'
#' @param data_dir A string, the directory containing FASTA files.
#' @param species A string, the name of species that will be used when running TCRconvert with these lookup tables.
#'
#' @return A string, path to new lookup directory
#' @importFrom rappdirs user_data_dir
#' @autoglobal
#' @export
#' @examples
#' # For the example, create and use a temporary folder
#' fastadir <- file.path(tempdir(), "TCRconvertR_tmp")
#' dir.create(fastadir)
#' trav <- get_example_path("fasta_dir/test_trav.fa")
#' trbv <- get_example_path("fasta_dir/test_trbv.fa")
#' file.copy(c(trav, trbv), fastadir)
#'
#' # Build lookup tables
#' build_lookup_from_fastas(fastadir, "rabbit")
#'
#' # Clean up temporary folder
#' unlink(fastadir, recursive = TRUE)
build_lookup_from_fastas <- function(data_dir, species) {
  # Check that species can be a valid folder name
  forbidden_char <- "[/\\\\:*?\"<>|~`\n\t]"
  if (grepl(forbidden_char, species)) {
    sanitized <- gsub("[/\\\\:*?\"<>|~`\n\t]", "_", species)
    stop(paste(
      "Proposed folder name", species, "contains invalid characters.\n",
      " Suggestion:", sanitized
    ))
  }

  # Get the user data directory for saving lookup tables
  user_dir <- rappdirs::user_data_dir("TCRconvertR", "Emmma Bishop")
  save_dir <- file.path(user_dir, species)
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE) # Ensure species subdirectory exists

  lookup <- extract_imgt_genes(data_dir)

  # Create the 10X column
  lookup[["tenx"]] <- sub("/DV", "DV", substr(
    lookup[["imgt"]], 1,
    nchar(lookup[["imgt"]]) - 3
  ))

  # Create Adaptive columns
  lookup[["adaptive"]] <- lookup[["imgt"]]
  for (pattern in names(adaptive_replacements)) {
    replacement <- adaptive_replacements[[pattern]]
    lookup[["adaptive"]] <- gsub(pattern, replacement, lookup[["adaptive"]])
  }
  lookup[["adaptive"]] <- vapply(lookup[["adaptive"]], add_dash_one,
    FUN.VALUE = character(1), USE.NAMES = FALSE
  )
  lookup[["adaptive"]] <- vapply(lookup[["adaptive"]], pad_single_digit,
    FUN.VALUE = character(1), USE.NAMES = FALSE
  )
  lookup[["adaptivev2"]] <- lookup[["adaptive"]]
  lookup[grepl("C", lookup[["imgt"]]), c("adaptive", "adaptivev2")] <- "NoData"

  # If converting from 10X just need the *01 allele
  from_tenx <- stats::aggregate(. ~ tenx, data = lookup, FUN = function(x) x[1])

  # Make table for Adaptive genes with and without allele
  lookup2 <- subset(lookup, !grepl("NoData", lookup$adaptivev2))
  from_adapt <- lookup2[c("adaptivev2", "imgt", "tenx")]
  from_adapt["adaptive"] <- substr(
    from_adapt[["adaptivev2"]], 1,
    nchar(from_adapt[["adaptivev2"]]) - 3
  )
  from_adapt <- stats::aggregate(. ~ adaptive,
    data = subset(from_adapt, select = -adaptivev2),
    FUN = function(x) x[1]
  )
  from_adapt["adaptivev2"] <- from_adapt["adaptive"]
  from_adaptive <- rbind(lookup2, from_adapt)
  from_adaptive <- from_adaptive[, c("adaptive", "adaptivev2", "imgt", "tenx")]

  # Remove duplicate rows
  lookup <- lookup[!duplicated(lookup), ]
  from_tenx <- from_tenx[!duplicated(from_tenx), ]
  from_adaptive <- from_adaptive[!duplicated(from_adaptive), ]

  # Save to files
  message("Writing lookup tables to: ", save_dir)
  save_lookup(lookup, save_dir, "lookup.csv")
  save_lookup(from_tenx, save_dir, "lookup_from_tenx.csv")
  save_lookup(from_adaptive, save_dir, "lookup_from_adaptive.csv")

  save_dir
}
