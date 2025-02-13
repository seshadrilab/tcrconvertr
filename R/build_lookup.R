#' Extract gene names from a reference FASTA.
#'
#' @param infile A string, the path to FASTA file.
#'
#' @return A character vector of gene names.
#' @examples
#' # Given a FASTA file containing this header:
#' #   >SomeText|TRBV29/OR9-2*01|MoreText|
#' #   >SomeText|TRBVA/OR9-2*01|MoreText|
#'
#' fasta <- get_example_path("fasta_dir/test_trbv.fa")
#' TCRconvertR:::parse_imgt_fasta(fasta)
parse_imgt_fasta <- function(infile) {
  lines <- readLines(infile)
  
  # Extract gene names from headers
  imgt_list <- sapply(lines[grep("^>", lines)], function(line) {
    strsplit(line, "\\|")[[1]][2]
  }, USE.NAMES = FALSE)
  
  return(imgt_list)
}

#' Extract gene names from all reference FASTA files in a folder.
#'
#' @param data_dir A string, the path to directory containing FASTA files.
#'
#' @return A dataframe of gene names.
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
#' fastadir <- get_example_path('fasta_dir/')
#' TCRconvertR:::extract_imgt_genes(fastadir)
extract_imgt_genes <- function(data_dir) {
  # List all FASTA files
  fasta_files <- list.files(data_dir, pattern = "\\.(fa|fasta)$", full.names = TRUE)
  
  # Extract gene names from each file
  imgt <- unlist(lapply(fasta_files, parse_imgt_fasta))
  
  # Create and sort a data frame
  lookup <- data.frame(imgt = imgt, stringsAsFactors = FALSE)
  lookup_sorted <- lookup[order(lookup[["imgt"]]), , drop = FALSE]
  rownames(lookup_sorted) <- NULL
  
  return(lookup_sorted)
}

#' Add a `-01` to genes without IMGT gene-level designation.
#'
#' @param gene_str A string, the gene name.
#'
#' @return A string, the updated gene name.
#' @examples
#' TCRconvertR:::add_dash_one('TRBV2*01')
add_dash_one <- function(gene_str) {
  if (!length(gene_str) == 1) {
    stop("Attempting to add '-01' to more than one gene name at a time.")
  }
  if (!grepl("-", gene_str)) {
    return(sub("\\*", "-01*", gene_str))
  }
  return(gene_str)
}

#' Add a zero to single-digit gene-level designation in gene names.
#'
#' @param gene_str A string, the gene name.
#'
#' @return A string, the updated gene name.
#' @examples
#' TCRconvertR:::pad_single_digit('TCRBV1-2')
pad_single_digit <- function(gene_str) {
  return(gsub("([A-Za-z]+)(\\d)([-\\*])", "\\10\\2\\3", gene_str))
}

#' Create lookup tables from reference FASTAs
#'
#' Create lookup tables within in a given directory that contains FASTA files:
#'    - lookup.csv
#'    - lookup_from_tenx.csv
#'    - lookup_from_adaptive.csv
#'
#' @param data_dir A string, the directory containing FASTA files.
#'
#' @return Nothing.
#' @autoglobal
#' @export
#' @examples
#' # For the example, create and use a temporary folder
#' fastadir <- file.path(tempdir(), "tcrconvertr_tmp")
#' dir.create(fastadir)
#' trav <- get_example_path('fasta_dir/test_trav.fa')
#' trbv <- get_example_path('fasta_dir/test_trbv.fa')
#' file.copy(c(trav, trbv), fastadir)
#'
#' # Build lookup tables
#' build_lookup_from_fastas(fastadir)
#'
#' # Clean up temporary folder
#' unlink(fastadir, recursive = TRUE)
build_lookup_from_fastas <- function(data_dir) {
  lookup <- extract_imgt_genes(data_dir)
  
  # Create the 10X column by removing allele info (e.g. *01) and slash from "/DV".
  # Do this by substituting "DV" for "/DV" in what's left of the TCR gene name
  # after removing the last three characters
  lookup[["tenx"]] <- sub("/DV", "DV", substr(lookup[["imgt"]], 1,
                                              nchar(lookup[["imgt"]]) - 3))
  
  # Create Adaptive columns by adding letters, 0's, removing /DV and renaming /OR
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
  lookup[["adaptive"]] <- lookup[["imgt"]]
  for (pattern in names(adaptive_replacements)) {
    replacement <- adaptive_replacements[[pattern]]
    lookup[["adaptive"]] <- gsub(pattern, replacement, lookup[["adaptive"]])
  }
  lookup[["adaptive"]] <- sapply(lookup[["adaptive"]], add_dash_one)
  lookup[["adaptive"]] <- sapply(lookup[["adaptive"]], pad_single_digit)
  lookup[["adaptivev2"]] <- lookup[["adaptive"]]
  
  # Set Adaptive columns to 'NoData' for constant genes (Adaptive only captures VDJ)
  # 'NoData' gets converted to NA by convert.convert_gene()
  lookup[grepl("C", lookup[["imgt"]]), c("adaptive", "adaptivev2")] <- "NoData"
  
  # If converting from 10X will just need the first *01 allele
  from_tenx <- stats::aggregate(. ~ tenx, data = lookup, FUN = function(x) x[1])
  
  # Make table for Adaptive genes with or without allele
  lookup2 <- subset(lookup, !grepl("NoData", lookup$adaptivev2))
  from_adapt <- lookup2[c("adaptivev2", "imgt", "tenx")]
  from_adapt["adaptive"] <-  substr(from_adapt[["adaptivev2"]], 1,
                                    nchar(from_adapt[["adaptivev2"]]) - 3)
  from_adapt <- stats::aggregate(. ~ adaptive,
                                 data = subset(from_adapt, select = -adaptivev2),
                                 FUN = function(x) x[1])
  from_adapt["adaptivev2"] <- from_adapt["adaptive"]
  from_adaptive <- rbind(lookup2, from_adapt)
  from_adaptive <- from_adaptive[, c("adaptive", "adaptivev2", "imgt", "tenx")]

  # Remove any duplicate rows and save
  lookup <- lookup[!duplicated(lookup), ]
  from_tenx <- from_tenx[!duplicated(from_tenx), ]
  from_adaptive <- from_adaptive[!duplicated(from_adaptive), ]

  # Save to files
  utils::write.csv(lookup, file.path(data_dir, "lookup.csv"), row.names = FALSE)
  utils::write.csv(from_tenx, file.path(data_dir, "lookup_from_tenx.csv"), row.names = FALSE)
  utils::write.csv(from_adaptive, file.path(data_dir, "lookup_from_adaptive.csv"), row.names = FALSE)
}
