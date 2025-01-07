parse_imgt_fasta <- function(infile) {
  # Extract gene names from a reference FASTA file
  #
  # Args:
  #   infile: Path to the FASTA file
  #
  # Returns:
  #   A vector of gene names
  
  # Read lines from the file
  lines <- readLines(infile)
  
  # Extract gene names from headers
  imgt_list <- sapply(lines[grep("^>", lines)], function(line) {
    strsplit(line, "\\|")[[1]][2]
  })
  
  return(imgt_list)
}

extract_imgt_genes <- function(data_dir) {
  # Extract gene names from all FASTA files in a folder
  #
  # Args:
  #   data_dir: Path to the directory containing FASTA files
  #
  # Returns:
  #   A sorted data frame with gene names
  
  # List all FASTA files
  fasta_files <- list.files(data_dir, pattern = "\\.(fa|fasta)$", full.names = TRUE)
  
  # Extract gene names from each file
  imgt <- unlist(lapply(fasta_files, parse_imgt_fasta))
  
  # Create and sort a data frame
  lookup <- data.frame(imgt = imgt, stringsAsFactors = FALSE)
  lookup_sorted <- lookup[order(lookup$imgt), , drop = FALSE]
  rownames(lookup_sorted) <- NULL
  
  return(lookup_sorted)
}

add_dash_one <- function(gene_str) {
  # Add "-01" to genes without IMGT gene-level designation
  #
  # Args:
  #   gene_str: Gene name as a string
  #
  # Returns:
  #   Modified gene name
  
  if (!grepl("-", gene_str)) {
    return(sub("\\*", "-01*", gene_str))
  }
  return(gene_str)
}

pad_single_digit <- function(gene_str) {
  # Add a leading zero to single-digit gene-level designation
  #
  # Args:
  #   gene_str: Gene name as a string
  #
  # Returns:
  #   Modified gene name
  
  return(gsub("([A-Za-z]+)(\\d)([-\\*])", "\\1\\02\\3", gene_str))
}

build_lookup_from_fastas <- function(data_dir) {
  # Create lookup tables from FASTA files in a directory
  #
  # Args:
  #   data_dir: Directory containing FASTA files
  #
  # Returns:
  #   None (saves lookup tables to files)
  
  # Extract IMGT gene names
  lookup <- extract_imgt_genes(data_dir)
  
  # Create the 10X column
  lookup$tenx <- sub("/DV", "DV", substr(lookup$imgt, 1, nchar(lookup$imgt) - 3))
  
  # Create the Adaptive column
  adaptive_replacements <- c(
    "TRAV14/DV4" = "TRAV14-1", "TRAV23/DV6" = "TRAV23-1", "TRAV29/DV5" = "TRAV29-1",
    "TRAV36/DV7" = "TRAV36-1", "TRAV38-2/DV8" = "TRAV38-2", "TRAV4-4/DV10" = "TRAV4-4/",
    "TRAV6-7/DV9" = "TRAV6-7", "TRAV13-4/DV7" = "TRAV13-4", "TRAV14D-3/DV8" = "TRAV14D-3",
    "TRAV15D-1/DV6D-1" = "TRAV15D-1", "TRAV15-1/DV6-1" = "TRAV15-1", "TRAV16D/DV11" = "TRAV16D-1",
    "TRAV21/DV12" = "TRAV21-1", "TRAV15-2/DV6-2" = "TRAV15-2", "TRAV15D-2/DV6D-2" = "TRAV15D-2",
    "TR" = "TCR", "-" = "-0", "/OR9-02" = "-or09_02"
  )
  lookup$adaptive <- lookup$imgt
  for (pattern in names(adaptive_replacements)) {
    replacement <- adaptive_replacements[[pattern]]
    lookup$adaptive <- gsub(pattern, replacement, lookup$adaptive)
  }
  lookup$adaptive <- sapply(lookup$adaptive, add_dash_one)
  lookup$adaptive <- sapply(lookup$adaptive, pad_single_digit)
  lookup$adaptivev2 <- lookup$adaptive
  
  # Set Adaptive columns to NA for constant genes
  lookup[grepl("C", lookup$imgt), c("adaptive", "adaptivev2")] <- NA
  
  # Create from_tenx and from_adaptive tables
  from_tenx <- aggregate(. ~ tenx, data = lookup, FUN = function(x) x[1])
  from_adaptive <- aggregate(. ~ adaptive, data = subset(lookup, !is.na(adaptive)), FUN = function(x) x[1])
  
  # Save to files
  write.csv(lookup, file.path(data_dir, "lookup.csv"), row.names = FALSE)
  write.csv(from_tenx, file.path(data_dir, "lookup_from_tenx.csv"), row.names = FALSE)
  write.csv(from_adaptive, file.path(data_dir, "lookup_from_adaptive.csv"), row.names = FALSE)
}
