% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_lookup.R
\name{extract_imgt_genes}
\alias{extract_imgt_genes}
\title{Extract all gene names from a folder of FASTAs}
\usage{
extract_imgt_genes(data_dir)
}
\arguments{
\item{data_dir}{A string, the path to directory containing FASTA files.}
}
\value{
A dataframe of gene names.
}
\description{
\code{extract_imgt_genes()} first runs \code{parse_imgt_fasta()} on all FASTA files in
a given folder to pull out the gene names. Then it returns those names in an
alphabetically sorted dataframe.
}
\examples{
# Given a folder with FASTA files containing these headers:
#   >SomeText|TRAC*01|MoreText|
#   >SomeText|TRAV1-1*01|MoreText|
#   >SomeText|TRAV1-1*02|MoreText|
#   >SomeText|TRAV1-2*01|MoreText|
#   >SomeText|TRAV14/DV4*01|MoreText|
#   >SomeText|TRAV38-1*01|MoreText|
#   >SomeText|TRAV38-2/DV8*01|MoreText|
#   >SomeText|TRBV29-1*01|MoreText|
#   >SomeText|TRBV29-1*02|MoreText|
#   >SomeText|TRBV29/OR9-2*01|MoreText|

fastadir <- get_example_path("fasta_dir/")
extract_imgt_genes(fastadir)
}
\keyword{internal}
