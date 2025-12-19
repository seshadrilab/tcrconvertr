# Extract gene names from a reference FASTA

`parse_imgt_fasta()` extracts the second element from a "\|"-delimited
FASTA header, which will be the gene name for IMGT reference FASTAs.

## Usage

``` r
parse_imgt_fasta(infile)
```

## Arguments

- infile:

  A string, the path to FASTA file.

## Value

A character vector of gene names.

## Examples

``` r
# Given a FASTA file containing this header:
#   >SomeText|TRBV29-1*01|MoreText|
#   >SomeText|TRBV29-1*02|MoreText|
#   >SomeText|TRBV29/OR9-2*01|MoreText|

fasta <- get_example_path("fasta_dir/test_trbv.fa")
parse_imgt_fasta(fasta)
#> [1] "TRBV29-1*01"     "TRBV29-1*02"     "TRBV29/OR9-2*01"
```
