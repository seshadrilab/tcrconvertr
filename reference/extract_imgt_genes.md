# Extract all gene names from a folder of FASTAs

`extract_imgt_genes()` first runs
[`parse_imgt_fasta()`](https://seshadrilab.github.io/tcrconvertr/reference/parse_imgt_fasta.md)
on all FASTA files in a given folder to pull out the gene names. Then it
returns those names in an alphabetically sorted dataframe.

## Usage

``` r
extract_imgt_genes(data_dir)
```

## Arguments

- data_dir:

  A string, the path to directory containing FASTA files.

## Value

A dataframe of gene names.

## Examples

``` r
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
#>               imgt
#> 1          TRAC*01
#> 2       TRAV1-1*01
#> 3       TRAV1-1*02
#> 4       TRAV1-2*01
#> 5    TRAV14/DV4*01
#> 6      TRAV38-1*01
#> 7  TRAV38-2/DV8*01
#> 8      TRBV29-1*01
#> 9      TRBV29-1*02
#> 10 TRBV29/OR9-2*01
```
