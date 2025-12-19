# Add a `0` to single-digit gene-level designation

`pad_single_digit()` takes a gene name and ensures that any single-digit
number following a sequence of letters is padded with a leading zero.
This is to match the Adaptive format.

## Usage

``` r
pad_single_digit(gene_str)
```

## Arguments

- gene_str:

  A string, the gene name.

## Value

A string, the updated gene name.

## Examples

``` r
pad_single_digit("TCRBV1-2")
#> [1] "TCRBV01-2"
```
