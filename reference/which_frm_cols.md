# Determine input columns to use

`which_frm_cols()` determines the columns that are expected to hold gene
name information in the input file based on the input format (`frm`). It
returns a vector of those column names.

## Usage

``` r
which_frm_cols(df, frm, frm_cols = NULL, verbose = TRUE)
```

## Arguments

- df:

  Dataframe containing TCR gene names.

- frm:

  A string, the input format of TCR data. Must be one of `"tenx"`,
  `"adaptive"`, `"adaptivev2"`, or `"imgt"`.

- frm_cols:

  A character vector, the custom column names to use.

- verbose:

  A boolean, whether to show messages. Optional; defaults to `TRUE`

## Value

A character vector, column names to use.

## Examples

``` r
tcr_file <- get_example_path("tenx.csv")
df <- read.csv(tcr_file)
which_frm_cols(df, "tenx")
#> [1] "v_gene" "d_gene" "j_gene" "c_gene"
```
