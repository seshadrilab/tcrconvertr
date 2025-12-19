# Choose lookup table

`choose_lookup()` determines which CSV lookup table to use based on the
the input format (`frm`) and returns the path to that file.

## Usage

``` r
choose_lookup(frm, to, species = "human", verbose = TRUE)
```

## Arguments

- frm:

  A string, the input format of TCR data. Must be one of `"tenx"`,
  `"adaptive"`, `"adaptivev2"`, or `"imgt"`.

- to:

  A string, the output format of TCR data. Must be one of `"tenx"`,
  `"adaptive"`, `"adaptivev2"`, or `"imgt"`.

- species:

  A string, the species. Optional; defaults to `"human"`.

- verbose:

  A boolean, whether to show messages. Optional; defaults to `TRUE`

## Value

A string, the path to correct lookup table.

## Examples

``` r
choose_lookup("imgt", "adaptive")
#> [1] "/home/runner/work/_temp/Library/TCRconvertR/extdata/human/lookup.csv"
```
