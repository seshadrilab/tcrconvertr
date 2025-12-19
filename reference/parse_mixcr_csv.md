# Extract gene names from a reference CSV for a given species

`parse_mixcr_csv()` generates mixcr gene name data for a given species.

## Usage

``` r
parse_mixcr_csv(data_dir)
```

## Arguments

- data_dir:

  A string, the path to directory containing mixcr CSV files.

## Value

A character vector of gene names.

## Examples

``` r
mixcr_dir <- get_example_path("mixcr_dir/test_mixcr")
parse_mixcr_csv(mixcr_dir)
#>           mixcr
#> 1       TCRG-C3
#> 2    TCRG-C3*00
#> 3          TRAC
#> 4       TRAC*00
#> 5     TRAV12D-3
#> 6  TRAV12D-3*00
#> 7      TRAV14-2
#> 8   TRAV14-2*00
#> 9         TRBD2
#> 10     TRBD2*00
#> 11         TRDC
#> 12      TRDC*00
#> 13        TRDD2
#> 14     TRDD2*00
#> 15        TRDJ1
#> 16     TRDJ1*00
```
