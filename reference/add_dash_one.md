# Add `-01` to gene names lacking gene-level info

Some genes just have the IMGT subgroup (e.g. TRBV2) and allele (e.g.
\*01) designation. The Adaptive format always includes an IMGT gene
(e.g. -01) designation, with "-01" as the apparent default.
`add_dash_one()` adds a default gene-level designation if it's missing.

## Usage

``` r
add_dash_one(gene_str)
```

## Arguments

- gene_str:

  A string, the gene name.

## Value

A string, the updated gene name.

## Examples

``` r
add_dash_one("TRBV2*01")
#> [1] "TRBV2-01*01"
```
