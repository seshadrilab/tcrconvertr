# Convert gene names

`convert_gene()` converts T-cell receptor (TCR) gene names between the
IMGT, 10X, and Adaptive formats. It determines the columns to convert
based on the input format (`frm`) unless specified by the user
(`frm_cols`). It returns a modified version of the input data frame with
converted gene names while preserving row order.

## Usage

``` r
convert_gene(df, frm, to, species = "human", frm_cols = NULL, verbose = TRUE)
```

## Arguments

- df:

  A dataframe containing TCR gene names.

- frm:

  A string, the input format of TCR data. Must be one of `"imgt"`,
  `"tenx"`, `"adaptive"`, or `"adaptivev2"`.

- to:

  A string, the output format of TCR data. Must be one of `"imgt"`,
  `"tenx"`, `"adaptive"`, or `"adaptivev2"`.

- species:

  A string,the species. Optional; defaults to `"human"`.

- frm_cols:

  A character vector of custom gene column names. Optional; defaults to
  `NULL`.

- verbose:

  A boolean, whether to display messages. Optional; defaults to `TRUE`.

## Value

A dataframe with converted TCR gene names.

## Details

Gene names are converted by performing a `merge` between the relevant
input columns and a species-specific lookup table containing IMGT
reference genes in all three formats.

**Behavioral Notes**

- If a gene name cannot be mapped, it is replaced with `NA` and a
  warning is raised.

- If `frm` is `'imgt'` and `frm_cols` is not provided, 10X column names
  are assumed.

- Constant (C) genes are set to `NA` when converting to Adaptive
  formats, as Adaptive does not capture constant regions.

- The input does not need to include all gene types; partial inputs
  (e.g., only V genes) are supported.

- If no values in a custom column can be mapped (e.g. a CDR3 column) it
  is skipped and a warning is raised.

**Standard Column Names**

If `frm_cols` is not provided, these column names will be used if
present:

- **IMGT**: `"v_gene"`, `"d_gene"`, `"j_gene"`, `"c_gene"`

- **10X**: `"v_gene"`, `"d_gene"`, `"j_gene"`, `"c_gene"`

- **Adaptive**: `"v_resolved"`, `"d_resolved"`, `"j_resolved"`

- **Adaptive v2**: `"vMaxResolved"`, `"dMaxResolved"`, `"jMaxResolved"`

## Examples

``` r
tcr_file <- get_example_path("tenx.csv")
df <- read.csv(tcr_file)[c("barcode", "v_gene", "j_gene", "cdr3")]
df
#>              barcode       v_gene  j_gene            cdr3
#> 1 AAACCTGAGACCACGA-1   TRAV29/DV5  TRAJ12    CAVMDSSYKLIF
#> 2 AAACCTGAGACCACGA-1 TRBV20/OR9-2 TRBJ2-1 CASSGLAGGYNEQFF
#> 3 AAACCTGAGGCTCTTA-1        TRDV2   TRDJ3 CASSGVAGGTDTQYF
#> 4 AAACCTGAGGCTCTTA-1        TRGV9   TRGJ1    CAVKDSNYQLIW
convert_gene(df, "tenx", "adaptive", verbose = FALSE)
#> Warning: Adaptive only captures VDJ genes; C genes will be NA.
#>              barcode             v_gene        j_gene            cdr3
#> 1 AAACCTGAGACCACGA-1      TCRAV29-01*01 TCRAJ12-01*01    CAVMDSSYKLIF
#> 2 AAACCTGAGACCACGA-1 TCRBV20-or09_02*01 TCRBJ02-01*01 CASSGLAGGYNEQFF
#> 3 AAACCTGAGGCTCTTA-1      TCRDV02-01*01 TCRDJ03-01*01 CASSGVAGGTDTQYF
#> 4 AAACCTGAGGCTCTTA-1      TCRGV09-01*01 TCRGJ01-01*01    CAVKDSNYQLIW
```
