# Get started

## Installation

Install from
[CRAN](https://cran.r-project.org/web/packages/TCRconvertR/index.html):

``` r
install.packages("TCRconvertR")
```

## Basic usage

#### 1. Load TCRs into a data frame

Examples of files you may want to load:

- **10X**: `filtered_contig_annotations.csv`
- **Adaptive**: `Sample_TCRB.tsv`
- **IMGT**: Output from `MiXCR` or other tools

``` r
library(TCRconvertR)

tcr_file <- get_example_path("tenx.csv") # Using built-in example file
tcrs <- read.csv(tcr_file)[c("barcode", "v_gene", "j_gene", "cdr3")]
tcrs
#>              barcode       v_gene  j_gene            cdr3
#> 1 AAACCTGAGACCACGA-1   TRAV29/DV5  TRAJ12    CAVMDSSYKLIF
#> 2 AAACCTGAGACCACGA-1 TRBV20/OR9-2 TRBJ2-1 CASSGLAGGYNEQFF
#> 3 AAACCTGAGGCTCTTA-1        TRDV2   TRDJ3 CASSGVAGGTDTQYF
#> 4 AAACCTGAGGCTCTTA-1        TRGV9   TRGJ1    CAVKDSNYQLIW
```

#### 2. Convert

``` r
new_tcrs <- convert_gene(tcrs, frm = "tenx", to = "adaptive")
#> Warning in convert_gene(tcrs, frm = "tenx", to = "adaptive"): Adaptive only
#> captures VDJ genes; C genes will be NA.
#> Converting from 10X. Using *01 as allele for all genes.
new_tcrs
#>              barcode             v_gene        j_gene            cdr3
#> 1 AAACCTGAGACCACGA-1      TCRAV29-01*01 TCRAJ12-01*01    CAVMDSSYKLIF
#> 2 AAACCTGAGACCACGA-1 TCRBV20-or09_02*01 TCRBJ02-01*01 CASSGLAGGYNEQFF
#> 3 AAACCTGAGGCTCTTA-1      TCRDV02-01*01 TCRDJ03-01*01 CASSGVAGGTDTQYF
#> 4 AAACCTGAGGCTCTTA-1      TCRGV09-01*01 TCRGJ01-01*01    CAVKDSNYQLIW
```

> Tip: Suppress messages by setting `verbose = FALSE`. Warnings and
> errors will still appear.

> Tip: If your Adaptive data lacks `x_resolved`/`xMaxResolved` columns,
> create them yourself by combining the `x_gene`/`xGeneName` and
> `x_allele`/`xGeneAllele` columns. See the FAQs.

## AIRR data

Supply the standard AIRR gene column names to `frm_cols`:

``` r
new_airr <- convert_gene(airr, frm = "imgt", to = "adaptive", 
                         frm_cols = c('v_call', 'd_call', 'j_call', 'c_call'))
```

## Custom column names

By default, `TCRconvertR` assumes these column names based on the input
nomenclature (`frm`):

- `frm = 'imgt'` : `c('v_gene', 'd_gene', 'j_gene', 'c_gene')`
- `frm = 'tenx'` : `c('v_gene', 'd_gene', 'j_gene', 'c_gene')`
- `frm = 'adaptive'` : `c('v_resolved', 'd_resolved', 'j_resolved')`
- `frm = 'adaptivev2'` :
  `c('vMaxResolved', 'dMaxResolved', 'jMaxResolved')`

You can override these columns using `frm_cols`:

**1. Load 10X data with custom column names**

``` r
custom_file <- get_example_path("customcols.csv")

custom <- read.csv(custom_file)
custom
#>   myVgene myDgene myJgene myCgene          myCDR3 antigen
#> 1 TRAV1-2   TRBD1  TRAJ12    TRAC    CAVMDSSYKLIF     Flu
#> 2 TRBV6-1   TRBD2 TRBJ2-1   TRBC2 CASSGLAGGYNEQFF     Flu
#> 3 TRBV6-4   TRBD2 TRBJ2-3   TRBC2 CASSGVAGGTDTQYF     CMV
#> 4 TRAV1-2   TRBD1  TRAJ33    TRAC    CAVKDSNYQLIW     CMV
#> 5   TRBV2   TRBD1 TRBJ1-2   TRBC1   CASNQGLNYGYTF     CMV
```

**2. Specify names using `frm_cols` and convert to IMGT**

``` r
custom_new <- convert_gene(
  custom,
  frm = "tenx",
  to = "imgt",
  verbose = FALSE,
  frm_cols = c("myVgene", "myDgene", "myJgene", "myCgene"),
)
custom_new
#>      myVgene  myDgene    myJgene  myCgene          myCDR3 antigen
#> 1 TRAV1-2*01 TRBD1*01  TRAJ12*01  TRAC*01    CAVMDSSYKLIF     Flu
#> 2 TRBV6-1*01 TRBD2*01 TRBJ2-1*01 TRBC2*01 CASSGLAGGYNEQFF     Flu
#> 3 TRBV6-4*01 TRBD2*01 TRBJ2-3*01 TRBC2*01 CASSGVAGGTDTQYF     CMV
#> 4 TRAV1-2*01 TRBD1*01  TRAJ33*01  TRAC*01    CAVKDSNYQLIW     CMV
#> 5   TRBV2*01 TRBD1*01 TRBJ1-2*01 TRBC1*01   CASNQGLNYGYTF     CMV
```

## Rhesus or mouse data

Use `species = "rhesus"` or `species = "mouse"`

``` r
new_tcrs <- convert_gene(
  tcrs,
  frm = "tenx",
  to = "imgt",
  species = "rhesus", # or 'mouse'
  verbose = FALSE
)
#> Warning in convert_gene(tcrs, frm = "tenx", to = "imgt", species = "rhesus", : These genes are not in IMGT for this species and will be replaced with NA:
#>  TRAV29/DV5, TRBV20/OR9-2, TRGJ1
new_tcrs
#>              barcode   v_gene     j_gene            cdr3
#> 1 AAACCTGAGACCACGA-1     <NA>  TRAJ12*01    CAVMDSSYKLIF
#> 2 AAACCTGAGACCACGA-1     <NA> TRBJ2-1*01 CASSGLAGGYNEQFF
#> 3 AAACCTGAGGCTCTTA-1 TRDV2*01   TRDJ3*01 CASSGVAGGTDTQYF
#> 4 AAACCTGAGGCTCTTA-1 TRGV9*01       <NA>    CAVKDSNYQLIW
```
