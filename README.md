
<!-- badges: start -->

[![codecov](https://codecov.io/gh/seshadrilab/tcrconvertr/graph/badge.svg?token=JVURVQO10D)](https://app.codecov.io/gh/seshadrilab/tcrconvertr)
[![R-CMD-check](https://github.com/seshadrilab/tcrconvertr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/seshadrilab/tcrconvertr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# Convert TCR gene names

`TCRconvertR` converts V, D, J, and/or C gene names between the 10X
Genomics, Adaptive Biotechnologies, and IMGT nomenclatures. It supports
alpha-beta and gamma-delta T cell receptors (TCRs) for human, mouse, and
rhesus macaque. Users can also define custom species, see:
`vignette("custom-species")`. A [Python
version](https://github.com/seshadrilab/tcrconvert) with command-line
support is also available.

## Background

TCR annotation tools use different gene naming conventions, making
cross-dataset searches difficult (e.g., identifying 10X-annotated TCRs
in Adaptive data). Manual conversion is complex and error-prone due to
inconsistencies in naming rules.

`TCRconvertR` automates this process efficiently and accurately. Our
approach is based on analyzing multiple 10X and Adaptive data sets to
capture their naming variations.

## Installation

Install from GitHub:

``` r
# install.packages("pak")
pak::pak("seshadrilab/tcrconvertr")
```

## Usage

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
#> 1 AAACCTGAGACCACGA-1    TRAV29DV5  TRAJ12    CAVMDSSYKLIF
#> 2 AAACCTGAGACCACGA-1 TRBV20/OR9-2 TRBJ2-1 CASSGLAGGYNEQFF
#> 3 AAACCTGAGGCTCTTA-1        TRDV2   TRDJ3 CASSGVAGGTDTQYF
#> 4 AAACCTGAGGCTCTTA-1        TRGV9   TRGJ1    CAVKDSNYQLIW
```

#### 2. Convert

``` r
new_tcrs <- convert_gene(tcrs, frm = "tenx", to = "adaptive")
#> Warning in convert_gene(tcrs, frm = "tenx", to = "adaptive"): Adaptive captures
#> only VDJ genes; C genes will be NA.
#> Converting from 10X. Using *01 as allele for all genes.
new_tcrs
#>              barcode             v_gene        j_gene            cdr3
#> 1 AAACCTGAGACCACGA-1      TCRAV29-01*01 TCRAJ12-01*01    CAVMDSSYKLIF
#> 2 AAACCTGAGACCACGA-1 TCRBV20-or09_02*01 TCRBJ02-01*01 CASSGLAGGYNEQFF
#> 3 AAACCTGAGGCTCTTA-1      TCRDV02-01*01 TCRDJ03-01*01 CASSGVAGGTDTQYF
#> 4 AAACCTGAGGCTCTTA-1      TCRGV09-01*01 TCRGJ01-01*01    CAVKDSNYQLIW
```

## Contributing

Contributions are welcome! To contribute, submit a pull request. See the
documentation for details.

## Issues

To report a bug or request a feature please open an
[issue](https://github.com/seshadrilab/tcrconvertr/issues).

## Contact

For other inquiries, contact Emma Bishop: emmab5 at uw dot edu.

## Acknowledgments

This project was supported by the Fred Hutchinson Cancer Center
Translational Data Science Integrated Research Center (TDS IRC) through
the 2024 Data Scientist Collaboration Grant. Special thanks to Scott
Chamberlain for development support and Shashidhar Ravishankar for gene
name curation.
