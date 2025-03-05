
<!-- badges: start -->

[![codecov](https://codecov.io/gh/seshadrilab/tcrconvertr/graph/badge.svg?token=JVURVQO10D)](https://codecov.io/gh/seshadrilab/tcrconvertr)
[![R-CMD-check](https://github.com/seshadrilab/tcrconvertr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/seshadrilab/tcrconvertr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Convert T cell receptor (TCR) gene names

`TCRconvertR` converts V, D, J, and/or C gene names between the 10X
Genomics, Adaptive Biotechnologies, and IMGT nomenclatures. It supports
alpha-beta and gamma-delta TCRs for human, mouse, and rhesus macaque.
Users can also define custom species (see docs). A Python version with
command-line support is also available.

## Background

TCR annotation tools use different gene naming conventions, making
cross-dataset searches difficult (e.g., identifying 10X-annotated TCRs
in Adaptive data). Manual conversion is complex and error-prone due to
inconsistencies in naming rules.

`TCRconvert` automates this process efficiently and accurately. Our
approach is based on analyzing multiple 10X and Adaptive data sets to
capture their naming variations. Additionally, we provide a guide (see
docs) for generating the necessary gene name columns in Adaptive data
when missing.

## Installation

Install from [GitHub](https://github.com/):

``` r
# install.packages("pak")
pak::pak("seshadrilab/tcrconvertr")
```

## Usage

**1. Load your TCRs into a data frame**

Example input files:

- **10X**: `filtered_contig_annotations.csv`
- **Adaptive**: `Sample_TCRB.tsv`
- **IMGT**: Output from `MiXCR` or other tools

``` r
library(TCRconvertR)

tcr_file <- get_example_path("tenx.csv")
tcrs <- read.csv(tcr_file)[c("barcode", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3")]
tcrs
#>              barcode  v_gene d_gene  j_gene c_gene            cdr3
#> 1 AAACCTGAGACCACGA-1 TRAV1-2  TRBD1  TRAJ12   TRAC    CAVMDSSYKLIF
#> 2 AAACCTGAGACCACGA-1 TRBV6-1  TRBD2 TRBJ2-1  TRBC2 CASSGLAGGYNEQFF
#> 3 AAACCTGAGGCTCTTA-1 TRBV6-4  TRBD2 TRBJ2-3  TRBC2 CASSGVAGGTDTQYF
#> 4 AAACCTGAGGCTCTTA-1 TRAV1-2  TRBD1  TRAJ33   TRAC    CAVKDSNYQLIW
#> 5 AAACCTGAGTGAACGC-1   TRBV2  TRBD1 TRBJ1-2  TRBC1   CASNQGLNYGYTF
```

**2. Convert**

``` r
new_tcrs <- convert_gene(tcrs, frm = "tenx", to = "adaptive")
#> Warning in convert_gene(tcrs, frm = "tenx", to = "adaptive"): Adaptive only
#> captures VDJ genes, any C genes will become NA.
#> Converting from 10X which lacks allele info. Choosing *01 as allele for all genes.
new_tcrs
#>              barcode        v_gene        d_gene        j_gene c_gene
#> 1 AAACCTGAGACCACGA-1 TCRAV01-02*01 TCRBD01-01*01 TCRAJ12-01*01   <NA>
#> 2 AAACCTGAGACCACGA-1 TCRBV06-01*01 TCRBD02-01*01 TCRBJ02-01*01   <NA>
#> 3 AAACCTGAGGCTCTTA-1 TCRBV06-04*01 TCRBD02-01*01 TCRBJ02-03*01   <NA>
#> 4 AAACCTGAGGCTCTTA-1 TCRAV01-02*01 TCRBD01-01*01 TCRAJ33-01*01   <NA>
#> 5 AAACCTGAGTGAACGC-1 TCRBV02-01*01 TCRBD01-01*01 TCRBJ01-02*01   <NA>
#>              cdr3
#> 1    CAVMDSSYKLIF
#> 2 CASSGLAGGYNEQFF
#> 3 CASSGVAGGTDTQYF
#> 4    CAVKDSNYQLIW
#> 5   CASNQGLNYGYTF
```

## Contributing

Feedback is welcome! To contribute, submit a pull request.

## Issues

For problems or questions, file an issue on GitHub.

## Contact

For other inquiries, contact Emma Bishop: emmab5 at uw dot edu.

## Acknowledgments

This project was supported by the Fred Hutchinson Cancer Center
Translational Data Science Integrated Research Center (TDS IRC) through
the 2024 Data Scientist Collaboration Grant. Special thanks to Scott
Chamberlain for development support and Shashidhar Ravishankar for gene
name curation.
