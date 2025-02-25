
# TCRconvertR

<!-- badges: start -->

[![R-CMD-check](https://github.com/seshadrilab/tcrconvertr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/seshadrilab/tcrconvertr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

> **Warning**: This project is in **beta stage**. It is under active
> development and may be unstable.

**Rename T-cell receptor genes between 10X, Adaptive, and IMGT formats**

TCRconvertR takes T-cell receptor (TCR) data containing V, D, J, and/or
C genes from 10X, Adaptive, or other sequencing platforms and renames
them from any of these formats to any other one:

- **10X**: TRAV1-2
- **Adaptive**: TCRAV01-02\*01
- **IMGT**: TRAV1-2\*01

TCRconvertR works with human, mouse, and rhesus macaque data
out-of-the-box, but users can also add their own species.

TCRconvertR helps researchers unify TCR datasets by converting them to a
standard naming convention. It is fast, reliable, and prevents errors
from manual conversions. Unlike other tools that require custom objects,
TCRconvertR works directly with DataFrames and CSV/TSV files.

# Installation

You can install TCRconvertR from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("seshadrilab/tcrconvertr")
```

# Usage

**Load some 10X data**

``` r
library(TCRconvertR)

tcr_file <- "/home/emmabishop/workspace/tcrconvertr/inst/extdata/examples/tenx.csv"

tcrs <- read.csv(tcr_file)
tcrs
#>              barcode is_cell                   contig_id high_confidence length
#> 1 AAACCTGAGACCACGA-1    TRUE AAACCTGAGACCACGA-1_contig_1            TRUE    521
#> 2 AAACCTGAGACCACGA-1    TRUE AAACCTGAGACCACGA-1_contig_3            TRUE    584
#> 3 AAACCTGAGGCTCTTA-1    TRUE AAACCTGAGGCTCTTA-1_contig_1            TRUE    551
#> 4 AAACCTGAGGCTCTTA-1    TRUE AAACCTGAGGCTCTTA-1_contig_2            TRUE    518
#> 5 AAACCTGAGTGAACGC-1    TRUE AAACCTGAGTGAACGC-1_contig_1            TRUE    674
#>   chain  v_gene d_gene  j_gene c_gene full_length productive            cdr3
#> 1   TRA TRAV1-2  TRBD1  TRAJ12   TRAC        TRUE       TRUE    CAVMDSSYKLIF
#> 2   TRB TRBV6-1  TRBD2 TRBJ2-1  TRBC2        TRUE       TRUE CASSGLAGGYNEQFF
#> 3   TRB TRBV6-4  TRBD2 TRBJ2-3  TRBC2        TRUE       TRUE CASSGVAGGTDTQYF
#> 4   TRA TRAV1-2  TRBD1  TRAJ33   TRAC        TRUE       TRUE    CAVKDSNYQLIW
#> 5   TRB   TRBV2  TRBD1 TRBJ1-2  TRBC1        TRUE       TRUE   CASNQGLNYGYTF
#>                                         cdr3_nt reads umis raw_clonotype_id
#> 1          TGTGCTGTGATGGATAGCAGCTATAAATTGATCTTC  1569    2      clonotype16
#> 2 TGTGCCAGCAGTGGACTAGCGGGGGGCTACAATGAGCAGTTCTTC  5238    7      clonotype16
#> 3 TGTGCCAGCAGTGGGGTAGCGGGAGGCACAGATACGCAGTATTTT  3846    4      clonotype26
#> 4          TGTGCTGTGAAGGATAGCAACTATCAGTTAATCTGG  2019    2      clonotype26
#> 5       TGTGCCAGCAATCAGGGCCTTAACTATGGCTACACCTTC  3002    6      clonotype81
#>          raw_consensus_id
#> 1 clonotype16_consensus_1
#> 2 clonotype16_consensus_2
#> 3 clonotype26_consensus_2
#> 4 clonotype26_consensus_1
#> 5 clonotype81_consensus_2
```

**Convert gene names from the 10X format to the Adaptive format**

``` r
new_tcrs <- convert_gene(tcrs, frm = "tenx", to = "adaptive")
#> Warning in convert_gene(tcrs, frm = "tenx", to = "adaptive"): Adaptive only
#> captures VDJ genes, any C genes will become NA.
#> Converting from 10X which lacks allele info. Choosing *01 as allele for all genes.
new_tcrs
#>              barcode is_cell                   contig_id high_confidence length
#> 1 AAACCTGAGACCACGA-1    TRUE AAACCTGAGACCACGA-1_contig_1            TRUE    521
#> 2 AAACCTGAGACCACGA-1    TRUE AAACCTGAGACCACGA-1_contig_3            TRUE    584
#> 3 AAACCTGAGGCTCTTA-1    TRUE AAACCTGAGGCTCTTA-1_contig_1            TRUE    551
#> 4 AAACCTGAGGCTCTTA-1    TRUE AAACCTGAGGCTCTTA-1_contig_2            TRUE    518
#> 5 AAACCTGAGTGAACGC-1    TRUE AAACCTGAGTGAACGC-1_contig_1            TRUE    674
#>   chain        v_gene        d_gene        j_gene c_gene full_length productive
#> 1   TRA TCRAV01-02*01 TCRBD01-01*01 TCRAJ12-01*01   <NA>        TRUE       TRUE
#> 2   TRB TCRBV06-01*01 TCRBD02-01*01 TCRBJ02-01*01   <NA>        TRUE       TRUE
#> 3   TRB TCRBV06-04*01 TCRBD02-01*01 TCRBJ02-03*01   <NA>        TRUE       TRUE
#> 4   TRA TCRAV01-02*01 TCRBD01-01*01 TCRAJ33-01*01   <NA>        TRUE       TRUE
#> 5   TRB TCRBV02-01*01 TCRBD01-01*01 TCRBJ01-02*01   <NA>        TRUE       TRUE
#>              cdr3                                       cdr3_nt reads umis
#> 1    CAVMDSSYKLIF          TGTGCTGTGATGGATAGCAGCTATAAATTGATCTTC  1569    2
#> 2 CASSGLAGGYNEQFF TGTGCCAGCAGTGGACTAGCGGGGGGCTACAATGAGCAGTTCTTC  5238    7
#> 3 CASSGVAGGTDTQYF TGTGCCAGCAGTGGGGTAGCGGGAGGCACAGATACGCAGTATTTT  3846    4
#> 4    CAVKDSNYQLIW          TGTGCTGTGAAGGATAGCAACTATCAGTTAATCTGG  2019    2
#> 5   CASNQGLNYGYTF       TGTGCCAGCAATCAGGGCCTTAACTATGGCTACACCTTC  3002    6
#>   raw_clonotype_id        raw_consensus_id
#> 1      clonotype16 clonotype16_consensus_1
#> 2      clonotype16 clonotype16_consensus_2
#> 3      clonotype26 clonotype26_consensus_2
#> 4      clonotype26 clonotype26_consensus_1
#> 5      clonotype81 clonotype81_consensus_2
```

# Contributing

I welcome feedback! If you would like to resolve an issue or add
improvements please submit a pull request.

# Issues

If you run into problems or need help running TCRconvertR please file an
issue on GitHub.

# Contact

For other questions please contact Emma Bishop: `emmab5` at `uw` dot
`edu`

# Acknowledgments

This project was created with support from the Fred Hutchinson Cancer
Center Translational Data Science Integrated Research Center (TDS IRC)
through the 2024 Data Scientist Collaboration Grant.
