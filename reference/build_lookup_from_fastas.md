# Create lookup tables

`build_lookup_from_fastas()` processes IMGT reference FASTA files in a
given folder to generate lookup tables used for making gene name
conversions. It extracts all gene names and transforms them into 10X and
Adaptive formats following predefined conversion rules. The resulting
files are created:

- `lookup.csv`: IMGT gene names and their 10X and Adaptive equivalents.

- `lookup_from_tenx.csv`: Gene names aggregated by their 10X
  identifiers, with one representative allele (`*01`) for each.

- `lookup_from_adaptive.csv`: Adaptive gene names, with or without
  alleles and gene designations, and their IMGT and 10X equivalents.

The files are stored in a given subfolder (`species`) within the
appropriate application folder via `rappdirs`. For example:

- MacOS: `~/Library/Application Support/<AppName>`

- Windows:
  `C:\Documents and Settings\<User>\Application Data\Local Settings\<AppAuthor>\<AppName>`

- Linux: `~/.local/share/<AppName>`

If a folder named `species` already exists in that location, it will be
replaced.

## Usage

``` r
build_lookup_from_fastas(data_dir, species)
```

## Arguments

- data_dir:

  A string, the directory containing FASTA files.

- species:

  A string, the name of species that will be used when running
  TCRconvert with these lookup tables.

## Value

A string, path to new lookup directory

## Details

Key transformations from IMGT:

- **10X:**

  - Remove allele information (e.g., `*01`) and modify `/DV`
    occurrences.

- **Adaptive:**

  - Apply renaming rules, such as adding gene-level designations and
    zero-padding single-digit numbers.

  - Convert constant genes to `"NoData"` (Adaptive only captures VDJ)
    which become `NA` after the merge in
    [`convert_gene()`](https://seshadrilab.github.io/tcrconvertr/reference/convert_gene.md).

## Examples

``` r
# For the example, create and use a temporary folder
fastadir <- file.path(tempdir(), "TCRconvertR_tmp")
dir.create(fastadir, showWarnings = FALSE, recursive = TRUE)
trav <- get_example_path("fasta_dir/test_trav.fa")
trbv <- get_example_path("fasta_dir/test_trbv.fa")
file.copy(c(trav, trbv), fastadir)
#> [1] TRUE TRUE

# Build lookup tables
build_lookup_from_fastas(fastadir, "rabbit")
#> Writing lookup tables to: ~/.local/share/TCRconvertR/rabbit
#> [1] "~/.local/share/TCRconvertR/rabbit"

# Clean up temporary folder
unlink(fastadir, recursive = TRUE)
```
