# Save a lookup table to a CSV file

`save_lookup()` saves a data frame as a CSV file (without row names) in
the specified directory.

## Usage

``` r
save_lookup(df, savedir, name)
```

## Arguments

- df:

  A data frame containing the lookup table data.

- savedir:

  A string, the path to the save directory.

- name:

  A string, the file name (should end in `.csv`).

## Value

Nothing

## Examples

``` r
# Create a temp save directory and load an example
save_dir <- file.path(tempdir(), "TCRconvertR_tmp")
dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
dat <- read.csv(get_example_path("fasta_dir/lookup.csv"))

save_lookup(dat, save_dir, "newlookup.csv")

# Clean up temporary folder
unlink(save_dir, recursive = TRUE)
```
