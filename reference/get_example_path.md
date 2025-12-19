# Get full path to an example file or directory

`get_example_path()` takes a file or folder name that is expected to be
located under the `TCRconvertR` `examples` directory and gets the full
path to that item.

## Usage

``` r
get_example_path(file_name)
```

## Arguments

- file_name:

  A string, the name of the example file or directory.

## Value

A string, the path to example file or directory.

## Examples

``` r
# Will probably be in a temp folder for the function example
get_example_path("tenx.csv")
#> [1] "/home/runner/work/_temp/Library/TCRconvertR/extdata/examples/tenx.csv"
```
