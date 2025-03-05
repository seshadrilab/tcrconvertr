#' Get full path to an example file or directory
#'
#' `get_example_path()` takes a file or folder name that is expected to be
#' located under the `TCRconvertR` `examples` directory and gets the full path
#' to that item.
#'
#' @param file_name A string, the name of the example file or directory.
#'
#' @return A string, the path to example file or directory.
#' @export
#' @examples
#' get_example_path("tenx.csv")
get_example_path <- function(file_name) {
  # Use system.file to locate the directory
  example_dir <- system.file("extdata/examples", package = "TCRconvertR")

  # Construct the full path
  file_path <- file.path(example_dir, file_name)

  return(file_path)
}
