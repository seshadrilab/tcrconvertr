get_example_path <- function(file_name) {
  # Get full path to a given example file or directory
  #
  # Args:
  #   file_name: Name of the example file or directory (string)
  #
  # Returns:
  #   Full path to the example file or directory (string)
  
  # Assuming the "tcrconvert" package has an "examples" folder
  # Use system.file to locate the directory
  example_dir <- system.file("examples", package = "tcrconvert")
  
  # Construct the full path
  file_path <- file.path(example_dir, file_name)
  
  return(file_path)
}

# Example usage:
# Assuming the "tcrconvert" package contains an "examples" folder
# print(get_example_path("tenx.csv"))
