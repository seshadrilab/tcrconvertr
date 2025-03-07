test_that("can parse an imgt fasta", {
  fasta <- get_example_path("fasta_dir/test_trav.fa")
  expected_out <- c("TRAV1-1*01", "TRAV1-1*02", "TRAV14/DV4*01", "TRAV38-2/DV8*01", "TRAC*01")
  expect_equal(parse_imgt_fasta(fasta), expected_out)
})

test_that("can extract imgt genes", {
  fastadir <- get_example_path("fasta_dir")
  expected_out <- data.frame(imgt = c(
    "TRAC*01",
    "TRAV1-1*01",
    "TRAV1-1*02",
    "TRAV14/DV4*01",
    "TRAV38-2/DV8*01",
    "TRBV29/OR9-2*01",
    "TRBVA/OR9-2*01"
  ))
  expect_equal(extract_imgt_genes(fastadir), expected_out)
})

test_that("will add dash one", {
  gene_str1 <- "TRBV2*01"
  gene_str2 <- "TRBV1-01*01"
  expect_equal(add_dash_one(gene_str1), "TRBV2-01*01")
  expect_equal(add_dash_one(gene_str2), gene_str2)
})

test_that("will pad single digits", {
  gene_str1 <- "TCRBV1-2"
  gene_str2 <- "TCRBV11-2"
  expect_equal(pad_single_digit(gene_str1), "TCRBV01-2")
  expect_equal(pad_single_digit(gene_str2), gene_str2)
})


test_that("can build lookup tables from fastas", {
  skip_if_not_installed("mockery")

  # Create mock folder in temporary directory to write to
  mock_path <- file.path(tempdir(), "TCRconvertR_tmp")
  dir.create(mock_path, showWarnings = FALSE, recursive = TRUE)

  # Generate a temporary directory to put input files in
  fastadir <- file.path(mock_path, "fastas")
  dir.create(fastadir)
  # Copy test files into it
  example_dir <- get_example_path("fasta_dir")
  example_fastas <- list.files(example_dir, pattern = "\\.fa$", full.names = TRUE)
  file.copy(example_fastas, fastadir)

  # Stub `rappdirs::user_data_dir` within `build_lookup_from_fastas`
  mockery::stub(build_lookup_from_fastas, "rappdirs::user_data_dir", function(...) mock_path)

  # Create lookup tables
  suppressMessages(
    build_lookup_from_fastas(fastadir, species = "rabbit")
  )

  # Check adaptive lookup table
  adapt <- read.csv(file.path(mock_path, "rabbit/lookup_from_adaptive.csv"))
  expected_adapt <- data.frame(
    adaptive = c(
      "TCRAV01-01*01",
      "TCRAV01-01*02",
      "TCRAV14-01*01",
      "TCRAV38-02*01",
      "TCRBV29-or09_02*01",
      "TCRBVA-or09_02*01",
      "TCRAV01-01",
      "TCRAV14-01",
      "TCRAV38-02",
      "TCRBV29-or09_02",
      "TCRBVA-or09_02"
    ),
    adaptivev2 = c(
      "TCRAV01-01*01",
      "TCRAV01-01*02",
      "TCRAV14-01*01",
      "TCRAV38-02*01",
      "TCRBV29-or09_02*01",
      "TCRBVA-or09_02*01",
      "TCRAV01-01",
      "TCRAV14-01",
      "TCRAV38-02",
      "TCRBV29-or09_02",
      "TCRBVA-or09_02"
    ),
    imgt = c(
      "TRAV1-1*01",
      "TRAV1-1*02",
      "TRAV14/DV4*01",
      "TRAV38-2/DV8*01",
      "TRBV29/OR9-2*01",
      "TRBVA/OR9-2*01",
      "TRAV1-1*01",
      "TRAV14/DV4*01",
      "TRAV38-2/DV8*01",
      "TRBV29/OR9-2*01",
      "TRBVA/OR9-2*01"
    ),
    tenx = c(
      "TRAV1-1",
      "TRAV1-1",
      "TRAV14DV4",
      "TRAV38-2DV8",
      "TRBV29/OR9-2",
      "TRBVA/OR9-2",
      "TRAV1-1",
      "TRAV14DV4",
      "TRAV38-2DV8",
      "TRBV29/OR9-2",
      "TRBVA/OR9-2"
    )
  )
  expect_equal(adapt, expected_adapt)

  # Check 10X lookup table
  tenx <- read.csv(file.path(mock_path, "rabbit/lookup_from_tenx.csv"))
  expected_tenx <- data.frame(
    tenx = c(
      "TRAC",
      "TRAV1-1",
      "TRAV14DV4",
      "TRAV38-2DV8",
      "TRBV29/OR9-2",
      "TRBVA/OR9-2"
    ),
    imgt = c(
      "TRAC*01",
      "TRAV1-1*01",
      "TRAV14/DV4*01",
      "TRAV38-2/DV8*01",
      "TRBV29/OR9-2*01",
      "TRBVA/OR9-2*01"
    ),
    adaptive = c(
      "NoData",
      "TCRAV01-01*01",
      "TCRAV14-01*01",
      "TCRAV38-02*01",
      "TCRBV29-or09_02*01",
      "TCRBVA-or09_02*01"
    ),
    adaptivev2 = c(
      "NoData",
      "TCRAV01-01*01",
      "TCRAV14-01*01",
      "TCRAV38-02*01",
      "TCRBV29-or09_02*01",
      "TCRBVA-or09_02*01"
    )
  )
  expect_equal(tenx, expected_tenx)

  # Check regular lookup table
  lookup <- read.csv(file.path(mock_path, "rabbit/lookup.csv"))
  expected_lookup <- data.frame(
    imgt = c(
      "TRAC*01",
      "TRAV1-1*01",
      "TRAV1-1*02",
      "TRAV14/DV4*01",
      "TRAV38-2/DV8*01",
      "TRBV29/OR9-2*01",
      "TRBVA/OR9-2*01"
    ),
    tenx = c(
      "TRAC",
      "TRAV1-1",
      "TRAV1-1",
      "TRAV14DV4",
      "TRAV38-2DV8",
      "TRBV29/OR9-2",
      "TRBVA/OR9-2"
    ),
    adaptive = c(
      "NoData",
      "TCRAV01-01*01",
      "TCRAV01-01*02",
      "TCRAV14-01*01",
      "TCRAV38-02*01",
      "TCRBV29-or09_02*01",
      "TCRBVA-or09_02*01"
    ),
    adaptivev2 = c(
      "NoData",
      "TCRAV01-01*01",
      "TCRAV01-01*02",
      "TCRAV14-01*01",
      "TCRAV38-02*01",
      "TCRBV29-or09_02*01",
      "TCRBVA-or09_02*01"
    )
  )
  expect_equal(lookup, expected_lookup)

  # Delete temp directory
  unlink(mock_path, recursive = TRUE)
})

test_that("will reject invalid species name", {
  skip_if_not_installed("mockery")

  # Create mock folder in temporary directory to write to
  mock_path <- file.path(tempdir(), "TCRconvertR_tmp")
  dir.create(mock_path, showWarnings = FALSE, recursive = TRUE)

  # Generate a temporary directory to put input files in
  fastadir <- file.path(mock_path, "fastas")
  dir.create(fastadir)
  # Copy test files into it
  example_dir <- get_example_path("fasta_dir")
  example_fastas <- list.files(example_dir, pattern = "\\.fa$", full.names = TRUE)
  file.copy(example_fastas, fastadir)

  # Stub `rappdirs::user_data_dir` within `build_lookup_from_fastas`
  mockery::stub(build_lookup_from_fastas, "rappdirs::user_data_dir", function(...) mock_path)

  # Attempt to create lookup tables with invalid species folder names
  expect_error(build_lookup_from_fastas(fastadir, species = "rabbit\\"))
  expect_error(build_lookup_from_fastas(fastadir, species = "rabbit/"))
  expect_error(build_lookup_from_fastas(fastadir, species = "rabbit:"))
  expect_error(build_lookup_from_fastas(fastadir, species = "rabbit*"))
  expect_error(build_lookup_from_fastas(fastadir, species = "rabbit?"))
  expect_error(build_lookup_from_fastas(fastadir, species = "rabbit<"))
  expect_error(build_lookup_from_fastas(fastadir, species = "rabbit>"))
  expect_error(build_lookup_from_fastas(fastadir, species = "rabbit|"))
  expect_error(build_lookup_from_fastas(fastadir, species = "rabbit~"))
  expect_error(build_lookup_from_fastas(fastadir, species = "rabbit`"))
  expect_error(build_lookup_from_fastas(fastadir, species = "rabbit\n"))
  expect_error(build_lookup_from_fastas(fastadir, species = "rabbi\t"))

  # Delete temp directories
  unlink(mock_path, recursive = TRUE)
  unlink(fastadir, recursive = TRUE)
})

test_that("will save lookup tables", {
  save_dir <- file.path(tempdir(), "TCRconvertR_tmp")
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  dat <- read.csv(get_example_path("fasta_dir/lookup.csv"))

  save_lookup(dat, save_dir, "newlookup.csv")

  # Check that the file exists
  expect_true(file.exists(file.path(save_dir, "newlookup.csv")))

  # Delete temp directory
  unlink(save_dir, recursive = TRUE)
})
