test_that("can convert genes", {
  # First, define all the dataframes we need
  imgt_df <- data.frame(
    v_gene = c("TRAV12-1*01", "TRBV15*01"),
    d_gene = c(NA, "TRBD1*01"),
    j_gene = c("TRAJ16*01", "TRBJ2-5*01"),
    c_gene = c("TRAC*01", "TRBC2*01"),
    cdr3 = c("CAVLIF", "CASSGF")
  )

  tenx_df <- data.frame(
    v_gene = c("TRAV12-1", "TRBV15"),
    d_gene = c(NA, "TRBD1"),
    j_gene = c("TRAJ16", "TRBJ2-5"),
    c_gene = c("TRAC", "TRBC2"),
    cdr3 = c("CAVLIF", "CASSGF")
  )

  adapt_df <- data.frame(
    v_resolved = c("TCRAV12-01*01", "TCRBV15-01*01"),
    d_resolved = c(NA, "TCRBD01-01*01"),
    j_resolved = c("TCRAJ16-01*01", "TCRBJ02-05*01"),
    cdr3_amino_acid = c("CAVLIF", "CASSGF")
  )

  adapt_v2_df <- data.frame(
    vMaxResolved = c("TCRAV12-01*01", "TCRBV15-01*01"),
    dMaxResolved = c(NA, "TCRBD01-01*01"),
    jMaxResolved = c("TCRAJ16-01*01", "TCRBJ02-05*01"),
    aminoAcid = c("CAVLIF", "CASSGF")
  )

  custom_df <- imgt_df
  colnames(custom_df) <- c("myV", "myD", "myJ", "myC", "myCDR3")

  custom_vj_tenx_df <- data.frame(
    myV = c("TRAV12-1", "TRBV15"),
    myD = c(NA, "TRBD1*01"),
    myJ = c("TRAJ16", "TRBJ2-5"),
    myC = c("TRAC*01", "TRBC2*01"),
    myCDR3 = c("CAVLIF", "CASSGF")
  )

  tenx_to_adapt_df <- adapt_df
  colnames(tenx_to_adapt_df) <- c("v_gene", "d_gene", "j_gene", "cdr3")
  # Insert a 'c_gene' column before 'cdr3'
  tenx_to_adapt_df <- cbind(
    tenx_to_adapt_df[, 1:3],
    c_gene = NA,
    tenx_to_adapt_df["cdr3"]
  )
  tenx_to_adapt_df$c_gene <- as.character(tenx_to_adapt_df$c_gene)

  adapt_to_tenx_df <- tenx_df
  adapt_to_tenx_df$c_gene <- NULL
  colnames(adapt_to_tenx_df) <- c("v_resolved", "d_resolved", "j_resolved", "cdr3_amino_acid")

  adapt_to_imgt_df <- imgt_df
  adapt_to_imgt_df$c_gene <- NULL
  colnames(adapt_to_imgt_df) <- c("v_resolved", "d_resolved", "j_resolved", "cdr3_amino_acid")

  adaptv2_to_tenx_df <- tenx_df
  adaptv2_to_tenx_df$c_gene <- NULL
  colnames(adaptv2_to_tenx_df) <- c("vMaxResolved", "dMaxResolved", "jMaxResolved", "aminoAcid")

  adaptv2_to_imgt_df <- imgt_df
  adaptv2_to_imgt_df$c_gene <- NULL
  colnames(adaptv2_to_imgt_df) <- c("vMaxResolved", "dMaxResolved", "jMaxResolved", "aminoAcid")

  custom_to_tenx_df <- tenx_df
  colnames(custom_to_tenx_df) <- c("myV", "myD", "myJ", "myC", "myCDR3")

  adapt_no_allele_df <- data.frame(
    v_resolved = c("TCRAV12-01", "TCRBV15-01*01"),
    d_resolved = c(NA, "TCRBD01-01"),
    j_resolved = c("TCRAJ16-01*01", "TCRBJ02-05"),
    cdr3_amino_acid = c("CAVLIF", "CASSGF")
  )
  # 10X <-> Adaptive
  suppressMessages(suppressWarnings({
    expect_equal(convert_gene(tenx_df, "tenx", "adaptive"), tenx_to_adapt_df)
    expect_equal(convert_gene(tenx_df, "tenx", "adaptivev2"), tenx_to_adapt_df)
    expect_equal(convert_gene(adapt_df, "adaptive", "tenx"), adapt_to_tenx_df)
    expect_equal(convert_gene(adapt_v2_df, "adaptivev2", "tenx"), adaptv2_to_tenx_df)
    # 10X <-> IMGT
    expect_equal(convert_gene(tenx_df, "tenx", "imgt"), imgt_df)
    expect_equal(convert_gene(imgt_df, "imgt", "tenx"), tenx_df)
    # IMGT <-> Adaptive
    expect_equal(convert_gene(imgt_df, "imgt", "adaptive"), tenx_to_adapt_df)
    expect_equal(convert_gene(imgt_df, "imgt", "adaptivev2"), tenx_to_adapt_df)
    expect_equal(convert_gene(adapt_df, "adaptive", "imgt"), adapt_to_imgt_df)
    expect_equal(convert_gene(adapt_v2_df, "adaptivev2", "imgt"), adaptv2_to_imgt_df)
    # Custom column names
    expect_equal(convert_gene(custom_df, "imgt", "tenx",
      frm_cols = c("myV", "myD", "myJ", "myC")
    ), custom_to_tenx_df)
    # MOUSE
    expect_equal(convert_gene(tenx_df, "tenx", "adaptive", species = "mouse"), tenx_to_adapt_df)
    expect_equal(convert_gene(tenx_df, "tenx", "adaptivev2", species = "mouse"), tenx_to_adapt_df)
    expect_equal(convert_gene(adapt_df, "adaptive", "tenx", species = "mouse"), adapt_to_tenx_df)
    expect_equal(convert_gene(adapt_v2_df, "adaptivev2", "tenx", species = "mouse"), adaptv2_to_tenx_df)
    expect_equal(convert_gene(tenx_df, "tenx", "imgt", species = "mouse"), imgt_df)
    expect_equal(convert_gene(imgt_df, "imgt", "tenx", species = "mouse"), tenx_df)
    expect_equal(convert_gene(imgt_df, "imgt", "adaptive", species = "mouse"), tenx_to_adapt_df)
    expect_equal(convert_gene(imgt_df, "imgt", "adaptivev2", species = "mouse"), tenx_to_adapt_df)
    expect_equal(convert_gene(adapt_df, "adaptive", "imgt", species = "mouse"), adapt_to_imgt_df)
    expect_equal(convert_gene(adapt_v2_df, "adaptivev2", "imgt", species = "mouse"), adaptv2_to_imgt_df)
    expect_equal(convert_gene(custom_df, "imgt", "tenx",
      species = "mouse",
      frm_cols = c("myV", "myD", "myJ", "myC")
    ), custom_to_tenx_df)
    # RHESUS MACAQUE
    expect_equal(convert_gene(tenx_df, "tenx", "adaptive", species = "rhesus"), tenx_to_adapt_df)
    expect_equal(convert_gene(tenx_df, "tenx", "adaptivev2", species = "rhesus"), tenx_to_adapt_df)
    expect_equal(convert_gene(adapt_df, "adaptive", "tenx", species = "rhesus"), adapt_to_tenx_df)
    expect_equal(convert_gene(adapt_v2_df, "adaptivev2", "tenx", species = "rhesus"), adaptv2_to_tenx_df)
    expect_equal(convert_gene(tenx_df, "tenx", "imgt", species = "rhesus"), imgt_df)
    expect_equal(convert_gene(imgt_df, "imgt", "tenx", species = "rhesus"), tenx_df)
    expect_equal(convert_gene(imgt_df, "imgt", "adaptive", species = "rhesus"), tenx_to_adapt_df)
    expect_equal(convert_gene(imgt_df, "imgt", "adaptivev2", species = "rhesus"), tenx_to_adapt_df)
    expect_equal(convert_gene(adapt_df, "adaptive", "imgt", species = "rhesus"), adapt_to_imgt_df)
    expect_equal(convert_gene(adapt_v2_df, "adaptivev2", "imgt", species = "rhesus"), adaptv2_to_imgt_df)
    expect_equal(convert_gene(custom_df, "imgt", "tenx",
      species = "rhesus",
      frm_cols = c("myV", "myD", "myJ", "myC")
    ), custom_to_tenx_df)
    # Some Adaptive genes without allele
    expect_equal(convert_gene(adapt_no_allele_df, "adaptive", "imgt"), adapt_to_imgt_df)
    # Confirm won't convert non-VDJC gene column to NAs
    expect_equal(convert_gene(custom_df, "imgt", "tenx",
      frm_cols = c("myV", "myJ", "myCDR3")
    ), custom_vj_tenx_df)
  }))
})


test_that("bad input raises error", {
  tenx_df <- data.frame(
    v_gene = c("TRAV12-1", "TRBV15"),
    d_gene = c(NA, "TRBD1"),
    j_gene = c("TRAJ16", "TRBJ2-5"),
    c_gene = c("TRAC", "TRBC2"),
    cdr3 = c("CAVLIF", "CASSGF")
  )

  # Same input and output format
  expect_error(convert_gene(tenx_df, "tenx", "tenx"), '"frm" and "to" formats should be different.')
  # Empty input
  expect_error(convert_gene(data.frame(), "tenx", "imgt"), "Input data is empty.")
})


test_that("chooses correct lookup", {
  lookup_tenx <- system.file("extdata/human", "lookup_from_tenx.csv", package = "TCRconvertR")
  lookup_adapt <- system.file("extdata/human", "lookup_from_adaptive.csv", package = "TCRconvertR")
  lookup_imgt <- system.file("extdata/human", "lookup.csv", package = "TCRconvertR")

  # Species we don't have lookups for
  expect_error(
    choose_lookup("tenx", "imgt", "non-existent-species", verbose = FALSE),
    "Lookup table not found, please ensure reference files are available."
  )
  # From different 'frm' formats
  expect_equal(choose_lookup("tenx", "imgt", verbose = FALSE), lookup_tenx)
  expect_equal(choose_lookup("adaptivev2", "imgt", verbose = FALSE), lookup_adapt)
  expect_equal(choose_lookup("imgt", "tenx", verbose = FALSE), lookup_imgt)
})


test_that("chooses correct frm_cols", {
  col_ref <- list(
    adaptive = c("v_resolved", "d_resolved", "j_resolved"),
    adaptivev2 = c("vMaxResolved", "dMaxResolved", "jMaxResolved"),
    imgt = c("v_gene", "d_gene", "j_gene", "c_gene"),
    tenx = c("v_gene", "d_gene", "j_gene", "c_gene")
  )

  adapt_df <- data.frame(
    v_resolved = c("TCRAV12-01*01", "TCRBV15-01*01"),
    d_resolved = c(NA, "TCRBD01-01*01"),
    j_resolved = c("TCRAJ16-01*01", "TCRBJ02-05*01"),
    cdr3_amino_acid = c("CAVLIF", "CASSGF")
  )

  adapt_v2_df <- data.frame(
    vMaxResolved = c("TCRAV12-01*01", "TCRBV15-01*01"),
    dMaxResolved = c(NA, "TCRBD01-01*01"),
    jMaxResolved = c("TCRAJ16-01*01", "TCRBJ02-05*01"),
    aminoAcid = c("CAVLIF", "CASSGF")
  )

  imgt_df <- data.frame(
    v_gene = c("TRAV12-1*01", "TRBV15*01"),
    d_gene = c(NA, "TRBD1*01"),
    j_gene = c("TRAJ16*0", "TRBJ2-5*01"),
    c_gene = c("TRAC*01", "TRBC2*01"),
    cdr3 = c("CAVLIF", "CASSGF")
  )

  tenx_df <- data.frame(
    v_gene = c("TRAV12-1", "TRBV15"),
    d_gene = c(NA, "TRBD1"),
    j_gene = c("TRAJ16", "TRBJ2-5"),
    c_gene = c("TRAC", "TRBC2"),
    cdr3 = c("CAVLIF", "CASSGF")
  )

  custom_df <- data.frame(
    myV = c("TRAV12-1*01", "TRBV15*01"),
    myD = c(NA, "TRBD1*01"),
    myJ = c("TRAJ16*0", "TRBJ2-5*01"),
    myC = c("TRAC*01", "TRBC2*01"),
    myCDR3 = c("CAVLIF", "CASSGF")
  )

  suppressWarnings({
    expect_equal(which_frm_cols(adapt_df, "adaptive", verbose = FALSE), col_ref$adaptive)
    expect_equal(which_frm_cols(adapt_v2_df, "adaptivev2", verbose = FALSE), col_ref$adaptivev2)
    expect_equal(which_frm_cols(imgt_df, "imgt", verbose = FALSE), col_ref$imgt)
    expect_equal(which_frm_cols(tenx_df, "tenx", verbose = FALSE), col_ref$tenx)
    # Custom columns
    custom_col <- c("myV", "myD", "myJ", "myC")
    expect_equal(which_frm_cols(custom_df, "tenx", frm_cols = custom_col, verbose = FALSE), custom_col)
    # Non-existent column
    expect_error(
      which_frm_cols(tenx_df, "tenx", frm_cols = c("v_gene", "j_gene", "x_gene"), verbose = FALSE),
      "These columns are not in the input dataframe: x_gene"
    )
  })
})


test_that("choose_lookup verbose flag works", {
  # Capture messages when verbose = TRUE
  expect_message(choose_lookup("tenx", "adaptive"),
    fixed = TRUE,
    "Converting from 10X. Using *01 as allele for all genes."
  )
  expect_message(choose_lookup("adaptive", "imgt"),
    fixed = TRUE,
    "Converting from Adaptive to IMGT. Using *01 for genes lacking alleles."
  )
  # Ensure no messages when verbose = FALSE
  expect_silent(choose_lookup("tenx", "adaptive", verbose = FALSE))
  expect_silent(choose_lookup("adaptive", "imgt", verbose = FALSE))
})


test_that("which_frm_cols verbose flag works", {
  custom_df <- data.frame(
    myV = c("TRAV12-1*01", "TRBV15*01"),
    myD = c(NA, "TRBD1*01"),
    myJ = c("TRAJ16*0", "TRBJ2-5*01"),
    myC = c("TRAC*01", "TRBC2*01"),
    myCDR3 = c("CAVLIF", "CASSGF")
  )
  # Custom columns
  expect_message(
    which_frm_cols(custom_df, "imgt", frm_cols = c("myV", "myJ")),
    "Using custom column names: myV, myJ"
  )
  # Custom columns with verbose = FALSE
  expect_silent(
    which_frm_cols(custom_df, "imgt", frm_cols = c("myV", "myJ"), verbose = FALSE)
  )
  # IMGT format without column names, still expect warnings with verbose = FALSE
  expect_warning(
    which_frm_cols(custom_df, "imgt", verbose = FALSE),
    "No column names for IMGT data. Using 10X columns: v_gene, d_gene, j_gene, c_gene"
  )
})



test_that("convert_gene verbose flag works", {
  # Note, the flag mostly affects downstream functions called by convert_gene.
  tenx_df_bad <- data.frame(
    v_gene = c("TRAV12-1", "TRBV15"),
    d_gene = c(NA, "TRBD1"),
    j_gene = c("TRAJ16", "BAD_J_GENE"),
    c_gene = c("TRAC", "TRBC2"),
    cdr3 = c("CAVLIF", "CASSGF")
  )

  captured_warnings <- character()

  # Capture warnings without printing
  result <- withCallingHandlers(
    convert_gene(tenx_df_bad, "tenx", "adaptive", verbose = FALSE),
    warning = function(w) {
      # Add each one as it comes up to our vector above
      captured_warnings <<- c(captured_warnings, conditionMessage(w))
      # Prevent warnings from coming up and being shown
      invokeRestart("muffleWarning")
    }
  )

  # Verify expected warnings
  expect_true(any(grepl(
    "Adaptive captures only VDJ genes; C genes will be NA.",
    captured_warnings
  )))
  expect_true(any(grepl(
    "These genes are not in IMGT for this species and will be replaced with NA",
    captured_warnings
  )))
})
