

test_that("can convert genes", {
  # First, define all the dataframes we need
  imgt_df <- data.frame(v_gene = c('TRAV12-1*01', 'TRBV15*01'),
                        d_gene = c(NA, 'TRBD1*01'),
                        j_gene = c('TRAJ16*0', 'TRBJ2-5*01'),
                        c_gene = c('TRAC*01', 'TRBC2*01'),
                        cdr3 = c('CAVLIF', 'CASSGF'))

  tenx_df <- data.frame(v_gene = c('TRAV12-1', 'TRBV15'),
                        d_gene = c(NA, 'TRBD1'),
                        j_gene = c('TRAJ16', 'TRBJ2-5'),
                        c_gene = c('TRAC', 'TRBC2'),
                        cdr3 = c('CAVLIF', 'CASSGF'))

  adapt_df <- data.frame(v_resolved = c('TCRAV12-01*01', 'TCRBV15-01*01'),
                         d_resolved = c(NA, 'TCRBD01-01*01'),
                         j_resolved = c('TCRAJ16-01*01', 'TCRBJ02-05*01'),
                         cdr3_amino_acid = c('CAVLIF', 'CASSGF'))

  adapt_v2_df <- data.frame(vMaxResolved = c('TCRAV12-01*01', 'TCRBV15-01*01'),
                            dMaxResolved = c(NA, 'TCRBD01-01*01'),
                            jMaxResolved = c('TCRAJ16-01*01', 'TCRBJ02-05*01'),
                            aminoAcid = c('CAVLIF', 'CASSGF'))
  
  custom_df <- imgt_df
  colnames(custom_df) <- c("myV", "myD", "myJ", "myC", "myCDR3")

  tenx_to_adapt_df <- adapt_df
  colnames(tenx_to_adapt_df) <- c("v_gene", "d_gene", "j_gene", "cdr3")
  # Insert a 'c_gene' column before 'cdr3'
  tenx_to_adapt_df <- cbind(
    tenx_to_adapt_df[, 1:3],
    c_gene = NA,
    tenx_to_adapt_df["cdr3"])

  adapt_to_tenx_df <- tenx_df
  colnames(adapt_to_tenx_df) <- c("v_resolved", "d_resolved", "j_resolved", "cdr3_amino_acid")
  adapt_to_tenx_df$c_gene <- NULL

  adapt_to_imgt_df <- imgt_df
  colnames(adapt_to_imgt_df) <- c("v_resolved", "d_resolved", "j_resolved", "cdr3_amino_acid")
  adapt_to_imgt_df$c_gene <- NULL

  adaptv2_to_tenx_df <- tenx_df
  colnames(adaptv2_to_tenx_df) <- c("vMaxResolved", "dMaxResolved", "jMaxResolved", "aminoAcid")
  adaptv2_to_tenx_df <- NULL

  adaptv2_to_imgt_df <- imgt_df
  colnames(adaptv2_to_imgt_df) <- c("vMaxResolved", "dMaxResolved", "jMaxResolved", "aminoAcid")
  adaptv2_to_imgt_df <- NULL

  custom_to_tenx_df <- tenx_df
  colnames(custom_to_tenx_df) <- c("myV", "myD", "myJ", "myC", "myCDR3")

  adapt_no_allele_df <- data.frame(v_resolved = c("TCRAV12-01", "TCRBV15-01*01"),
                                   d_resolved = c(NA, "TCRBD01-01"),
                                   j_resolved = c("TCRAJ16-01*01", "TCRBJ02-05"),
                                   cdr3_amino_acid = c("CAVLIF", "CASSGF"))
  # TODO: Fix this and finish building out all possible parameter combos
  # 10X <-> Adaptive
  expect_equal(convert_gene(tenx_df, 'tenx', 'adaptive'), tenx_to_adapt_df)
})


test_that("bad input raises error", {
  tenx_df <- data.frame(v_gene = c('TRAV12-1', 'TRBV15'),
                        d_gene = c(NA, 'TRBD1'),
                        j_gene = c('TRAJ16', 'TRBJ2-5'),
                        c_gene = c('TRAC', 'TRBC2'),
                        cdr3 = c('CAVLIF', 'CASSGF'))

  # Same input and output format
  expect_error(convert_gene(tenx_df, 'tenx', 'tenx'), '"frm" and "to" formats should be different.')
  # Empty input
  expect_error(convert_gene(data.frame(), 'tenx', 'imgt'), 'Input data is empty.')
})


test_that("chooses correct lookup", {
  lookup_tenx <- system.file("extdata/human", "lookup_from_tenx.csv", package = "TCRconvertR")
  lookup_adapt <- system.file("extdata/human", "lookup_from_adaptive.csv", package = "TCRconvertR")
  lookup_imgt <- system.file("extdata/human", "lookup.csv", package = "TCRconvertR")
  
  # Species we don't have lookups for
  expect_error(choose_lookup('tenx', 'imgt', 'non-existent-species'), 
               "Lookup table not found, please ensure reference files are available.")
  # From different 'frm' formats
  expect_equal(choose_lookup('tenx', 'imgt'), lookup_tenx)
  expect_equal(choose_lookup('adaptivev2', 'imgt'), lookup_adapt)
  expect_equal(choose_lookup('imgt', 'tenx'), lookup_imgt)
})


test_that("chooses correct frm_cols", {
  col_ref <- list(
    adaptive = c("v_resolved", "d_resolved", "j_resolved"),
    adaptivev2 = c("vMaxResolved", "dMaxResolved", "jMaxResolved"),
    imgt = c("v_gene", "d_gene", "j_gene", "c_gene"),
    tenx = c("v_gene", "d_gene", "j_gene", "c_gene"))

  adapt_df <- data.frame(v_resolved = c('TCRAV12-01*01', 'TCRBV15-01*01'),
                         d_resolved = c(NA, 'TCRBD01-01*01'),
                         j_resolved = c('TCRAJ16-01*01', 'TCRBJ02-05*01'),
                         cdr3_amino_acid = c('CAVLIF', 'CASSGF'))

  adapt_v2_df <- data.frame(vMaxResolved = c('TCRAV12-01*01', 'TCRBV15-01*01'),
                            dMaxResolved = c(NA, 'TCRBD01-01*01'),
                            jMaxResolved = c('TCRAJ16-01*01', 'TCRBJ02-05*01'),
                            aminoAcid = c('CAVLIF', 'CASSGF'))

  imgt_df <- data.frame(v_gene = c('TRAV12-1*01', 'TRBV15*01'),
                        d_gene = c(NA, 'TRBD1*01'),
                        j_gene = c('TRAJ16*0', 'TRBJ2-5*01'),
                        c_gene = c('TRAC*01', 'TRBC2*01'),
                        cdr3 = c('CAVLIF', 'CASSGF'))

  tenx_df <- data.frame(v_gene = c('TRAV12-1', 'TRBV15'),
                        d_gene = c(NA, 'TRBD1'),
                        j_gene = c('TRAJ16', 'TRBJ2-5'),
                        c_gene = c('TRAC', 'TRBC2'),
                        cdr3 = c('CAVLIF', 'CASSGF'))

  custom_df <- data.frame(myV = c('TRAV12-1*01', 'TRBV15*01'),
                        myD = c(NA, 'TRBD1*01'),
                        myJ = c('TRAJ16*0', 'TRBJ2-5*01'),
                        myC = c('TRAC*01', 'TRBC2*01'),
                        myCDR3 = c('CAVLIF', 'CASSGF'))

  expect_equal(which_frm_cols(adapt_df, 'adaptive'), col_ref$adaptive)
  expect_equal(which_frm_cols(adapt_v2_df, 'adaptivev2'), col_ref$adaptivev2)
  expect_equal(which_frm_cols(imgt_df, 'imgt'), col_ref$imgt)
  expect_equal(which_frm_cols(tenx_df, 'tenx'), col_ref$tenx)
  # Custom columns
  custom_col <- c('myV', 'myD', 'myJ', 'myC')
  expect_equal(which_frm_cols(custom_df, 'tenx', frm_cols = custom_col), custom_col)
  # Non-existent column
  expect_error(which_frm_cols(tenx_df, 'tenx', frm_cols = c('v_gene', 'j_gene', 'x_gene')),
               "These columns are not in the input dataframe: x_gene")
})

