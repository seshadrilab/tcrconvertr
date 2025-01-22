
test_that("can parse an imgt fasta", {
  fasta <- system.file("extdata/examples/fasta_dir", "test_trav.fa", package = "TCRconvertR")
  expected_out <- c("TRAV1-1*01", "TRAV1-1*02", "TRAV14/DV4*01", "TRAV38-2/DV8*01", "TRAC*01")
  names(expected_out) <- c('>AE000658|TRAV1-1*01|Homo sapiens|F|V-REGION|128090..128364|275 nt|1| | | | |275+0=275| | |',
                           '>X04939|TRAV1-1*02|Homo sapiens|(F)|V-REGION|52..320|269 nt|1| | | | |269+0=269| | |',
                           '>M21626|TRAV14/DV4*01|Homo sapiens|F|V-REGION|226..515|290 nt|1| | | | |290+0=290| | |', 
                           '>AE000661|TRAV38-2/DV8*01|Homo sapiens|F|V-REGION|34309..34597|289 nt|1| | | | |289+0=289| | |',
                           '>X02883|TRAC*01|Homo sapiens|F|EX1|n,273..544|273 nt|1|+1|-1| | |273+0=273| | |')
  expect_equal(parse_imgt_fasta(fasta), expected_out)
})

test_that("can extract imgt genes", {
  fastadir <- system.file("extdata/examples/fasta_dir", package = "TCRconvertR")
  expected_out <- data.frame(imgt = c('TRAC*01',
                                      'TRAV1-1*01',
                                      'TRAV1-1*02',
                                      'TRAV14/DV4*01',
                                      'TRAV38-2/DV8*01',
                                      'TRBV29/OR9-2*01',
                                      'TRBVA/OR9-2*01'))
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


# TODO: Make it so files are written to a temp dir (withr::local_tempfile() ?)
test_that("can build lookup tables from fastas", {
  fastadir <- system.file("extdata/examples/fasta_dir", package = "TCRconvertR")

  # Create lookup tables
  build_lookup_from_fastas(fastadir)

  # Check adaptive lookup table
  adapt <- read.csv(system.file("extdata/examples/fasta_dir", "lookup_from_adaptive.csv", package = "TCRconvertR"))
  expected_adapt <- data.frame(adaptive = c("TCRAV01-01*01",
                                            "TCRAV01-01*02",
                                            "TCRAV14-01*01",
                                            "TCRAV38-02*01",
                                            "TCRBV29-or09_02*01",
                                            "TCRBVA-or09_02*01",
                                            "TCRAV01-01",
                                            "TCRAV14-01",
                                            "TCRAV38-02",
                                            "TCRBV29-or09_02",
                                            "TCRBVA-or09_02"),
                               adaptivev2 = c("TCRAV01-01*01",
                                              "TCRAV01-01*02",
                                              "TCRAV14-01*01",
                                              "TCRAV38-02*01",
                                              "TCRBV29-or09_02*01",
                                              "TCRBVA-or09_02*01",
                                              "TCRAV01-01",
                                              "TCRAV14-01",
                                              "TCRAV38-02",
                                              "TCRBV29-or09_02",
                                              "TCRBVA-or09_02"),
                               imgt = c("TRAV1-1*01",
                                        "TRAV1-1*02",
                                        "TRAV14/DV4*01",
                                        "TRAV38-2/DV8*01",
                                        "TRBV29/OR9-2*01",
                                        "TRBVA/OR9-2*01",
                                        "TRAV1-1*01",
                                        "TRAV14/DV4*01",
                                        "TRAV38-2/DV8*01",
                                        "TRBV29/OR9-2*01",
                                        "TRBVA/OR9-2*01"),
                               tenx = c("TRAV1-1",
                                        "TRAV1-1",
                                        "TRAV14DV4",
                                        "TRAV38-2DV8",
                                        "TRBV29/OR9-2",
                                        "TRBVA/OR9-2",
                                        "TRAV1-1",
                                        "TRAV14DV4",
                                        "TRAV38-2DV8",
                                        "TRBV29/OR9-2",
                                        "TRBVA/OR9-2"))
  expect_equal(adapt, expected_adapt)
  
  # Check 10X lookup table
  tenx <- read.csv(system.file("extdata/examples/fasta_dir", "lookup_from_tenx.csv", package = "TCRconvertR"))
  expected_tenx <- data.frame(tenx = c("TRAC",
                                       "TRAV1-1",
                                       "TRAV14DV4",
                                       "TRAV38-2DV8",
                                       "TRBV29/OR9-2",
                                       "TRBVA/OR9-2"),
                              imgt = c("TRAC*01",
                                       "TRAV1-1*01",
                                       "TRAV14/DV4*01",
                                       "TRAV38-2/DV8*01",
                                       "TRBV29/OR9-2*01",
                                       "TRBVA/OR9-2*01"),
                              adaptive = c("NoData",
                                           "TCRAV01-01*01",
                                           "TCRAV14-01*01",
                                           "TCRAV38-02*01",
                                           "TCRBV29-or09_02*01",
                                           "TCRBVA-or09_02*01"),
                              adaptivev2 = c("NoData",
                                             "TCRAV01-01*01",
                                             "TCRAV14-01*01",
                                             "TCRAV38-02*01",
                                             "TCRBV29-or09_02*01",
                                             "TCRBVA-or09_02*01"))
  expect_equal(tenx, expected_tenx)
  
  # Check regular lookup table
  lookup <- read.csv(system.file("extdata/examples/fasta_dir", "lookup.csv", package = "TCRconvertR"))
  expected_lookup <- data.frame(imgt = c("TRAC*01",
                                         "TRAV1-1*01",
                                         "TRAV1-1*02",
                                         "TRAV14/DV4*01",
                                         "TRAV38-2/DV8*01",
                                         "TRBV29/OR9-2*01",
                                         "TRBVA/OR9-2*01"),
                                tenx = c("TRAC",
                                         "TRAV1-1",
                                         "TRAV1-1",
                                         "TRAV14DV4",
                                         "TRAV38-2DV8",
                                         "TRBV29/OR9-2",
                                         "TRBVA/OR9-2"),
                                adaptive = c("NoData",
                                             "TCRAV01-01*01",
                                             "TCRAV01-01*02",
                                             "TCRAV14-01*01",
                                             "TCRAV38-02*01",
                                             "TCRBV29-or09_02*01",
                                             "TCRBVA-or09_02*01"),
                                adaptivev2 = c("NoData",
                                               "TCRAV01-01*01",
                                               "TCRAV01-01*02",
                                               "TCRAV14-01*01",
                                               "TCRAV38-02*01",
                                               "TCRBV29-or09_02*01",
                                               "TCRBVA-or09_02*01"))
  expect_equal(lookup, expected_lookup)
})
