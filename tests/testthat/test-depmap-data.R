hgnc_dm <- data.frame(symbol = c("A1BG", "KRAS"), entrez_id = c(1, 3845), ensembl_gene_id = c("E1", "E2"))

test_that("DepMap gene columns map to HGNC symbols", {
  expect_equal(depmapSymbols(c("A1BG (1)", "KRAS (3845)", "FAKE (999999)", "NOSYM"), hgnc_dm),
               c("A1BG", "KRAS", "FAKE", "NOSYM"))
})

test_that("readDepmapEffect returns genes x models", {
  f <- tempfile(fileext = ".csv")
  writeLines(c(",A1BG (1),KRAS (3845)", "ACH-1,0.1,-1.2", "ACH-2,0.2,-0.8"), f)
  E <- readDepmapEffect(f, hgnc_dm)
  expect_equal(dim(E), c(2, 2))
  expect_equal(rownames(E), c("A1BG", "KRAS"))
  expect_equal(colnames(E), c("ACH-1", "ACH-2"))
  expect_equal(E["KRAS", "ACH-2"], -0.8)
  g <- tempfile(fileext = ".csv")
  writeLines(c("gene,ORG-1,ORG-2", "A1BG,0.1,0.2", "KRAS,-1.2,-0.8"), g)
  E2 <- readDepmapEffect(g, hgnc_dm, genes_in_rows = TRUE)
  expect_equal(rownames(E2), c("A1BG", "KRAS")); expect_equal(colnames(E2), c("ORG-1", "ORG-2"))
  expect_equal(E2["KRAS", "ORG-2"], -0.8)
})

test_that("readNonessential parses the control gene list", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("Gene", "A1BG (1)", "KRAS (3845)"), f)
  expect_no_warning(out <- readNonessential(f, hgnc_dm))   # "SYMBOL (ENTREZ)" must not split on the space
  expect_equal(out, c("A1BG", "KRAS"))
})

test_that("hcmiTokenFromRow finds HCMI ids in any text column", {
  df <- data.frame(ModelID = c("ACH-1", "ACH-2", "ACH-3"),
                   CellLineName = c("HCM-BROD-0416-C71", "HT29", NA),
                   Alias = c(NA, NA, "HCM-SANG-0270-C20-A"), stringsAsFactors = FALSE)
  expect_equal(hcmiTokenFromRow(df), c("HCM-BROD-0416-C71", NA, "HCM-SANG-0270-C20"))
})

test_that("combine3dSources combines agreeing sources and otherwise keeps the larger one", {
  withr::with_seed(16, {
    base <- stats::rnorm(100)
    E1 <- matrix(base + stats::rnorm(100 * 4, sd = 0.2), 100, dimnames = list(paste0("g", 1:100), paste0("dm", 1:4)))
    E2 <- matrix(base + stats::rnorm(100 * 3, sd = 0.2), 100, dimnames = list(paste0("g", 1:100), paste0("sg", 1:3)))
    E3 <- matrix(stats::rnorm(100 * 6), 100, dimnames = list(paste0("g", 1:100), paste0("sg", 1:6)))
  })
  a <- combine3dSources(E1, E2)
  expect_equal(ncol(a$effect), 7); expect_match(a$decision, "combined")
  b <- combine3dSources(E1, E3)
  expect_equal(colnames(b$effect), colnames(E3)); expect_match(b$decision, "Sanger alone")
  expect_equal(combine3dSources(E1)$decision, "DepMap 3D only")
})
