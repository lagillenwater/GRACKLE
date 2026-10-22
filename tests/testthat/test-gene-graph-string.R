toy_links <- function() {
  data.frame(protein1 = c("P1", "P1", "P2", "P3", "P1", "P5"),
             protein2 = c("P2", "P3", "P1", "P4", "P4", "P1"),
             neighborhood = 0, neighborhood_transferred = 0, fusion = 0, cooccurence = 0, homology = 0,
             coexpression = c(0, 500, 0, 0, 900, 0), coexpression_transferred = 0,
             experiments = c(400, 0, 400, 300, 0, 800), experiments_transferred = 0,
             database = c(0, 0, 0, 900, 0, 0), database_transferred = 0,
             textmining = c(200, 300, 200, 0, 0, 0), textmining_transferred = 0,
             combined_score = c(500, 600, 500, 930, 900, 800))
}
## P2 and P3 both map to B; P4 maps to C and C2; P5 maps to nothing.
toy_map <- data.frame(string_protein_id = c("P1", "P2", "P3", "P4", "P4"),
                      symbol = c("A", "B", "B", "C", "C2"), stringsAsFactors = FALSE)

test_that("stringCombinedAll reproduces STRING's combined_score", {
  x <- readStringLinks(test_path("fixtures", "string_links_full_sample.txt"))
  expect_gt(nrow(x), 2000)
  expect_lte(max(abs(stringCombinedAll(x) - x$combined_score)), 1)
})

test_that("the mechanistic score never exceeds combined_score and is 0 without evidence", {
  x <- readStringLinks(test_path("fixtures", "string_links_full_sample.txt"))
  m <- stringRecombine(x)
  expect_true(all(m * 1000 <= x$combined_score + 1))
  none <- x$experiments + x$experiments_transferred + x$database + x$database_transferred == 0
  expect_true(all(m[none] == 0))
})

test_that("stringRecombine follows STRING's prior-corrected combination", {
  m <- stringRecombine(toy_links())
  e <- (0.3 - 0.041) / 0.959; d <- (0.9 - 0.041) / 0.959
  expect_equal(m[4], (1 - (1 - e) * (1 - d)) * 0.959 + 0.041)
  expect_equal(m[1], 0.4)
  expect_equal(m[2], 0)
})

test_that("stringGraph maps proteins to genes, aggregates by max, and keeps unmapped genes empty", {
  genes <- c("A", "B", "C", "D")
  mech <- stringGraph(toy_links(), toy_map, genes, "mech")
  expect_equal(dim(mech), c(4, 4))
  expect_equal(mech["A", "B"], 0.4)
  expect_equal(mech["B", "C"], (1 - (1 - (0.3 - 0.041) / 0.959) * (1 - (0.9 - 0.041) / 0.959)) * 0.959 + 0.041)
  expect_equal(mech["A", "C"], 0)
  expect_equal(sum(mech["D", ]), 0)
  expect_true(isSymmetric(mech)); expect_false(anyNA(mech)); expect_equal(unname(diag(mech)), rep(0, 4))
  all_ch <- stringGraph(toy_links(), toy_map, genes, "all")
  expect_equal(all_ch["A", "B"], 0.6)
  expect_equal(all_ch["A", "C"], 0.9)
})

test_that("stringProteinToSymbol maps through Ensembl gene ids", {
  aliases <- data.frame(string_protein_id = c("9606.P1", "9606.P1", "9606.P2"),
                        alias = c("ENSG01", "TP53X", "ENSG02"),
                        source = c("Ensembl_gene", "Ensembl_HGNC_symbol", "Ensembl_gene"))
  hgnc <- data.frame(symbol = c("TP53", "KRAS"), ensembl_gene_id = c("ENSG01", "ENSG03"))
  m <- stringProteinToSymbol(aliases, hgnc)
  expect_equal(m$string_protein_id, "9606.P1")
  expect_equal(m$symbol, "TP53")
})

test_that("graphSummary counts genes, edges and mean degree", {
  G <- matrix(0, 3, 3); G[1, 2] <- G[2, 1] <- 0.5
  s <- graphSummary(G)
  expect_equal(s$n_genes_with_edges, 2); expect_equal(s$n_edges, 1); expect_equal(s$mean_degree, 0.5)
})

## The database channel imports curated pathways (Reactome among them), and a graph that keeps every
## weak edge is one tangle; a modular graph from experimental evidence avoids both.
test_that("stringGraph's exp mode uses only the experiments channel, and min_score drops weak edges", {
  genes <- c("A", "B", "C", "C2", "D")
  exp <- stringGraph(toy_links(), toy_map, genes, "exp")
  expect_equal(exp["A", "B"], 0.4, tolerance = 1e-9)   # experiments 400
  expect_equal(exp["B", "C"], 0.3, tolerance = 1e-9)   # experiments 300; its database 900 is ignored
  expect_equal(exp["A", "C"], 0)                       # coexpression only
  strong <- stringGraph(toy_links(), toy_map, genes, "exp", min_score = 0.35)
  expect_equal(strong["A", "B"], 0.4, tolerance = 1e-9)
  expect_equal(strong["B", "C"], 0)
  expect_true(isSymmetric(strong))
})
