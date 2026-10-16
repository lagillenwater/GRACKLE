test_that("enrichHypergeom finds the planted gene set", {
  genes <- sprintf("g%03d", 1:200)
  H <- matrix(0.01, 2, 200, dimnames = list(c("LV1", "LV2"), genes))
  H["LV1", 1:10] <- 1
  H["LV2", 101:110] <- 1
  sets <- list(SET_A = genes[1:20], SET_B = genes[101:120], SET_C = genes[150:170])
  e <- enrichHypergeom(H, sets, frac = 0.05)
  expect_named(e, c("program", "gene_set", "overlap", "set_size", "n_top", "p", "fdr"))
  top <- e[order(e$p), ]
  expect_equal(top$gene_set[top$program == "LV1"][1], "SET_A")
  expect_equal(top$gene_set[top$program == "LV2"][1], "SET_B")
  expect_equal(e$overlap[e$program == "LV1" & e$gene_set == "SET_A"], 10)
})

test_that("programAssociation is NA-safe and finds a negative purity association", {
  withr::with_seed(10, {
    purity <- stats::runif(50, 0.2, 0.9)
    act <- cbind(LV1 = 1 - purity + stats::rnorm(50, sd = 0.05), LV2 = stats::rnorm(50), LV3 = 0)
  })
  act[1, "LV2"] <- NA
  a <- programAssociation(act, purity, "purity_dna")
  expect_lt(a$rho[1], -0.8); expect_lt(a$fdr[1], 0.01)
  expect_equal(a$rho[3], 0); expect_equal(a$p[3], 1)
  expect_false(anyNA(a$p))
})

test_that("mutationMatrix keeps non-silent calls and marks unsequenced samples NA", {
  mut <- data.frame(sample_id = c("s1", "s1", "s2", "s3"),
                    Hugo_Symbol = c("KRAS", "TP53", "KRAS", "KRAS"),
                    Variant_Classification = c("Missense_Mutation", "Silent", "Nonsense_Mutation", "Missense_Mutation"))
  M <- mutationMatrix(mut, row_ids = c("p1", "p2", "p3", "p4"), row_sample = c("s1", "s2", "s3", "s4"),
                      genes = c("KRAS", "TP53"), profiled_samples = c("s1", "s2", "s3"))
  expect_equal(M[, "KRAS"], c(p1 = TRUE, p2 = TRUE, p3 = TRUE, p4 = NA))
  expect_equal(M[, "TP53"], c(p1 = FALSE, p2 = FALSE, p3 = FALSE, p4 = NA))
})

test_that("genotypeAssociation tests only genes with enough mutated and wild-type tumors", {
  withr::with_seed(11, {
    m <- c(rep(TRUE, 15), rep(FALSE, 35))
    act <- cbind(LV1 = ifelse(m, 2, 0) + stats::rnorm(50, sd = 0.2), LV2 = stats::rnorm(50))
  })
  mut <- cbind(KRAS = m, RARE = c(rep(TRUE, 3), rep(FALSE, 47)))
  g <- genotypeAssociation(act, mut, min_mut = 10)
  expect_equal(unique(g$gene), "KRAS")
  expect_lt(g$fdr[g$program == "LV1"], 0.01)
  expect_gt(g$delta[g$program == "LV1"], 1)
})

test_that("estimateScores returns stromal and immune scores per sample", {
  skip_if_not_installed("tidyestimate")
  genes <- unique(tidyestimate::common_genes$hgnc_symbol)
  withr::with_seed(12, expr <- matrix(stats::rnorm(length(genes) * 4, 5), length(genes),
                                     dimnames = list(genes, paste0("s", 1:4))))
  sc <- estimateScores(expr)
  expect_named(sc, c("profile_id", "stromal", "immune", "estimate"))
  expect_equal(sc$profile_id, paste0("s", 1:4))
})
