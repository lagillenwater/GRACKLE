toy_effect <- function(seed = 1, n_models = 30) {
  withr::with_seed(seed, {
    driver <- stats::rnorm(n_models)
    module <- t(vapply(1:5, function(i) driver + stats::rnorm(n_models, sd = 0.3), numeric(n_models)))
    noise <- matrix(stats::rnorm(12 * n_models, sd = 0.05), 12)
    other <- matrix(stats::rnorm(183 * n_models), 183)
    quiet <- matrix(stats::rnorm(n_models, sd = 0.01), 1)   # far below every control's spread
  })
  E <- rbind(module, noise, other, quiet)
  dimnames(E) <- list(c(paste0("G", 1:5), paste0("N", 1:12), paste0("V", 1:183), "Q1"), paste0("M", seq_len(n_models)))
  E
}

test_that("geneArm parses HGNC cytogenetic locations", {
  expect_equal(geneArm(c("17q21.31", "1p36.33", "Xp22.33 and Yp11.2", "mitochondria", "13cen-q12", NA)),
               c("17q", "1p", "Xp", NA, NA, NA))
})

test_that("hcmiToken extracts the case-level HCMI id", {
  expect_equal(hcmiToken(c("HCM-BROD-0577-C15-85M", "x HCM-CSHL-0143-C20 y", "ACH-000001", NA)),
               c("HCM-BROD-0577-C15", "HCM-CSHL-0143-C20", NA, NA))
})

test_that("informativeGenes keeps variable genes and drops noise-level genes", {
  E <- toy_effect()
  inf <- informativeGenes(E, nonessential = paste0("N", 1:12))
  expect_true(all(paste0("G", 1:5) %in% inf))
  expect_lte(sum(paste0("N", 1:12) %in% inf), 1)   # the controls' own top 5% pass by construction
  expect_false("Q1" %in% inf)
  expect_error(informativeGenes(E, nonessential = "N1"), "fewer than 10")
})

test_that("glsWhiten makes the model covariance the identity", {
  withr::with_seed(2, {
    A <- matrix(stats::rnorm(20 * 20), 20)
    E <- matrix(stats::rnorm(500 * 20), 500) %*% A
  })
  Ew <- glsWhiten(E)
  expect_lt(max(abs(stats::cov(Ew) - diag(20))), 1e-6)
  expect_length(attr(Ew, "intercept"), 20)
})

test_that("GLS correlation removes each gene's mean effect", {
  withr::with_seed(3, {
    A <- matrix(stats::rnorm(25 * 25), 25)                 # correlated models
    E <- matrix(stats::rnorm(300 * 25), 300) %*% A
  })
  E[1:2, ] <- E[1:2, ] - 3   # two strongly essential genes whose profiles are otherwise independent
  rownames(E) <- paste0("g", 1:300)
  C <- glsCorrelation(glsWhiten(E))
  expect_lt(abs(C["g1", "g2"]), 0.5)
})

test_that("coessGraph links module genes, skips same-arm pairs and leaves other genes empty", {
  E <- toy_effect()
  arm <- c(G1 = "1p", G2 = "1p", G3 = "2q", G4 = "3q", G5 = "4p")
  genes <- c(paste0("G", 1:5), "Q1", "V1", "Z")
  G <- coessGraph(E, genes, nonessential = paste0("N", 1:12), arm = arm, k_nn = 3)
  expect_equal(dim(G), c(8, 8))
  expect_true(isSymmetric(G)); expect_false(anyNA(G)); expect_true(all(G >= 0))
  expect_equal(unname(diag(G)), rep(0, 8))
  expect_equal(G["G1", "G2"], 0)
  expect_gt(G["G3", "G4"], 0)
  expect_equal(sum(G["Q1", ]), 0)
  expect_equal(sum(G["Z", ]), 0)
})

test_that("dropHeldOutModels removes only the held-out cases' models", {
  E <- toy_effect()
  model_case <- stats::setNames(c("caseA", "caseB", rep(NA, 28)), colnames(E))
  out <- dropHeldOutModels(E, model_case, heldout_cases = "caseA")
  expect_false("M1" %in% colnames(out))
  expect_true(all(c("M2", "M3") %in% colnames(out)))
  expect_equal(ncol(out), 29)   # 30 models minus caseA's one
})
