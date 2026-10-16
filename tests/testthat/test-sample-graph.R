test_that("knnSparsify keeps each row's k largest positive similarities, then symmetrizes", {
  S <- matrix(c(0, .9, .8, .1,
                .9, 0, .2, 0,
                .8, .2, 0, .3,
                .1, 0, .3, 0), 4, byrow = TRUE, dimnames = list(letters[1:4], letters[1:4]))
  K <- knnSparsify(S, 1)
  expect_equal(K["a", "b"], .9); expect_equal(K["b", "a"], .9)
  expect_equal(K["a", "c"], .8); expect_equal(K["c", "d"], .3)
  expect_equal(K["a", "d"], 0)
  expect_true(isSymmetric(K)); expect_equal(unname(diag(K)), rep(0, 4))
})

test_that("pair edges link each organoid only to its own parent tumor", {
  P <- pairGraph(toy_samples())
  expect_equal(P["o1", "t1"], 1)
  expect_equal(P["o2", "t2b"], 1); expect_equal(P["o2", "t2a"], 0)
  expect_equal(P["o3", "t3a"], 1); expect_equal(P["o3", "t3b"], 1)
  expect_equal(P["o6", "t6"], 1)
  expect_equal(sum(P[upper.tri(P)]), 5)
  org <- toy_samples()$role == "organoid"
  expect_equal(sum(P[org, org]), 0)
})

test_that("pairTable lists the same pairs as pairGraph", {
  pt <- pairTable(toy_samples())
  expect_setequal(paste(pt$tumor, pt$organoid), c("t1 o1", "t2b o2", "t3a o3", "t3b o3", "t6 o6"))
})

test_that("purity edges never touch organoids or tumors without purity", {
  s <- toy_samples()
  Pu <- purityGraph(s, k_nn = 2)
  expect_equal(sum(Pu[s$role == "organoid", ]), 0)
  expect_equal(sum(Pu["t3b", ]), 0)
  expect_gt(sum(Pu), 0)
  expect_true(isSymmetric(Pu))
})

test_that("the annotation graph handles missing and constant annotations", {
  A <- annotGraph(toy_samples(), k_nn = 3)
  expect_false(anyNA(A))
  expect_true(isSymmetric(A))
  expect_equal(unname(diag(A)), rep(0, nrow(A)))
  expect_gt(sum(A), 0)
})

test_that("buildSampleGraph is symmetric, NA-free and scaled per part", {
  g <- buildSampleGraph(toy_samples(), k_nn = 3)
  expect_true(isSymmetric(g$S))
  expect_false(anyNA(g$S))
  expect_equal(unname(diag(g$S)), rep(0, 12))
  expect_equal(g$summary$part, c("pair", "annot", "purity"))
  expect_equal(g$summary$mean_degree, c(1, 1, 1))
})

test_that("the graph does not depend on row order", {
  s <- toy_samples()
  g1 <- buildSampleGraph(s, k_nn = 3)$S
  g2 <- buildSampleGraph(s[12:1, ], k_nn = 3)$S
  expect_equal(g2[rownames(g1), colnames(g1)], g1)
})

test_that("a graph built from a training subset contains only those samples", {
  s <- toy_samples()
  train <- s[!s$case_id %in% c("c2", "c6"), ]
  g <- buildSampleGraph(train, k_nn = 3)$S
  expect_equal(rownames(g), train$profile_id)
  expect_false(any(c("o2", "o6", "t2a") %in% rownames(g)))
})

test_that("knnSparsify breaks ties by the byte order of the ids, whatever the collation locale", {
  ids <- c("b", "a", "C")
  S <- matrix(1, 3, 3, dimnames = list(ids, ids))
  K <- withr::with_collate("C", knnSparsify(S, 1))
  ## byte order is "C" < "a" < "b": b and a keep C, C keeps a
  expect_equal(unname(K[cbind(c("b", "a", "C"), c("C", "C", "a"))]), c(1, 1, 1))
  expect_equal(K["a", "b"], 0)
  skip_if_not(collationAvailable(), "en_US.UTF-8 collation not available")
  expect_identical(withr::with_collate("en_US.UTF-8", knnSparsify(S, 1)), K)
})
