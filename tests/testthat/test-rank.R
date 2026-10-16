planted <- function(r = 3, n = 60, g = 200, seed = 4) {
  withr::with_seed(seed, {
    U <- matrix(stats::runif(n * r), n)
    V <- matrix(0, r, g)
    blocks <- split(seq_len(r * 50), rep(seq_len(r), each = 50))
    V[cbind(rep(seq_len(r), each = 50), unlist(blocks))] <- 5
    Y <- pmax(U %*% V + matrix(stats::rnorm(n * g, sd = 0.3), n), 0)
  })
  Y
}

test_that("parallel analysis recovers a planted rank", {
  pa <- parallelAnalysisRank(planted(3), n_perm = 30, nv = 20)
  expect_equal(pa$k, 3L)
  expect_length(pa$observed, 20)
  expect_length(pa$threshold, 20)
  expect_equal(dim(pa$null), c(30, 20))
})

test_that("parallel analysis is reproducible with the seed", {
  a <- parallelAnalysisRank(planted(2), n_perm = 10, nv = 10)
  b <- parallelAnalysisRank(planted(2), n_perm = 10, nv = 10)
  expect_identical(a$threshold, b$threshold)
})
