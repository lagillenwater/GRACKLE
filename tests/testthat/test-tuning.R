test_that("traceLaplacian equals the explicit trace", {
  withr::with_seed(5, {
    S <- matrix(stats::runif(36), 6); S <- (S + t(S)) / 2; diag(S) <- 0
    W <- matrix(stats::runif(18), 6); H <- matrix(stats::runif(18), 3)
  })
  L <- diag(rowSums(S)) - S
  expect_equal(traceLaplacian(W, S, "row"), sum(diag(t(W) %*% L %*% W)))
  expect_equal(traceLaplacian(H, S, "col"), sum(diag(H %*% L %*% t(H))))
})

test_that("lambdaFromAlpha makes each penalty alpha times the starting reconstruction", {
  withr::with_seed(6, {
    Y <- matrix(stats::runif(60), 6)
    S_p <- matrix(stats::runif(36), 6); S_p <- (S_p + t(S_p)) / 2; diag(S_p) <- 0
    S_g <- matrix(stats::runif(100), 10); S_g <- (S_g + t(S_g)) / 2; diag(S_g) <- 0
  })
  init <- svdInit(Y, 2)
  lam <- lambdaFromAlpha(Y, init$W, init$H, S_p, S_g, alpha_1 = 0.1, alpha_2 = 0.3)
  recon <- sum((Y - init$W %*% init$H)^2)
  expect_equal(lam[["lambda_1"]] * traceLaplacian(init$W, S_p, "row"), 0.1 * recon)
  expect_equal(lam[["lambda_2"]] * traceLaplacian(init$H, S_g, "col"), 0.3 * recon)
  zero <- lambdaFromAlpha(Y, init$W, init$H, S_p * 0, S_g, 0.1, 0.3)
  expect_equal(zero[["lambda_1"]], 0)
})

test_that("selectSetting applies the one-SE rule, then reproducibility, then the tie-break", {
  tab <- data.frame(
    setting_id = c("a", "b", "c", "d", "e"),
    alpha_1 = c(0, 0.1, 0.3, 0.1, 0.01), alpha_2 = c(0, 0.1, 0.3, 0.3, 0.1),
    all_converged = c(TRUE, TRUE, TRUE, FALSE, TRUE),
    recon_mean = c(0.100, 0.104, 0.120, 0.090, 0.104),
    recon_se = c(0.005, 0.005, 0.005, 0.005, 0.005),
    reproducibility = c(0.80, 0.90, 0.99, 0.99, 0.90))
  ## d is best but not converged; the band is recon <= 0.105 -> a, b, e; b and e tie on
  ## reproducibility; e has the smaller alpha_1.
  expect_equal(selectSetting(tab), "e")
  expect_error(selectSetting(transform(tab, all_converged = FALSE)), "no setting")
})
