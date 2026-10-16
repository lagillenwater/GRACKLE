skip_if_no_tf <- function() {
  testthat::skip_if_not(reticulate::py_module_available("tensorflow"), "TensorFlow not available")
}

toy_problem <- function(n = 20, g = 30, k = 3, seed = 1) {
  withr::with_seed(seed, {
    Y <- matrix(stats::runif(n * g), n, g,
                dimnames = list(paste0("s", seq_len(n)), paste0("g", seq_len(g))))
    S_p <- matrix(stats::rbinom(n * n, 1, 0.2), n, n)
    S_g <- matrix(stats::rbinom(g * g, 1, 0.1), g, g)
  })
  S_p <- pmax(S_p, t(S_p)); diag(S_p) <- 0
  S_g <- pmax(S_g, t(S_g)); diag(S_g) <- 0
  list(Y = Y, S_p = S_p, S_g = S_g, init = svdInit(Y, k), k = k)
}

mu_nmf <- function(Y, W, H, iterations, eps = 1e-10) {
  for (i in seq_len(iterations)) {
    H <- H * (t(W) %*% Y) / (t(W) %*% W %*% H + eps)
    W <- W * (Y %*% t(H)) / (W %*% H %*% t(H) + eps)
  }
  list(W = W, H = H)
}
