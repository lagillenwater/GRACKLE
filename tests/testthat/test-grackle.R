test_that("svdInit returns nonnegative factors scaled to [0, 1]", {
  p <- toy_problem()
  expect_equal(dim(p$init$W), c(20, 3))
  expect_equal(dim(p$init$H), c(3, 30))
  expect_true(all(p$init$W >= 0) && all(p$init$H >= 0))
  expect_equal(unname(apply(p$init$W, 2, max)), rep(1, 3))
  expect_equal(unname(apply(p$init$H, 1, max)), rep(1, 3))
})

test_that("the objective does not increase across iterations", {
  skip_if_no_tf()
  p <- toy_problem()
  traces <- lapply(c(0, 0.1), function(lam) {
    GRACKLE(p$Y, p$S_g, p$S_p, diff_threshold = 0, lambda_1 = lam, lambda_2 = lam, k = p$k,
            iterations = 50, svd_W = p$init$W, svd_H = p$init$H, dtype = "float64")$loss_trace$total
  })
  lapply(traces, function(tr) expect_true(all(diff(tr) <= 1e-9 * tr[-1])))
})

test_that("loss_trace reports the true objective parts", {
  skip_if_no_tf()
  p <- toy_problem()
  fit <- GRACKLE(p$Y, p$S_g, p$S_p, diff_threshold = 0, lambda_1 = 0.1, lambda_2 = 0.1, k = p$k,
                 iterations = 5, svd_W = p$init$W, svd_H = p$init$H, dtype = "float64")
  lt <- fit$loss_trace
  expect_named(lt, c("iteration", "reconstruction", "sample_graph", "gene_graph", "total"))
  expect_equal(lt$total, lt$reconstruction + lt$sample_graph + lt$gene_graph, tolerance = 1e-10)
  L_p <- diag(rowSums(p$S_p)) - p$S_p
  L_g <- diag(rowSums(p$S_g)) - p$S_g
  W <- unname(fit$W); H <- unname(fit$H)
  expect_equal(lt$reconstruction[5], sum((p$Y - W %*% H)^2), tolerance = 1e-8)
  expect_equal(lt$sample_graph[5], 0.1 * sum(diag(t(W) %*% L_p %*% W)), tolerance = 1e-8)
  expect_equal(lt$gene_graph[5], 0.1 * sum(diag(H %*% L_g %*% t(H))), tolerance = 1e-8)
})

test_that("iterations is honored and convergence is reported", {
  skip_if_no_tf()
  p <- toy_problem()
  fit <- GRACKLE(p$Y, p$S_g, p$S_p, diff_threshold = 0, lambda_1 = 0.1, lambda_2 = 0.1, k = p$k,
                 iterations = 7, svd_W = p$init$W, svd_H = p$init$H)
  expect_equal(fit$n_iter, 7L)
  expect_false(fit$converged)
  expect_equal(nrow(fit$loss_trace), 7L)
  fit2 <- GRACKLE(p$Y, p$S_g, p$S_p, diff_threshold = 1, lambda_1 = 0.1, lambda_2 = 0.1, k = p$k,
                  iterations = 50, svd_W = p$init$W, svd_H = p$init$H)
  expect_equal(fit2$n_iter, 1L)
  expect_true(fit2$converged)
})

test_that("dtype is honored", {
  skip_if_no_tf()
  p <- toy_problem()
  run <- function(dt) GRACKLE(p$Y, p$S_g, p$S_p, diff_threshold = 0, lambda_1 = 0.1, lambda_2 = 0.1,
                              k = p$k, iterations = 20, svd_W = p$init$W, svd_H = p$init$H, dtype = dt)
  f64 <- run("float64"); f32 <- run("float32"); f16 <- run("float16")
  expect_equal(f32$W, f64$W, tolerance = 1e-3)
  expect_false(isTRUE(all.equal(f32$W, f64$W, tolerance = 1e-12)))
  expect_true(all(is.finite(f16$loss_trace$total)))
  expect_error(run("int8"))
})

test_that("with lambda = 0 GRACKLE equals plain multiplicative-update NMF", {
  skip_if_no_tf()
  p <- toy_problem()
  fit <- GRACKLE(p$Y, p$S_g, p$S_p, diff_threshold = 0, lambda_1 = 0, lambda_2 = 0, k = p$k,
                 iterations = 25, svd_W = p$init$W, svd_H = p$init$H, dtype = "float64")
  ref <- mu_nmf(p$Y, p$init$W, p$init$H, 25)
  ## the multiplicative updates are scale-equivariant, so unit-norm rows of H leave W %*% H unchanged
  expect_equal(unname(fit$W %*% fit$H), unname(ref$W %*% ref$H), tolerance = 1e-8)
})

test_that("dimnames carry sample and gene names", {
  skip_if_no_tf()
  p <- toy_problem()
  fit <- GRACKLE(p$Y, p$S_g, p$S_p, k = p$k, iterations = 2, svd_W = p$init$W, svd_H = p$init$H)
  expect_equal(rownames(fit$W), rownames(p$Y))
  expect_equal(colnames(fit$H), colnames(p$Y))
  expect_equal(colnames(fit$W), paste0("LV", 1:3))
})

## W %*% H is unchanged by W * c, H / c, so a penalty on one factor alone was escaped by shrinking that
## factor (HCMI run 20261016_6419dd7). Rows of H are held at unit norm, so the penalty must smooth H.
test_that("a gene-graph penalty alone smooths H instead of shrinking it", {
  skip_if_no_tf()
  p <- toy_problem()
  fit <- function(l2) GRACKLE(p$Y, p$S_g, p$S_p, diff_threshold = 0, lambda_1 = 0, lambda_2 = l2, k = p$k,
                              iterations = 500, svd_W = p$init$W, svd_H = p$init$H, dtype = "float64")
  L <- diag(rowSums(p$S_g)) - p$S_g
  rough <- function(H) sum(diag(H %*% L %*% t(H))) / sum(H^2)
  plain <- fit(0)
  smooth <- fit(5)
  expect_equal(unname(sqrt(rowSums(smooth$H^2))), rep(1, p$k), tolerance = 1e-8)
  expect_lt(rough(unname(smooth$H)), 0.8 * rough(unname(plain$H)))
})
