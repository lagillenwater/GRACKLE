small_fixture <- function() {
  s <- toy_samples()
  s <- s[rep(seq_len(nrow(s)), 3), ]
  s$profile_id <- paste0(s$profile_id, "_", rep(1:3, each = 12))
  s$case_id <- paste0(s$case_id, "_", rep(1:3, each = 12))
  s$sample_id <- paste0(s$sample_id, "_", rep(1:3, each = 12))
  s$parent_sample_id <- ifelse(is.na(s$parent_sample_id), NA, paste0(s$parent_sample_id, "_", rep(1:3, each = 12)))
  s <- makeFolds(s, seed = 1)
  withr::with_seed(13, Y <- matrix(stats::runif(nrow(s) * 40), nrow(s),
                                   dimnames = list(s$profile_id, paste0("g", 1:40))))
  G <- matrix(0, 40, 40, dimnames = list(colnames(Y), colnames(Y)))
  G[cbind(1:39, 2:40)] <- 1; G <- pmax(G, t(G))
  list(Y = Y, samples = s, G = G)
}

test_that("fitTask holds out the fold and projects held-out samples", {
  skip_if_no_tf()
  fx <- small_fixture()
  r <- fitTask(fx$Y, fx$samples, fx$G, k = 2, alpha_1 = 0.1, alpha_2 = 0.1, fold = 1, iterations = 30)
  held <- fx$samples$profile_id[fx$samples$fold == 1]
  expect_setequal(r$test_ids, held)
  expect_length(intersect(r$train_ids, held), 0)
  expect_equal(rownames(r$W_test), held)
  expect_equal(unname(sqrt(rowSums(r$H^2))), c(1, 1))
  expect_equal(rownames(r$H), c("LV1", "LV2"))
  expect_true(is.finite(r$recon_test) && r$recon_test > 0)
  expect_named(r$lambda, c("lambda_1", "lambda_2"))
})

test_that("the reference fit uses every sample and the transfer roles are respected", {
  skip_if_no_tf()
  fx <- small_fixture()
  ref <- fitTask(fx$Y, fx$samples, fx$G, k = 2, alpha_1 = 0.1, alpha_2 = 0.1, fold = 0, iterations = 10)
  expect_setequal(ref$train_ids, fx$samples$profile_id)
  expect_null(ref$W_test)
  tr <- fitTask(fx$Y, fx$samples, fx$G, k = 2, alpha_1 = 0.1, alpha_2 = 0.1, fold = 2,
                train_role = "tumor", test_role = "organoid", iterations = 10)
  role <- stats::setNames(fx$samples$role, fx$samples$profile_id)
  expect_true(all(role[tr$train_ids] == "tumor"))
  expect_true(all(role[tr$test_ids] == "organoid"))
})

test_that("misaligned inputs are rejected", {
  fx <- small_fixture()
  expect_error(fitTask(fx$Y[, 40:1], fx$samples, fx$G, 2, 0.1, 0.1, 1), "aligned")
})

## Run 20261016_d53b86d: at a 2,000-iteration cap every fit stopped unconverged (relative change
## 4e-6 to 1.2e-4), so the spec's exclusion rule would have removed every setting. Probes on the
## real data converged at 1e-6 within 3,480-10,417 iterations; the cap leaves twice that margin.
test_that("fitTask's iteration cap leaves room for the slowest settings to converge", {
  expect_equal(formals(fitTask)$iterations, 20000)
  expect_equal(formals(fitTask)$diff_threshold, 1e-6)
})

## coess_3d_a04_k2_f4 stopped at the 20,000 cap; rerunning it with a higher cap must refit that
## result and leave every finished or converged one alone.
test_that("needsFit refits only missing, failed, or cap-limited unconverged results", {
  done <- list(fit = list(converged = TRUE, n_iter = 500L), error = NULL)
  capped <- list(fit = list(converged = FALSE, n_iter = 20000L), error = NULL)
  failed <- list(fit = NULL, error = "boom")
  expect_true(needsFit(NULL))
  expect_true(needsFit(failed))
  expect_false(needsFit(done))
  expect_false(needsFit(capped))
  expect_true(needsFit(capped, iterations = 40000L))
  expect_false(needsFit(capped, iterations = 20000L))
  expect_false(needsFit(done, iterations = 40000L))
})

test_that("a task result is written whole, and a truncated one reads as missing", {
  d <- withr::local_tempdir()
  out <- file.path(d, "t.rds")
  saveResult(list(meta = list(task_id = "t"), fit = list(n_iter = 5L, converged = TRUE)), out)
  expect_identical(list.files(d), "t.rds")
  expect_identical(readResult(out)$fit$n_iter, 5L)
  expect_null(readResult(file.path(d, "absent.rds")))
  bytes <- readBin(out, "raw", file.size(out))
  writeBin(bytes[seq_len(length(bytes) %/% 2)], out)   # a task killed while writing
  expect_warning(r <- readResult(out), "unreadable")
  expect_null(r)
  expect_true(needsFit(r))
})
