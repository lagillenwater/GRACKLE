make_pairs <- function(n = 60, shared = TRUE, seed = 9) {
  withr::with_seed(seed, {
    a <- stats::rgamma(n, 2, 1)
    t_act <- a + stats::rnorm(n, sd = 0.1)
    o_act <- if (shared) a + stats::rnorm(n, sd = 0.1) else stats::rgamma(n, 2, 1)
  })
  data.frame(tumor_act = pmax(t_act, 0), organoid_act = pmax(o_act, 0),
             cancer_type = rep(c("COAD", "PAAD", "BRCA"), length.out = n),
             fold = rep(1:5, each = 1, length.out = n))
}

test_that("shared activity is concordant and independent activity is not", {
  hit <- persistenceConcordance(make_pairs(shared = TRUE), n_perm = 999, n_boot = 100)
  miss <- persistenceConcordance(make_pairs(shared = FALSE), n_perm = 999, n_boot = 100)
  expect_lt(hit$conc_p, 0.01); expect_gt(hit$conc_rho, 0.8)
  expect_gt(miss$conc_p, 0.01)
  expect_true(hit$conc_ci_low < hit$conc_rho && hit$conc_rho < hit$conc_ci_high)
  expect_lt(hit$detect_limit, 0.5)
})

test_that("concordance p-values are calibrated under the null", {
  p <- vapply(1:200, function(s) persistenceConcordance(make_pairs(shared = FALSE, seed = s),
                                                        n_perm = 199, n_boot = 2, seed = s)$conc_p, numeric(1))
  expect_true(mean(p < 0.05) > 0.01 && mean(p < 0.05) < 0.10)
})

test_that("all-zero organoid activity gives rho 0 and p 1, never NA", {          # Review Focus 1
  pr <- make_pairs(); pr$organoid_act <- 0
  r <- persistenceConcordance(pr, n_perm = 99, n_boot = 20)
  expect_equal(r$conc_rho, 0); expect_equal(r$conc_p, 1)
  expect_false(anyNA(unlist(r)))
})

test_that("a stratum of size one stays fixed and shuffles never cross strata", { # Review Focus 2
  strata <- c("a", "b", "b", "c", "b")
  idx <- withr::with_seed(1, permuteWithin(strata))
  expect_equal(idx[c(1, 4)], c(1L, 4L))
  expect_setequal(idx[c(2, 3, 5)], c(2L, 3L, 5L))
  big <- sample(letters[1:6], 300, replace = TRUE)
  expect_identical(big[withr::with_seed(2, permuteWithin(big))], big)
  pr <- make_pairs(); pr$cancer_type[1] <- "UNIQUE"
  expect_false(anyNA(unlist(persistenceConcordance(pr, n_perm = 99, n_boot = 20))))
})

test_that("retention detects dropout and is calibrated when there is none", {
  pr <- make_pairs(); pr$organoid_act <- pr$organoid_act * 0.1
  r <- persistenceRetention(pr, n_perm = 999, n_boot = 100)
  expect_lt(r$ret_median_log2, -2); expect_lt(r$ret_p, 0.01)
  p <- vapply(1:200, function(s) persistenceRetention(make_pairs(shared = FALSE, seed = s),
                                                      n_perm = 199, n_boot = 2, seed = s)$ret_p, numeric(1))
  expect_true(mean(p < 0.05) > 0.01 && mean(p < 0.05) < 0.10)
})

test_that("concordance power is near alpha at rho 0 and high for strong concordance", {
  strata <- data.frame(cancer_type = rep(c("COAD", "PAAD", "BRCA"), length.out = 60), fold = rep(1:5, length.out = 60))
  pw <- concordancePower(strata, rhos = c(0, 0.6), n_sims = 40, n_perm = 99)
  expect_equal(pw$rho, c(0, 0.6))
  expect_lt(pw$power[1], 0.2)
  expect_gt(pw$power[2], 0.8)
  expect_equal(unique(pw$n_pairs), 60L)
})

test_that("classifyPrograms maps the two calls to the four groups", {
  expect_equal(classifyPrograms(conc_fdr = c(0.5, 0.01, 0.01, 0.5),
                                ret_fdr = c(0.01, 0.01, 0.5, 0.5),
                                ret_median_log2 = c(-2, -2, 0, 0)),
               c("drops_out", "reduced_concordant", "persists", "present_not_patient_specific"))
  expect_equal(classifyPrograms(0.5, 0.01, -0.5), "present_not_patient_specific")
})

test_that("seeded nulls and bootstraps do not depend on the collation locale", {
  skip_if_not(collationAvailable(), "en_US.UTF-8 collation not available")
  strata <- rep(c("b", "a", "C"), each = 4)
  for_both <- function(f) lapply(c("C", "en_US.UTF-8"), function(loc) withr::with_collate(loc, withr::with_seed(1, f())))
  expect_identical(for_both(function() permuteWithin(strata))[[1]], for_both(function() permuteWithin(strata))[[2]])
  expect_identical(for_both(function() resampleWithin(strata))[[1]], for_both(function() resampleWithin(strata))[[2]])
  pairs <- withr::with_seed(2, data.frame(tumor_act = runif(24), organoid_act = runif(24),
                                          cancer_type = rep(c("b", "a", "C"), each = 8), fold = rep(1:2, 12)))
  conc <- for_both(function() persistenceConcordance(pairs, n_perm = 200, n_boot = 200, seed = 3))
  ret <- for_both(function() persistenceRetention(pairs, n_perm = 200, n_boot = 200, seed = 3))
  expect_identical(conc[[1]], conc[[2]])
  expect_identical(ret[[1]], ret[[2]])
})
