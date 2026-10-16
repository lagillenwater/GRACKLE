test_that("assignMax solves the case where greedy matching fails", {
  sim <- matrix(c(0.90, 0.80,
                  0.85, 0.10), 2, byrow = TRUE)   # rows: ref A, B; columns: fold programs 1, 2
  expect_equal(assignMax(sim), c(2L, 1L))
})

test_that("matchPrograms recovers a permutation and flags weak matches", {
  withr::with_seed(7, H_ref <- matrix(stats::runif(4 * 50), 4))
  perm <- c(3, 1, 4, 2)
  H_fold <- H_ref[perm, ]
  m <- matchPrograms(H_ref, H_fold, tau = 0.5)
  expect_equal(m$fold_prog, match(1:4, perm))
  expect_true(all(m$matched))
  expect_equal(m$cosine, rep(1, 4))
  m2 <- matchPrograms(H_ref, H_fold, tau = 1.01)
  expect_false(any(m2$matched))
})

test_that("the gene-permutation null sits below true matches", {
  withr::with_seed(8, H_ref <- matrix(stats::rexp(4 * 200), 4))
  null <- matchNull(H_ref, H_ref + 0.01, n_perm = 50, seed = 1)
  expect_length(null, 4 * 50)
  expect_true(all(null >= 0 & null <= 1))
  expect_lt(stats::quantile(null, 0.99), 0.99)
})

test_that("programStability counts folds matched", {
  matches <- data.frame(fold = rep(1:5, each = 2), ref = rep(1:2, 5),
                        cosine = c(.9, .5, .9, .5, .9, .9, .9, .5, .8, .5),
                        matched = c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE))
  st <- programStability(matches)
  expect_equal(st$n_matched, c(5L, 1L))
  expect_equal(st$stable, c(TRUE, FALSE))
  expect_equal(st$mean_cosine, c(mean(c(.9, .9, .9, .9, .8)), mean(c(.5, .5, .9, .5, .5))))
})
