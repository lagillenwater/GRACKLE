test_that("transferTable pairs folds within each arm", {
  mk <- function(arm, variant, fold, recon) list(meta = list(arm = arm, variant = variant, fold = fold),
                                                 fit = list(recon_test = recon), error = NULL)
  res <- c(lapply(1:5, function(f) mk("a", "selected", f, 0.10 + f / 1000)),
           lapply(1:5, function(f) mk("a", "baseline", f, 0.12 + f / 1000)),
           lapply(1:5, function(f) mk("b", "selected", f, 0.20)),
           lapply(1:5, function(f) mk("b", "baseline", f, 0.20)))
  tt <- transferTable(res)
  expect_equal(tt$arm, c("a", "b"))
  expect_equal(tt$diff_mean, c(-0.02, 0))
  expect_equal(tt$n_folds, c(5L, 5L))
})

test_that("cohenKappa matches a hand calculation", {
  expect_equal(cohenKappa(c("x", "x", "y", "y"), c("x", "y", "y", "y")), 0.5)
  expect_equal(cohenKappa(c("x", "y"), c("x", "y")), 1)
})

test_that("mapAgreement matches programs across arms and scores agreement", {
  withr::with_seed(15, H <- normalizeRows(matrix(stats::rexp(3 * 80), 3, dimnames = list(paste0("LV", 1:3), NULL))))
  pa <- data.frame(program = paste0("LV", 1:3), stable = TRUE, group = c("persists", "drops_out", "persists"), conc_rho = c(.8, .1, .6))
  Hb <- H[c(2, 3, 1), ]; rownames(Hb) <- paste0("LV", 1:3)
  pb <- data.frame(program = paste0("LV", 1:3), stable = TRUE, group = c("drops_out", "persists", "persists"), conc_rho = c(.1, .6, .8))
  ag <- mapAgreement(H, Hb, pa, pb, tau = 0.5)
  expect_equal(ag$summary$n_matched_stable, 3L)
  expect_equal(ag$summary$agreement, 1)
  expect_equal(ag$summary$kappa, 1)
  expect_equal(ag$pairs$program_b, c("LV3", "LV1", "LV2"))
})
