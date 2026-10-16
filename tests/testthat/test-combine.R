build_case <- function() {
  w <- fake_world()
  res <- c(fake_results(w, "s1", 0.01, noise = 0.08, recon = c(.100, .104, .096, .102, .098), seed = 1),
           fake_results(w, "s2", 0.10, noise = 0.02, recon = c(.101, .103, .099, .102, .100), seed = 2),
           fake_results(w, "s3", 0.30, noise = 0.02, recon = rep(.090, 5), converged = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE), seed = 3),
           fake_results(w, "s4", 1.00, noise = 0.02, recon = rep(.080, 5), missing_folds = 3L, seed = 4))
  expected <- data.frame(setting_id = c("s1", "s2", "s3", "s4", "s5"), alpha_1 = c(0.01, 0.1, 0.3, 1, 3), alpha_2 = 0)
  list(w = w, res = res, expected = expected)
}

test_that("settingTable flags missing folds, failed folds and non-convergence", {      # Review Focus 5
  b <- build_case()
  failed <- b$res[[1]]; failed$fit <- NULL; failed$error <- "boom"; failed$meta$setting_id <- "s5"
  tab <- settingTable(c(b$res, list(failed)), b$expected)
  expect_equal(tab$complete, c(TRUE, TRUE, TRUE, FALSE, FALSE))
  expect_equal(tab$all_converged, c(TRUE, TRUE, FALSE, FALSE, FALSE))
  expect_equal(tab$missing_or_failed[tab$setting_id == "s4"], "3")
  expect_true(is.na(tab$recon_mean[tab$setting_id == "s4"]))
})

test_that("assembleOutOfFold undoes each fold's program order", {
  b <- build_case()
  fits <- Filter(function(r) r$meta$setting_id == "s2", b$res)
  m <- matchSetting(fits, n_perm = 100)
  oof <- assembleOutOfFold(Filter(function(r) r$meta$fold > 0, fits), m$matches)
  expect_equal(oof[rownames(oof), "LV1"], b$w$act[rownames(oof), 1])
  expect_equal(oof[rownames(oof), "LV2"], b$w$act[rownames(oof), 2])
})

test_that("combineArm selects, scores and excludes incomplete settings without erroring", {
  b <- build_case()
  out <- combineArm(b$res, b$expected, b$w$samples, b$w$pairs, gene_sets = list(),
                    n_perm = 199, n_boot = 50, n_perm_match = 100, seed = 1)
  expect_equal(out$selected, "s2")
  expect_false(out$selection$all_converged[out$selection$setting_id == "s3"])
  expect_false(out$selection$complete[out$selection$setting_id == "s4"])
  grp <- stats::setNames(out$programs$group, out$programs$program)
  expect_equal(grp[["LV1"]], "persists")
  expect_equal(grp[["LV2"]], "drops_out")
  expect_true(all(c("conc_rho", "ret_median_log2", "conc_fdr", "stable", "purity_rho", "microenv") %in% names(out$programs)))
  expect_setequal(unique(out$maps_long$setting_id), c("s1", "s2", "s3"))
  expect_true(all(out$per_type$n_pairs >= 10))
  expect_equal(nrow(out$convergence), length(b$res))
})

## Combine job 34347703: three string_all settings had no stable program, and paste0("LV", integer(0))
## is "LV", not character(0), so scoring indexed a column that does not exist.
test_that("scoreSetting returns no rows when no program is stable", {
  oof <- matrix(c(1, 2, 3, 4), 2, dimnames = list(c("t1", "o1"), c("LV1", "LV2")))
  samples <- data.frame(profile_id = c("t1", "o1"), cancer_type = "X", fold = 1)
  stability <- data.frame(ref = 1:2, stable = c(FALSE, FALSE))
  out <- scoreSetting(oof, data.frame(tumor = "t1", organoid = "o1"), samples, stability, n_perm = 9, n_boot = 9)
  expect_equal(nrow(out), 0)
})
