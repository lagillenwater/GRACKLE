make_cases <- function() {
  data.frame(profile_id = paste0("p", 1:40),
             case_id = rep(paste0("c", 1:20), each = 2),
             cancer_type = rep(rep(c("COAD", "PAAD", "BRCA", "SKCM"), times = c(8, 6, 4, 2)), each = 2),
             stringsAsFactors = FALSE)
}

test_that("every case lands in exactly one fold", {
  s <- makeFolds(make_cases())
  per_case <- tapply(s$fold, s$case_id, function(f) length(unique(f)))
  expect_true(all(per_case == 1))
  expect_true(all(s$fold %in% 1:5))
})

test_that("folds are balanced within each cancer type", {
  s <- makeFolds(make_cases())
  cases <- unique(s[, c("case_id", "cancer_type", "fold")])
  balance <- tapply(cases$fold, cases$cancer_type, function(f) diff(range(tabulate(f, nbins = 5))))
  expect_true(all(balance <= 1))
})

test_that("folds are reproducible with the seed and change with it", {
  expect_identical(makeFolds(make_cases())$fold, makeFolds(make_cases())$fold)
  expect_false(identical(makeFolds(make_cases())$fold, makeFolds(make_cases(), seed = 7)$fold))
})

test_that("a case with two cancer types is an error", {
  s <- make_cases(); s$cancer_type[2] <- "PAAD"
  expect_error(makeFolds(s), "more than one cancer_type")
})

test_that("makeFolds leaves the caller's RNG state alone", {
  set.seed(99); before <- runif(1)
  set.seed(99); invisible(makeFolds(make_cases())); after <- runif(1)
  expect_identical(before, after)
})
