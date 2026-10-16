## httr2's req_retry() retries only transient HTTP statuses unless retry_on_failure = TRUE, so a
## dropped connection ("Connection reset by peer" from GDC, job 34344334) ended step 01 at once.
test_that("every download retry in the pipeline also retries dropped connections", {
  dir <- testthat::test_path("..", "..", "analysis", "hcmi_persistence")
  skip_if_not(dir.exists(dir))
  files <- list.files(dir, pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
  calls <- unlist(lapply(files, function(f) {
    text <- paste(readLines(f, warn = FALSE), collapse = "\n")
    regmatches(text, gregexpr("req_retry\\([^)]*\\)", text))[[1]]
  }))
  expect_gt(length(calls), 0)
  expect_true(all(grepl("retry_on_failure = TRUE", calls, fixed = TRUE)), info = paste(calls, collapse = "; "))
})
