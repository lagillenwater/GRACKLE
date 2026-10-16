## A failed sensitivity combine wrote ERROR.txt, and the later successful one wrote its results
## beside it (coess_3d k2, run 20261016_6419dd7, job 34351486), so each run clears the directory
## first and leaves one outcome.
test_that("the sensitivity combine clears its directory before writing an outcome", {
  f <- testthat::test_path("..", "..", "analysis", "hcmi_persistence", "06_combine.R")
  skip_if_not(file.exists(f))
  text <- paste(readLines(f, warn = FALSE), collapse = "\n")
  block <- regmatches(text, regexpr("args\\$stage == \"sensitivity\".*?\n}\n", text))
  expect_length(block, 1)
  clear <- regexpr("unlink(d, recursive = TRUE)", block, fixed = TRUE)
  combine <- regexpr("combineArm(", block, fixed = TRUE)
  expect_gt(clear, 0)
  expect_lt(clear, combine)
})
