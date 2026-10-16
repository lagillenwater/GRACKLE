ralpine_path <- function() testthat::test_path("..", "..", "scripts", "alpine", "ralpine")
ralpine_text <- function() paste(readLines(ralpine_path(), warn = FALSE), collapse = "\n")
first_group <- function(pattern, text) {
  m <- regmatches(text, regexec(pattern, text, perl = TRUE))[[1]]
  if (length(m) < 2) NA_character_ else m[2]
}
case_block <- function(verb) first_group(sprintf("(?s)\\n  %s\\)\\n(.*?)\\n    ;;\\n", verb), ralpine_text())

test_that("READ_ONLY excludes every escape-capable command", {
  skip_if_not(file.exists(ralpine_path()))
  block <- first_group("(?s)READ_ONLY=\\(\\n(.*?)\\n\\)", ralpine_text())
  expect_false(is.na(block))
  tokens <- strsplit(trimws(block), "\\s+")[[1]]
  escape <- c("find", "scontrol", "file", "nvidia-smi", "git", "bash", "sh", "xargs", "python", "perl", "R", "Rscript")
  expect_equal(intersect(escape, tokens), character(0))
})

test_that("log validates its pattern before calling remote_fixed", {
  skip_if_not(file.exists(ralpine_path()))
  b <- case_block("log")
  expect_false(is.na(b))
  v <- regexpr("^[A-Za-z0-9._-]*$", b, fixed = TRUE)
  r <- regexpr('remote_fixed "', b, fixed = TRUE)
  expect_true(v > 0 && r > 0 && v < r)
})

test_that("jobinfo validates the job id before calling remote_fixed and only shows the job", {
  skip_if_not(file.exists(ralpine_path()))
  b <- case_block("jobinfo")
  expect_false(is.na(b))
  v <- regexpr("^[0-9]+$", b, fixed = TRUE)
  r <- regexpr('remote_fixed "', b, fixed = TRUE)
  expect_true(v > 0 && r > 0 && v < r)
  expect_true(grepl("scontrol show job", b, fixed = TRUE))
})

test_that("remote() enforces the allowlist and rejects metacharacters", {
  skip_if_not(file.exists(ralpine_path()))
  body <- first_group("(?s)\\nremote\\(\\) \\{(.*?)\\n\\}\\n", ralpine_text())
  expect_false(is.na(body))
  expect_true(grepl("is_read_only", body, fixed = TRUE))
  expect_true(grepl("reject_metacharacters", body, fixed = TRUE))
})

test_that("the default root is the organoid-transfer worktree, not the shared checkout", {
  skip_if_not(file.exists(ralpine_path()))
  expect_true(grepl("repositories/GRACKLE-organoid-transfer}", ralpine_text(), fixed = TRUE))
  expect_false(grepl("repositories/GRACKLE}", ralpine_text(), fixed = TRUE))
})
