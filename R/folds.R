#' makeFolds
#'
#' Assigns every case to one of `n_folds` cross-validation folds, stratified by cancer type.
#' Cases are shuffled, sorted by cancer type, and dealt round-robin, so fold sizes within a
#' type differ by at most one. All profiles of a case share its fold.
#'
#' @param samples Data frame with columns `case_id` and `cancer_type`.
#' @param n_folds Number of folds. Default 5.
#' @param seed Random seed. Default 20261015. The caller's RNG state is restored.
#' @return `samples` with an integer `fold` column.
#' @export
makeFolds <- function(samples, n_folds = 5, seed = 20261015) {
  cases <- unique(samples[, c("case_id", "cancer_type")])
  if (anyDuplicated(cases$case_id)) {
    stop("makeFolds: a case has more than one cancer_type: ",
         paste(unique(cases$case_id[duplicated(cases$case_id)]), collapse = ", "))
  }
  cases <- withr::with_seed(seed, cases[sample.int(nrow(cases)), , drop = FALSE])
  cases <- cases[order(cases$cancer_type, method = "radix"), , drop = FALSE]
  cases$fold <- ((seq_len(nrow(cases)) - 1L) %% n_folds) + 1L
  samples$fold <- cases$fold[match(samples$case_id, cases$case_id)]
  samples
}
