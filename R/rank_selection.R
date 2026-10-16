topSingularValues <- function(A, nv) {
  if (nv >= min(dim(A)) / 2) return(svd(A, nu = 0, nv = 0)$d[seq_len(nv)])
  irlba::irlba(A, nv = nv)$d
}

#' parallelAnalysisRank
#'
#' Horn's parallel analysis (spec Section 7.1). Null matrices permute each gene (column)
#' independently. k is the number of leading components whose observed singular value exceeds
#' the null `quantile` at the same index. If k reaches `nv`, the analysis reruns with `2 * nv`.
#' @param Y Samples x genes matrix (uncentered, as factorized).
#' @param n_perm Number of null permutations. Default 100.
#' @param quantile Null quantile. Default 0.95.
#' @param nv Number of singular values computed. Default 200.
#' @param seed Random seed. Default 20261015.
#' @return A list with `k`, `observed`, `threshold` and `null` (n_perm x nv).
#' @export
parallelAnalysisRank <- function(Y, n_perm = 100, quantile = 0.95, nv = 200, seed = 20261015) {
  nv <- min(nv, min(dim(Y)) - 1)
  withr::with_seed(seed, {
    observed <- topSingularValues(Y, nv)
    null <- t(vapply(seq_len(n_perm), function(b) topSingularValues(apply(Y, 2, sample), nv), numeric(nv)))
  })
  threshold <- apply(null, 2, stats::quantile, probs = quantile)
  above <- observed > threshold
  k <- if (all(above)) nv else which(!above)[1] - 1L
  if (k == nv && nv < min(dim(Y)) - 1) return(parallelAnalysisRank(Y, n_perm, quantile, 2 * nv, seed))
  list(k = as.integer(k), observed = observed, threshold = unname(threshold), null = null)
}
