#' cosineSim
#' @param A,B Matrices with the same number of columns.
#' @return nrow(A) x nrow(B) matrix of cosine similarities between rows.
#' @export
cosineSim <- function(A, B) {
  na <- sqrt(rowSums(A^2)); nb <- sqrt(rowSums(B^2))
  na[na == 0] <- 1; nb[nb == 0] <- 1
  (A %*% t(B)) / outer(na, nb)
}

#' assignMax
#'
#' One-to-one assignment maximizing total similarity (Hungarian algorithm via `clue`).
#' @param sim Square nonnegative similarity matrix.
#' @return Integer vector: the column assigned to each row.
#' @export
assignMax <- function(sim) {
  sim[is.na(sim) | sim < 0] <- 0
  as.integer(clue::solve_LSAP(sim, maximum = TRUE))
}

#' matchPrograms
#'
#' Matches a fold's programs to the reference programs by gene-loading cosine (spec Section 8.2).
#' @param H_ref Reference H (k x genes).
#' @param H_fold Fold H (k x genes).
#' @param tau Cosine below which a forced match counts as unmatched.
#' @return Data frame with ref, fold_prog, cosine, matched.
#' @export
matchPrograms <- function(H_ref, H_fold, tau) {
  sim <- cosineSim(H_ref, H_fold)
  a <- assignMax(sim)
  cos <- sim[cbind(seq_len(nrow(sim)), a)]
  data.frame(ref = seq_len(nrow(sim)), fold_prog = a, cosine = cos, matched = cos >= tau)
}

#' matchNull
#'
#' Best-match cosine of each fold program against all reference programs after permuting the
#' fold's gene loadings. Pool across folds and take the 99th percentile for `tau`.
#' @param H_ref,H_fold k x genes matrices.
#' @param n_perm Number of permutations.
#' @param seed Random seed.
#' @return Numeric vector of length nrow(H_fold) * n_perm.
#' @export
matchNull <- function(H_ref, H_fold, n_perm = 1000, seed = 1) {
  withr::with_seed(seed, as.vector(vapply(seq_len(n_perm), function(b) {
    Hp <- H_fold[, sample.int(ncol(H_fold)), drop = FALSE]
    apply(cosineSim(Hp, H_ref), 1, max)
  }, numeric(nrow(H_fold)))))
}

#' programStability
#' @param matches Data frame with fold, ref, cosine, matched (one row per fold x reference program).
#' @param min_folds Folds a program must match in to count as stable. Default 4.
#' @return Data frame with ref, n_matched, mean_cosine, stable.
#' @export
programStability <- function(matches, min_folds = 4) {
  n_matched <- tapply(matches$matched, matches$ref, sum)
  mean_cos <- tapply(matches$cosine, matches$ref, mean)
  data.frame(ref = as.integer(names(n_matched)), n_matched = as.integer(n_matched),
             mean_cosine = as.numeric(mean_cos), stable = as.integer(n_matched) >= min_folds)
}
