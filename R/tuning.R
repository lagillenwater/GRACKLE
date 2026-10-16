#' traceLaplacian
#'
#' `tr(X' L X)` (by = "row", X is samples x k) or `tr(X L X')` (by = "col", X is k x genes)
#' for `L = diag(rowSums(S)) - S`, without forming `L`.
#' @param X Factor matrix.
#' @param S Symmetric similarity matrix.
#' @param by "row" or "col".
#' @return Numeric scalar.
#' @export
traceLaplacian <- function(X, S, by = c("row", "col")) {
  by <- match.arg(by)
  d <- rowSums(S)
  if (by == "row") sum(d * rowSums(X^2)) - sum(X * (S %*% X)) else sum(d * colSums(X^2)) - sum(X * (X %*% S))
}

#' lambdaFromAlpha
#'
#' Dimensionless penalty weights (spec Section 7.2): each penalty equals `alpha` times the
#' reconstruction error at the starting point. A graph with no edges gives lambda = 0.
#' @param Y Training matrix.
#' @param W0,H0 Starting factors, for example from [svdInit()].
#' @param S_p Sample similarity.
#' @param S_g Gene similarity.
#' @param alpha_1,alpha_2 Penalty shares.
#' @return Named numeric vector with `lambda_1` and `lambda_2`.
#' @export
lambdaFromAlpha <- function(Y, W0, H0, S_p, S_g, alpha_1, alpha_2) {
  recon <- sum((Y - W0 %*% H0)^2)
  g1 <- traceLaplacian(W0, S_p, "row")
  g2 <- traceLaplacian(H0, S_g, "col")
  c(lambda_1 = if (g1 > 0) alpha_1 * recon / g1 else 0,
    lambda_2 = if (g2 > 0) alpha_2 * recon / g2 else 0)
}

#' selectSetting
#'
#' Spec Section 7.3. Settings with a non-converged or missing fold are excluded. Among
#' settings with held-out reconstruction within one SE of the best, the most reproducible
#' wins; ties go to the smaller alpha_1, then the smaller alpha_2.
#' @param tab Data frame with setting_id, alpha_1, alpha_2, all_converged, recon_mean, recon_se, reproducibility.
#' @return The selected setting_id.
#' @export
selectSetting <- function(tab) {
  ok <- tab[tab$all_converged %in% TRUE, , drop = FALSE]
  if (!nrow(ok)) stop("selectSetting: no setting has all folds present and converged")
  best <- which.min(ok$recon_mean)
  band <- ok[ok$recon_mean <= ok$recon_mean[best] + ok$recon_se[best], , drop = FALSE]
  band <- band[order(-band$reproducibility, band$alpha_1, band$alpha_2), , drop = FALSE]
  band$setting_id[1]
}
