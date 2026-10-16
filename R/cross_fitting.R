#' normalizeRows
#' @param H Matrix.
#' @return `H` with each row scaled to unit L2 norm (zero rows unchanged).
#' @export
normalizeRows <- function(H) {
  n <- sqrt(rowSums(H^2))
  n[n == 0] <- 1
  H / n
}

#' fitTask
#'
#' One cross-fitting task (spec Section 8.1). It builds the sample graph from the training
#' samples, scales lambda from alpha, fits GRACKLE, normalizes the rows of H, and projects the
#' held-out samples onto H with NNLS.
#' @param Y Samples x genes matrix; rownames are profile ids.
#' @param samples Sample table whose `profile_id` equals `rownames(Y)`, with a `fold` column.
#' @param gene_graph Genes x genes matrix aligned with `colnames(Y)`.
#' @param k Rank.
#' @param alpha_1,alpha_2 Penalty shares (spec Section 7.2).
#' @param fold Held-out fold; 0 is the reference fit on all samples.
#' @param train_role "all", or "tumor" (transfer test).
#' @param test_role "all", or "organoid".
#' @param iterations,diff_threshold,dtype Passed to [GRACKLE()].
#' @param k_nn Neighbors in the sample graph.
#' @return A list with H, W_test, lambda, n_iter, converged, final_loss, recon_test, train_ids, test_ids.
#' @export
fitTask <- function(Y, samples, gene_graph, k, alpha_1, alpha_2, fold,
                    train_role = c("all", "tumor"), test_role = c("all", "organoid"),
                    iterations = 2000, diff_threshold = 1e-6, dtype = "float32", k_nn = 10) {
  train_role <- match.arg(train_role)
  test_role <- match.arg(test_role)
  if (!identical(rownames(Y), samples$profile_id) || !identical(colnames(Y), rownames(gene_graph))) {
    stop("fitTask: Y, samples and gene_graph must be aligned (rownames(Y) == samples$profile_id, colnames(Y) == rownames(gene_graph))")
  }
  train <- if (fold == 0) rep(TRUE, nrow(samples)) else samples$fold != fold
  test <- if (fold == 0) rep(FALSE, nrow(samples)) else samples$fold == fold
  if (train_role == "tumor") train <- train & samples$role == "tumor"
  if (test_role == "organoid") test <- test & samples$role == "organoid"
  Ytr <- Y[train, , drop = FALSE]
  S <- buildSampleGraph(samples[train, , drop = FALSE], k_nn = k_nn)$S
  init <- svdInit(Ytr, k)
  lam <- lambdaFromAlpha(Ytr, init$W, init$H, S, gene_graph, alpha_1, alpha_2)
  fit <- GRACKLE(Ytr, net_similarity = gene_graph, patient_similarity = S,
                 diff_threshold = diff_threshold, lambda_1 = lam[["lambda_1"]], lambda_2 = lam[["lambda_2"]],
                 k = k, iterations = iterations, svd_W = init$W, svd_H = init$H, dtype = dtype)
  H <- normalizeRows(fit$H)
  W_test <- NULL
  recon_test <- NA_real_
  if (any(test)) {
    Yte <- Y[test, , drop = FALSE]
    W_test <- project_W(Yte, H)
    recon_test <- sum((Yte - W_test %*% H)^2) / length(Yte)
  }
  list(H = H, W_test = W_test, lambda = lam, n_iter = fit$n_iter, converged = fit$converged,
       final_loss = utils::tail(fit$loss_trace, 1), recon_test = recon_test,
       train_ids = samples$profile_id[train], test_ids = samples$profile_id[test])
}
