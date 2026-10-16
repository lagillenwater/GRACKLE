## GRACKLE algorithm

#' GRACKLE
#'
#' Graph Regularized non-negative matrix factorization Across Contextual KnowLedgE.
#' Factorizes `Y` (samples x genes) as `W %*% H` with a sample-graph penalty on `W` and a
#' gene-graph penalty on `H`, using multiplicative updates. The objective is
#' `||Y - WH||_F^2 + lambda_1 tr(W' L_p W) + lambda_2 tr(H L_g H')` with `L = D - S` for each
#' similarity matrix `S`. Fitting stops when the relative change in that objective falls
#' below `diff_threshold`, or after `iterations` iterations.
#'
#' @importFrom Matrix Matrix
#' @import tensorflow
#' @param Y Numeric matrix, samples x genes, nonnegative.
#' @param net_similarity Gene similarity (genes x genes), symmetric and nonnegative.
#' @param patient_similarity Sample similarity (samples x samples), symmetric and nonnegative.
#' @param diff_threshold Relative change in the objective that counts as converged. Default 1e-4.
#' @param lambda_1 Sample-graph penalty weight. Default 0.5.
#' @param lambda_2 Gene-graph penalty weight. Default 0.5.
#' @param k Number of latent variables. Default 5.
#' @param verbose Print the objective at each iteration. Default FALSE.
#' @param beta Unused; kept for backward compatibility.
#' @param error_terms Unused; kept for backward compatibility.
#' @param iterations Maximum number of iterations. Default 100.
#' @param svd_W Initial W (samples x k), for example from [svdInit()].
#' @param svd_H Initial H (k x genes), for example from [svdInit()].
#' @param dtype TensorFlow precision: "float32" (default), "float16" or "float64".
#' @return A list with `W`, `H`, `svd_W`, `svd_H`, `loss_trace` (a data frame with columns
#'   iteration, reconstruction, sample_graph, gene_graph, total), `n_iter` and `converged`.
#' @export
GRACKLE <- function(Y, net_similarity, patient_similarity, diff_threshold = 1e-4,
                    lambda_1 = .5, lambda_2 = .5, k = 5, verbose = FALSE, beta = 0,
                    error_terms = FALSE, iterations = 100, svd_W, svd_H,
                    dtype = c("float32", "float16", "float64")) {
  dtype <- match.arg(dtype)
  Sys.setenv(TF_ENABLE_ONEDNN_OPTS = "0", TF_CPP_MIN_LOG_LEVEL = "3")
  gpus <- tf$config$list_physical_devices("GPU")
  if (length(gpus) > 0) tf$config$experimental$set_memory_growth(gpus[[1]], TRUE)

  tf_dtype <- switch(dtype, float32 = tf$float32, float16 = tf$float16, float64 = tf$float64)
  as_tf <- function(x) tf$constant(unname(as.matrix(x)), dtype = tf_dtype)
  ## Bookkeeping in at least float32 so float16 sums cannot overflow.
  wide <- function(x) if (dtype == "float16") tf$cast(x, tf$float32) else x

  Y_tf <- as_tf(Y)
  W <- as_tf(svd_W)
  H <- as_tf(svd_H)
  S_p <- as_tf(patient_similarity)
  S_g <- as_tf(net_similarity)
  d_p <- tf$reshape(tf$reduce_sum(S_p, axis = 1L), c(-1L, 1L))  # sample degrees, n x 1
  d_g <- tf$reshape(tf$reduce_sum(S_g, axis = 0L), c(1L, -1L))  # gene degrees, 1 x g
  l1 <- tf$constant(lambda_1, dtype = tf_dtype)
  l2 <- tf$constant(lambda_2, dtype = tf_dtype)
  ## Guards 0/0 in the multiplicative updates; far below any real denominator.
  eps <- tf$constant(if (dtype == "float16") 1e-7 else 1e-10, dtype = tf_dtype)

  objective <- function(W, H) {
    recon <- tf$reduce_sum(tf$square(wide(Y_tf - tf$matmul(W, H))))
    LpW <- wide(d_p * W - tf$matmul(S_p, W))   # L_p W without forming L_p
    HLg <- wide(H * d_g - tf$matmul(H, S_g))   # H L_g without forming L_g
    parts <- c(reconstruction = as.numeric(recon$numpy()),
               sample_graph = lambda_1 * as.numeric(tf$reduce_sum(wide(W) * LpW)$numpy()),
               gene_graph = lambda_2 * as.numeric(tf$reduce_sum(wide(H) * HLg)$numpy()))
    c(parts, total = sum(parts))
  }

  trace <- matrix(NA_real_, nrow = iterations, ncol = 4,
                  dimnames = list(NULL, c("reconstruction", "sample_graph", "gene_graph", "total")))
  prev <- objective(W, H)
  if (!is.finite(prev[["total"]])) stop("GRACKLE: the initial objective is not finite")
  converged <- FALSE
  n_iter <- 0L
  for (i in seq_len(iterations)) {
    H <- H * ((tf$matmul(W, Y_tf, transpose_a = TRUE) + l2 * tf$matmul(H, S_g)) /
              (tf$matmul(tf$matmul(W, W, transpose_a = TRUE), H) + l2 * (H * d_g) + eps))
    W <- W * ((tf$matmul(Y_tf, H, transpose_b = TRUE) + l1 * tf$matmul(S_p, W)) /
              (tf$matmul(W, tf$matmul(H, H, transpose_b = TRUE)) + l1 * (d_p * W) + eps))
    curr <- objective(W, H)
    if (!is.finite(curr[["total"]])) stop("GRACKLE: the objective became non-finite at iteration ", i)
    trace[i, ] <- curr
    n_iter <- i
    if (verbose) message(sprintf("iteration %d: total %.6g", i, curr[["total"]]))
    if (abs(prev[["total"]] - curr[["total"]]) / prev[["total"]] < diff_threshold) {
      converged <- TRUE
      break
    }
    prev <- curr
  }

  W <- W$numpy()
  H <- H$numpy()
  lv <- paste0("LV", seq_len(k))
  dimnames(W) <- list(rownames(Y), lv)
  dimnames(H) <- list(lv, colnames(Y))
  loss_trace <- data.frame(iteration = seq_len(n_iter), trace[seq_len(n_iter), , drop = FALSE])
  list(W = W, H = H, svd_W = svd_W, svd_H = svd_H, loss_trace = loss_trace,
       n_iter = n_iter, converged = converged)
}

#' svdInit
#'
#' SVD-based nonnegative starting values for [GRACKLE()], matching the HTP analysis scripts:
#' absolute values of the leading singular vectors (W scaled by the singular values), each
#' factor min-max scaled.
#'
#' @param Y Numeric matrix, samples x genes.
#' @param k Number of latent variables.
#' @return A list with `W` (samples x k) and `H` (k x genes).
#' @export
svdInit <- function(Y, k) {
  s <- svd(Y, nu = k, nv = k)
  W <- abs(s$u[, seq_len(k), drop = FALSE] %*% diag(s$d[seq_len(k)], k))
  H <- t(abs(s$v[, seq_len(k), drop = FALSE]))
  W <- min_max_scale(W, apply(W, 2, min), apply(W, 2, max))
  H <- t(min_max_scale(t(H), apply(H, 1, min), apply(H, 1, max)))
  list(W = W, H = H)
}


#' colScale
#' 
#' @description colScale is a fast function for scaling the columns of a matrix. 
#'
#' @import matrixStats
#' @param x Matrix to scale
#' @param center Boolean for whether to center or not. (Default is FALSE)
#' @param  Boolean  for whether to scale or not. (Default is TRUE)
#' @return A list containing the following components:

colScale = function(x,
                    center = FALSE,
                    scale = TRUE) {
                                        
                        
        cm = colMeans(x, na.rm = TRUE)

    if (scale) {
        csd = colSds(x, center = cm)
    } else {
        csd = rep(1, length = length(cm))
    }
    if (!center) {

        cm = rep(0, length = length(cm))
    }
    x = t( (t(x) - cm) / csd )

    
    
    return(x)
}


## #'matMult
## #' @description matMult is equivalent to %*% using C++.
## #' @param A Matrix
## #' @param B Matrix
## #' @export
## #' 
## matMult <- function(A,B) {
##     .Call('_GRACKLE_matMult', PACKAGE = 'GRACKLE', A, B)
## }


