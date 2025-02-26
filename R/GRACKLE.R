## GRACKLE algorithm

#' GRACKLE
#' 
#'  Function for running Graph Regularization Across Contextual KnowLedgE
#'
#' @importFrom Matrix Matrix
#' @import dynutils
#' @import tensorflow
#' @param Y Input gene expression data
#' @param net_similarity Similarity matrix based on GRN
#' @param patient_similarity Similarity matrix based on patient metadata
#' @param diff_threshold A numeric value representing the threshold for the difference in W and H matrices between iterations. Default is 1e-4
#' @param lambda_1 A numeric value representing the first regularization parameter. Default is 0.5.
#' @param lambda_2 A numeric value representing the second regularization parameter. Default is 0.5.
#' @param k An integer specifying the number of clusters or components. Default is 5.
#' @param verbose A boolean value as to whether to print loss error or not. Default is FALSE. 
#' @param beta A numeric value representing the degree of l2 regularization. Default is 0.
#' @param error_terms A boolean to determine whether or not to calculate error terms. Default is False
#' @param iterations Number of iterations
#' @return A list containing the following components:
#' \item{residual}{A matrix representing the residuals, calculated as \code{Y - W \%*\% H}.}
#' \item{H}{The matrix \code{H} obtained from the NMF calculation.}
#' \item{W}{The matrix \code{W} obtained from the NMF calculation.}
#' \item{error}{A vector containing the reconstruction error at each iteration.}
#' \item{H_diff}{A vector containing the differences in \code{H} between iterations.}
#' \item{W_diff}{A vector containing the differences in \code{W} between iterations.}
#' \item{pat_sim_error_vec}{A vector containing the patient similarity error at each iteration.}
#' \item{grn_error_vec}{A vector containing the gene regulatory network error at each iteration.}
#' @export
GRACKLE <- function(
    Y,
    net_similarity,
    patient_similarity,
    diff_threshold = 1e-4,
    lambda_1 = .5,
    lambda_2 = .5,
    k = 5, 
    verbose = FALSE,
    beta = 0,
    error_terms = FALSE,
    iterations = 100) {

                                        # Ensure TensorFlow uses the GPU
    Sys.setenv(TF_ENABLE_ONEDNN_OPTS = "0")
    Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")
    physical_devices <- tf$config$list_physical_devices('GPU')
        if (length(physical_devices) > 0) {
        suppressMessages({
            suppressWarnings({
                tf$config$experimental$set_memory_growth(physical_devices[[1]], TRUE)
                tf$config$set_visible_devices(physical_devices[[1]], 'GPU')
            })
        })
        } else {
      print("No GPU found. Using CPU.")
    }

    ## Enable mixed precision
    policy <- tf$keras$mixed_precision$Policy('mixed_float16')
    tf$keras$mixed_precision$set_global_policy(policy)

    ## calc degree
    ##D_p <- diag(rowSums(patient_similarity)) + diag(colSums(patient_similarity))
    ##D_g <- diag(rowSums(net_similarity)) + diag(colSums(net_similarity))
    
    
    ## Initialize matrices
    ## random
    n <- nrow(Y)
    m <- ncol(Y)
    set.seed(42)
    W <- matrix(runif(n * k, min = min(Y), max = max(Y) ), nrow = n, ncol = k)
    set.seed(42)
    H <- matrix(runif(k * m, min = min(Y), max = max(Y)), nrow = k, ncol = m)

    ## Perform SVD
    ## svd_result <- svd(Y)
    ## ## Initialize W and H matrices using SVD
    ## W <- abs(svd_result$u[,1:k] %*% diag(svd_result$d[1:k]))
    ## H <- t(abs(svd_result$v[,1:k]))

    Y_tf <- tf$convert_to_tensor(Y, dtype = tf$float16)
    patient_similarity_tf <- tf$convert_to_tensor(patient_similarity, dtype = tf$float16)
    net_similarity_tf <- tf$convert_to_tensor(net_similarity, dtype = tf$float16)
    W_tf <- tf$convert_to_tensor(W, dtype = tf$float16)
    H_tf <- tf$convert_to_tensor(H, dtype = tf$float16)
    D_p_tf <- tf$linalg$diag(tf$reduce_sum(patient_similarity_tf, axis = as.integer(0)))
    D_g_tf <- tf$linalg$diag(tf$reduce_sum(net_similarity_tf, axis = as.integer(0)))
    
    ## normalize the Laplacians
    ## Calculate the inverse square root of D_p
    ## L_p_norm <- tf$matmul(tf$matmul(D_p_inv_sqrt, L_p), D_p_inv_sqrt)
    ## L_g_norm <- tf$matmul(tf$matmul(D_g_inv_sqrt, L_g), D_g_inv_sqrt)

    ## Normalize the degree matrix
    ## D_p_inv_sqrt <- tf$linalg$diag(tf$sqrt(tf$linalg$diag_part(D_p))^-1)
    ## D_g_inv_sqrt <- tf$linalg$diag(tf$sqrt(tf$linalg$diag_part(D_g))^-1)
    
##   for(i in 1:iterations) {
##        print(i)
     H_diff <- 1        
  while(H_diff > diff_threshold) {
##        ## oldH <- H
     oldH_tf <- H_tf
      ## Iteratively update matrices
      ## base code ######
      ## W <- W *  ((tcrossprod(Y,H) + lambda_1 * patient_similarity %*% W   )  / (W%*% tcrossprod(H, H) + lambda_1 * D_p %*%W))
      ## W <- apply(W,2,function(x) scale(x,center = F))
      ## H <-H *  ((crossprod(W ,Y) +  H %*% net_similarity*lambda_2) / (crossprod(W,W) %*% H + H%*% D_g * lambda_2 ) )
      ## H <- t(apply(H,1,function(x) scale(x,center = F)))
      ## tensor  #######
      ##      print(colSums(as.matrix(W_tf)))
      ## print(rowSums(as.matrix(H_tf)))
      ## Numerator and Denominator for W update
      numerator_W <- tf$add(tf$matmul(Y_tf, H_tf, transpose_b = TRUE), lambda_1 * tf$matmul(patient_similarity_tf, W_tf))
      denominator_W <- tf$add(tf$matmul(W_tf, tf$matmul(H_tf, H_tf, transpose_b = TRUE)), lambda_1 * tf$matmul(D_p_tf, W_tf))
      W_tf <- tf$multiply(W_tf, tf$divide(numerator_W, denominator_W))
      numerator_H <- tf$add(tf$matmul(W_tf, Y_tf, transpose_a = TRUE), tf$matmul(H_tf, net_similarity_tf) * lambda_2)
      denominator_H <- tf$add(tf$matmul(tf$matmul(W_tf, W_tf, transpose_a = TRUE), H_tf), tf$matmul(H_tf, D_g_tf) * lambda_2)
      H_tf <- tf$multiply(H_tf, tf$divide(numerator_H, denominator_H))
      ## H_tf_transposed <- tf$transpose(H_tf)
      ## H_tf_scaled_transposed <- tf$map_fn(scale_without_centering, H_tf_transposed, dtype = tf$float64)
      ## H_tf <- tf$transpose(H_tf_scaled_transposed)
      ## Calculate the difference
      diff <- tf$subtract(H_tf, oldH_tf)
      ## Calculate the numerator: sum((H - oldH)^2)
      numerator <- tf$reduce_sum(tf$square(diff))
      ## Calculate the denominator: sum(H^2)
      denominator <- tf$reduce_sum(tf$square(H_tf))
      ## Calculate _diff
      H_diff <- tf$math$divide(numerator, denominator)
      H_diff <- as.numeric(H_diff)
#      print(H_diff)
    }
    
    W <- as.matrix(W_tf)
    H <- as.matrix(H_tf)


    colnames(W) <- paste0("LV",1:k)
    rownames(H) <- paste0("LV",1:k)
    
    if(error_terms) {
        out <- list(residual=(Y-W%*%H), H=H, W=W,error = reconstruction_error_vec,H_diff = H_diff_vec, W_diff= W_diff_vec, pat_sim_error_vec = pat_sim_error_vec, grn_error_vec=grn_error_vec )
    } else {
        out <- list( H=H, W=W)
    } 
    
    return(out)

    ## if(error_terms) {
    ##             ## Error vectors
    ##     error_vec <- numeric()
    ##     reconstruction_error_vec <- numeric()
    ##     pat_sim_error_vec <- numeric()
    ##     grn_error_vec <- numeric()
    ##             reconstruction_error <- sum((Y-W%*%H)^2)
    ##     pat_sim_error <- lambda_1 * sum(diag(crossprod(W,L_p)%*%W))
    ##     grn_error <-lambda_2 * sum(diag(tcrossprod(H %*% L_g,H)))
    ##     reconstruction_error_vec <- c(reconstruction_error_vec, reconstruction_error)
    ##     pat_sim_error_vec <- c(pat_sim_error_vec, pat_sim_error)
    ##     grn_error_vec <- c(grn_error_vec, grn_error)
    ##     if(verbose) {
    ##         message("reconstruction error=", round(reconstruction_error,2), " | patient similarity error=", round(pat_sim_error,2), " | grn error=", round(grn_error,2))
    ##     }
    ## }

    
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


