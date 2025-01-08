## GRACKLE algorithm

#' GRACKLE
#' 
#'  Function for running Graph Regularization Across Contextual KnowLedgE
#'
#' @importFrom Matrix Matrix
#' @import gpuR
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
    error_terms = FALSE) {

    Y <- Matrix(Y, sparse = T)
    

    ## calc degree
    D_p <- diag(rowSums(patient_similarity))
    D_g <- diag(rowSums(net_similarity))

    D_p <- Matrix(D_p, sparse = T)
    D_g <- Matrix(D_g, sparse = T)
    
    ## calc graph laplacian
    L_p <- D_p - patient_similarity
    L_g <- D_g - net_similarity

    L_p <- Matrix(L_p, sparse = T)
    L_g <- Matrix(L_g, sparse = T)
    
    ## Initialize matrices
    n <- nrow(Y)
    m <- ncol(Y)
    set.seed(42)
    W <- Matrix(runif(n * k, min = 0, max = 1), nrow = n, ncol = k, sparse = T)
    H <- Matrix(runif(k * m, min = 0, max = 1), nrow = k, ncol = m, sparse = T)
    H_diff_vec <- numeric()
    W_diff_vec <- numeric()
    
    if(error_terms) {
        
        ## Error vectors
        error_vec <- numeric()
        reconstruction_error_vec <- numeric()
        pat_sim_error_vec <- numeric()
        grn_error_vec <- numeric()
        
        reconstruction_error <- sum((Y-W%*%H)^2)
        pat_sim_error <- lambda_1 * sum(diag(crossprod(W,L_p)%*%W))
        grn_error <-lambda_2 * sum(diag(tcrossprod(H %*% L_g,H)))
        reconstruction_error_vec <- c(reconstruction_error_vec, reconstruction_error)
        pat_sim_error_vec <- c(pat_sim_error_vec, pat_sim_error)
        grn_error_vec <- c(grn_error_vec, grn_error)
    }
    if(verbose) {
        message("reconstruction error=", round(reconstruction_error,2), " | patient similarity error=", round(pat_sim_error,2), " | grn error=", round(grn_error,2))
    }
    
    H_diff <- 1
    while(H_diff > diff_threshold) {
        oldH <- H
        ## Iteratively update matrices
        W <- W *  ((tcrossprod(Y,H) + lambda_1 * patient_similarity %*% W  )   )  / (W%*% tcrossprod(H, H) + lambda_1 * D_p %*%W)
        ##W <- apply(W,2,function(x) scale(x,center = F))
        ##     W <- min_max_scale(W, min_vals = min(W), max_vals = max(W))
        ## W <- Matrix(W,sparse = T)

        H <-H *  ((crossprod(W ,Y) + H %*% net_similarity * lambda_2) / (crossprod(W,W) %*% H + H%*% D_g * lambda_2 ) )
        ## H <- t(min_max_scale(t(as.matrix(H)), min_vals = min(H), max_vals = max(H)))
        ## H <- t(apply(H,1,function(x) scale(x,center = F)))
        ## H <- Matrix(H, sparse = T)
        if(error_terms) {
            reconstruction_error <- sum((Y-W%*%H)^2)
            pat_sim_error <- lambda_1 * sum(diag(t(W) %*% L_p %*%W))
            grn_error <-lambda_2 * sum(diag(H  %*% L_g %*% t(H)))
            reconstruction_error_vec <- c(reconstruction_error_vec, reconstruction_error)
            pat_sim_error_vec <- c(pat_sim_error_vec, pat_sim_error)
            grn_error_vec <- c(grn_error_vec, grn_error)
        }
        if(verbose) {
            message("reconstruction error=", round(reconstruction_error,2), " | patient similarity error=", round(pat_sim_error,2), " | grn error=", round(grn_error,2))
        }
        H_diff <- sum((H - oldH)^2) / sum(H^2)
    }
    
    
            
    colnames(W) <- paste0("LV",1:k)
    rownames(H) <- paste0("LV",1:k)
    
    if(error_terms) {
        out <- list(residual=(Y-W%*%H), H=H, W=W,error = reconstruction_error_vec,H_diff = H_diff_vec, W_diff= W_diff_vec, pat_sim_error_vec = pat_sim_error_vec, grn_error_vec=grn_error_vec )
    } else {
        out <- list( H=H, W=W)
    } 
    
    return(out)
    
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
