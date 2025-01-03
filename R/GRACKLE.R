## GRACKLE algorithm

#' GRACKLE
#' 
#'  Function for running Graph Regularization Across Contextual KnowLedgE
#'
#' @param Y Input gene expression data
#' @param net_similarity Similarity matrix based on GRN
#' @param patient_similarity Similarity matrix based on patient metadata
#' @param diff_threshold A numeric value representing the threshold for the difference in W and H matrices between iterations. Default is 1e-4
#' @param lambda_1 A numeric value representing the first regularization parameter. Default is 0.5.
#' @param lambda_2 A numeric value representing the second regularization parameter. Default is 0.5.
#' @param k An integer specifying the number of clusters or components. Default is 5.
#' @param verbose A boolean value as to whether to print loss error or not. Default is FALSE. 
#' @param beta A numeric value representing the degree of l2 regularization. Default is 0. 
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
    beta = 0) {
  
    set.seed(42)
    ## calc degree
    D_p <- diag(rowSums(patient_similarity))
    D_g <- diag(rowSums(net_similarity)) 
    ## calc graph laplacian
    L_p <- D_p - patient_similarity
    L_g <- D_g - net_similarity

    ## Initialize matrices
    n <- nrow(Y)
    m <- ncol(Y)
    W <- matrix(runif(n * k, min = 0, max = 1), nrow = n, ncol = k)
    H <- matrix(runif(k * m, min = 0, max = 1), nrow = k, ncol = m)
    
    ## Error vectors
    error_vec <- numeric()
    reconstruction_error_vec <- numeric()
    pat_sim_error_vec <- numeric()
    grn_error_vec <- numeric()
  
    H_diff_vec <- numeric()
    W_diff_vec <- numeric()
    
     
    reconstruction_error <- sum((Y-W%*%H)^2)
    pat_sim_error <- lambda_1 * sum(diag(t(W) %*% L_p %*%W))
    grn_error <-lambda_2 * sum(diag(H  %*% L_g %*% t(H)))
    
    if(verbose) {
        message("reconstruction error=", round(reconstruction_error,2), " | patient similarity error=", round(pat_sim_error,2), " | grn error=", round(grn_error,2))
    }
    
    reconstruction_error_vec <- c(reconstruction_error_vec, reconstruction_error)
    pat_sim_error_vec <- c(pat_sim_error_vec, pat_sim_error)
    grn_error_vec <- c(grn_error_vec, grn_error)
        
    H_diff <- 1
    W_diff <- 1
    
        
    while((H_diff > diff_threshold) | (W_diff >diff_threshold)) {
        oldW <- W
        oldH <- H
                ## Iteratively update matrices
        W <- W *  ((Y %*% t(H) + lambda_1 * patient_similarity %*% W  )   / (W%*% H %*% t(H) + lambda_1 * D_p %*% W ) )
        W <- apply(W,2,function(x) scale(x,center = F))
        H <-H *  ((t(W) %*% Y + H %*% net_similarity * lambda_2) / (t(W) %*% W %*% H + H%*% D_g * lambda_2 ) )
        H <- t(apply(H,1,function(x) scale(x,center = F)))
        reconstruction_error <- sum((Y-W%*%H)^2)
        pat_sim_error <- lambda_1 * sum(diag(t(W) %*% L_p %*%W))
        grn_error <-lambda_2 * sum(diag(H  %*% L_g %*% t(H)))
        if(verbose) {
            message("reconstruction error=", round(reconstruction_error,2), " | patient similarity error=", round(pat_sim_error,2), " | grn error=", round(grn_error,2))
        }
        reconstruction_error_vec <- c(reconstruction_error_vec, reconstruction_error)
        pat_sim_error_vec <- c(pat_sim_error_vec, pat_sim_error)
        grn_error_vec <- c(grn_error_vec, grn_error)
        W_diff <- round(sum((W - oldW)^2)/sum(W^2),7)
      #  print(W_diff)
        W_diff_vec <-c(W_diff_vec, W_diff)
        H_diff <- round(sum((H - oldH)^2)/sum(H^2),7)
       # print(H_diff)
        H_diff_vec <-c(H_diff_vec, H_diff)
         #       if(reconstruction_error > (1.25 * reconstruction_error_vec[length(reconstruction_error_vec) - 1])) {break}
    }
    
    colnames(W) <- paste0("LV",1:k)
    rownames(H) <- paste0("LV",1:k)
    
    
    out <- list(residual=(Y-W%*%H), H=H, W=W,error = reconstruction_error_vec,H_diff = H_diff_vec, W_diff= W_diff_vec, pat_sim_error_vec = pat_sim_error_vec, grn_error_vec=grn_error_vec )
    
    
    return(out)
    
}

