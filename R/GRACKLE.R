## GRACKLE algorithm

#' GRACKLE
#' 
#'  Function for running Graph Regularization Across Contextual KnowLedgE
#'
#' @useDynLib GRACKLE
#' @param Y Input gene expression data
#' @param net_similarity Similarity matrix based on GRN
#' @param patient_similarity Similarity matrix based on patient metadata
#' @param diff_threshold A numeric value representing the threshold for the difference in W and H matrices between iterations. Default is 1e-4
#' @param lambda_1 A numeric value representing the first regularization parameter. Default is 0.5.
#' @param lambda_2 A numeric value representing the second regularization parameter. Default is 0.5.
#' @param k An integer specifying the number of clusters or components. Default is 5.
#' @param verbose A boolean value as to whether to print loss error or not. Default is FALSE. 
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
    verbose = FALSE) {
  

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
    
                                        # norm_recon_error <- reconstruction_error_vec/max(reconstruction_error_vec)
                                        # norm_pat_sim_error <- pat_sim_error_vec/max(pat_sim_error_vec)
                                        # norm_grn_error <- grn_error_vec/max(grn_error_vec)
    
    H_diff <- 1
    W_diff <- 1
    
    
    
    while((H_diff > diff_threshold) | (W_diff >diff_threshold)) {
        oldW <- W
        oldH <- H
        
        ## Iteratively update matrices

        W <- W *  ((Y %*% t(H))  / (W%*% H %*% t(H) + lambda_1 * D_p %*% W) )
        ##W <- W  * (1 + learning_rate * ((Y %*% t(H) + lambda_1*patient_similarity %*%W)/(W %*% H %*% t(H) + lambda_1 * D_p %*% W + beta_1*W) -1))
        ##W[W<threshold] <- 0
        ## W[W<0] <- .Machine$double.eps
        ##        W <- W/max(W)
        
        H <-H *  ((t(W) %*% Y ) / (t(W) %*% W %*% H + H%*% D_g * lambda_2) )
        
        ##H <- H  * ( 1+  learning_rate *  ((t(W) %*% Y + lambda_2*H%*%net_similarity ) / (t(W) %*% W %*% H + lambda_2 * H %*% D_g + beta_2*H ) -1 ))
        ##H[H<threshold] <-0
        ## H[H<0] <- .Machine$double.eps
        ##       H <- H/max(H)

        ##W <- pmax(W - beta_1,0)
       ## H <- pmax(H - beta_2,0)
        
        reconstruction_error <- sum((Y-W%*%H)^2)
        pat_sim_error <- lambda_1 * sum(diag(t(W) %*% L_p %*%W))
        grn_error <-lambda_2 * sum(diag(H  %*% L_g %*% t(H)))

        if(verbose) {
            message("reconstruction error=", round(reconstruction_error,2), " | patient similarity error=", round(pat_sim_error,2), " | grn error=", round(grn_error,2))
        }
        
        reconstruction_error_vec <- c(reconstruction_error_vec, reconstruction_error)
        pat_sim_error_vec <- c(pat_sim_error_vec, pat_sim_error)
        grn_error_vec <- c(grn_error_vec, grn_error)
        
                                        # norm_recon_error <- reconstruction_error_vec/max(reconstruction_error_vec)
                                        # norm_pat_sim_error <- pat_sim_error_vec/max(pat_sim_error_vec)
                                        # norm_grn_error <- grn_error_vec/max(grn_error_vec)
        
        
                                        #   if(tail(norm_pat_sim_error,n = 1) < tail(norm_grn_error, n = 1)) {
                                        #     lambda_1 <- lambda_1 * (1-learning_rate)
                                        #   } else {
                                        #     
                                        #     lambda_2 <- lambda_2 * (1-learning_rate)
                                        # }
                                        #   print(paste("lambda_1", lambda_1))
                                        #   print(paste("lambda_2", lambda_2))
        
        W_diff <- round(sum((W - oldW)^2)/sum(W^2),7)
        W_diff_vec <-c(W_diff_vec, W_diff)

        H_diff <- round(sum((H - oldH)^2)/sum(H^2),7)
        H_diff_vec <-c(H_diff_vec, H_diff)
        
    }
    
    colnames(W) <- paste0("LV",1:k)
    rownames(H) <- paste0("LV",1:k)
    
    
    out <- list(residual=(Y-W%*%H), H=H, W=W,error = reconstruction_error_vec,H_diff = H_diff_vec, W_diff= W_diff_vec, pat_sim_error_vec = pat_sim_error_vec, grn_error_vec=grn_error_vec )
    
    
    return(out)
    
}

