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

    ## Perform SVD
    svd_result <- svd(Y)
    ## Initialize W and H matrices using SVD
    W <- abs(svd_result$u[,1:k] %*% diag(svd_result$d[1:k]))
    H <- t(abs(svd_result$v[,1:k]))
    W <- min_max_scale(W, min_vals = apply(W,2,min), max_vals = apply(W,2,max)) 
    H <- t(min_max_scale(t(H), min_vals = apply(H,1,min), max_vals = apply(H,1,max))) 
    W <- tf$convert_to_tensor(W, dtype = tf$float16)
    H <- tf$convert_to_tensor(H, dtype = tf$float16)
    
    Y <- tf$constant(Y, dtype = tf$float16)
    

    
    patient_similarity_tf <- tf$convert_to_tensor(patient_similarity, dtype = tf$float16)
    ## net_similarity_tf <- tf$convert_to_tensor(net_similarity, dtype = tf$float16)
    D_p_tf <- tf$linalg$diag(tf$reduce_sum(patient_similarity_tf, axis = as.integer(0)))
    L_p_tf <- D_p_tf - patient_similarity_tf
    
    
    ## calculate the loss value
    reconstruction_loss <- tf$norm(tf$cast(Y - tf$matmul(W, H), dtype = tf$float32), ord = 'euclidean')$numpy() 
    graph_loss <-  lambda_1 * tf$reduce_sum(tf$matmul(W,tf$matmul(L_p_tf,W),transpose_a = T))$numpy()
    prev_obj_value <- reconstruction_loss + graph_loss
    ## print(paste("graph loss:", graph_loss))
    ## print(paste("total_loss:",prev_obj_value))
    
    for(i in 1:iterations){
        ## Update H
        H <- H * (tf$matmul(tf$transpose(W),Y ) / (tf$matmul(tf$transpose(W), tf$matmul(W, H))))
        ## Update W
        W <- W * ((tf$matmul(Y, tf$transpose(H)) + (lambda_1 * tf$matmul(patient_similarity_tf,W))) / (tf$matmul(W, tf$matmul(H, tf$transpose(H))) + (lambda_1 * tf$matmul(D_p_tf,W))))
        ## Normalize W and H
        W <- W / (tf$reduce_sum(W, axis = as.integer(1), keepdims = TRUE) )
        H <- H / (tf$reduce_sum(H, axis = as.integer(0), keepdims = TRUE) )
        ## Check for convergence
        reconstruction_loss <- tf$norm(tf$cast(Y - tf$matmul(W, H), dtype = tf$float32), ord = 'euclidean')$numpy() 
        graph_loss <-  lambda_1 * tf$reduce_sum(tf$matmul(W,tf$matmul(L_p_tf,W),transpose_a = T))$numpy()
        curr_obj_value <- reconstruction_loss + graph_loss
        ##print(paste("graph loss:", graph_loss))
        ##print(paste("total_loss:",curr_obj_value))
            
        if((abs(prev_obj_value - curr_obj_value)/prev_obj_value)  < diff_threshold) {
             break
        }
        ## Update the previous objective value
        prev_obj_value <- curr_obj_value
    }

    W <- W$numpy()
    H <- H$numpy()
    colnames(W) <- paste0("LV",1:k)
    rownames(H) <- paste0("LV",1:k)
    
        
    return(list(W =W, H= H))


    
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


