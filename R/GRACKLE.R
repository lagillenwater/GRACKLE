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
    iterations = 100,
    svd_W,
    svd_H) {

    ## Ensure TensorFlow uses the GPU
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

    


    success <- FALSE
    while(!success) {
        tryCatch({
            W <- tf$convert_to_tensor(svd_W, dtype = tf$float16)
            H <- tf$convert_to_tensor(svd_H, dtype = tf$float16)
            Y <- tf$constant(Y, dtype = tf$float16)
            success <- TRUE
        }, error = function(e) {
            cat("Error occurred: ", e$message, "\nRetrying...\n")
        })
    }
    patient_similarity_tf <- tf$convert_to_tensor(patient_similarity, dtype = tf$float16)
    net_similarity_tf <- tf$convert_to_tensor(net_similarity, dtype = tf$float16)
    D_p_tf <- tf$linalg$diag(tf$reduce_sum(patient_similarity_tf, axis= as.integer(0)))
    L_p_tf <- D_p_tf - patient_similarity_tf
    D_g_tf <- tf$linalg$diag(tf$reduce_sum(net_similarity_tf, axis = as.integer(0)))
    L_g_tf= D_g_tf - net_similarity_tf
    ## calculate the loss value
    reconstruction_loss <- tf$norm(tf$cast(Y - tf$matmul(W, H), dtype = tf$float32), ord = 'euclidean')$numpy() 
    pat_graph_loss <-  lambda_1 * tf$reduce_sum(tf$matmul(W,tf$matmul(L_p_tf,W),transpose_a = T))$numpy()
    if(!is.finite(pat_graph_loss)) {
        tmp_W <- tf$cast(W,dtype=tf$float32)
        tmp_L_p_tf <- tf$cast(L_p_tf, dtype = tf$float32)
        pat_graph_loss <-  lambda_1 * tf$reduce_sum(tf$matmul(tmp_W,tf$matmul(tmp_L_p_tf,tmp_W),transpose_a = T))$numpy()
    }
    gene_graph_loss <-  lambda_2 * tf$reduce_sum(tf$matmul(H,tf$matmul(L_g_tf,H, transpose_b = T)))$numpy()
    if(!is.finite(gene_graph_loss)) {
        tmp_H <- tf$cast(H,dtype=tf$float32)
        tmp_L_g_tf <- tf$cast(L_g_tf, dtype = tf$float32)
        gene_graph_loss <-  lambda_2 * tf$reduce_sum(tf$matmul(tmp_H,tf$matmul(tmp_L_g_tf,tmp_H, transpose_b=T)))$numpy()
    }
    prev_obj_value <- gene_graph_loss + pat_graph_loss + reconstruction_loss
    if(verbose){
        print(paste("reconstruction loss:", reconstruction_loss))
        print(paste("patient graph loss:", pat_graph_loss))
        print(paste("gene graph loss:", gene_graph_loss))
        print(paste("total_loss:",prev_obj_value))
    }

    counter = 1
    for(i in 1:292){
        ## Update H
        H <- H * ((tf$matmul(tf$transpose(W),Y ) + (lambda_2 * tf$matmul(H,net_similarity_tf))) / (tf$matmul(tf$transpose(W), tf$matmul(W, H)) + (lambda_2* tf$matmul(H,D_g_tf))))
        ## H_norms <- tf$norm(H, axis = as.integer(0))
        ## H <- H/H_norms
        ## Update W
        W <- W * ((tf$matmul(Y, tf$transpose(H)) + (lambda_1 * tf$matmul(patient_similarity_tf,W))) / (tf$matmul(W, tf$matmul(H, tf$transpose(H))) + (lambda_1 * tf$matmul(D_p_tf,W))))
        ## W_norms <- tf$norm(W, axis = as.integer(0))
        ## W <- W/W_norms
        ## Check for convergence
        success <- FALSE
        while(!success) {
            tryCatch({
                tmp_W <- tf$cast(W,dtype=tf$float32)
                tmp_H <- tf$cast(H,dtype=tf$float32)
                tmp_Y <- tf$cast(Y, dtype=tf$float32)
                reconstruction_loss <- tf$norm(tf$cast(tmp_Y - tf$matmul(tmp_W, tmp_H), dtype = tf$float32), ord = 'euclidean')$numpy()
                success <- TRUE
            }, error = function(e) {
                cat("Error occurred: ", e$message, "\nRetrying...\n")
            })
        }
        pat_graph_loss <-  lambda_1 * tf$reduce_sum(tf$matmul(W,tf$matmul(L_p_tf,W),transpose_a = T))$numpy()
        if(!is.finite(pat_graph_loss)) {
            tmp_W <- tf$cast(W,dtype=tf$float32)
            tmp_L_p_tf <- tf$cast(L_p_tf, dtype = tf$float32)
            pat_graph_loss <-  lambda_1 * tf$reduce_sum(tf$matmul(tmp_W,tf$matmul(tmp_L_p_tf,tmp_W),transpose_a = T))$numpy()
        }
        gene_graph_loss <-  lambda_2 * tf$reduce_sum(tf$matmul(H,tf$matmul(L_g_tf,H, transpose_b = T)))$numpy()
        if(!is.finite(gene_graph_loss)) {
            tmp_H <- tf$cast(H,dtype=tf$float32)
            tmp_L_g_tf <- tf$cast(L_g_tf, dtype = tf$float32)
            gene_graph_loss <-  lambda_2 * tf$reduce_sum(tf$matmul(tmp_H,tf$matmul(tmp_L_g_tf,tmp_H, transpose_b=T)))$numpy()
        }
        curr_obj_value <- gene_graph_loss + pat_graph_loss + reconstruction_loss
        if(verbose) {
            print(paste("reconstruction loss:", reconstruction_loss))
            print(paste("pat graph loss:", pat_graph_loss))
            print(paste("gene graph loss:", gene_graph_loss))
            print(paste("total_loss:",curr_obj_value))
            print(counter)
        }
        if((abs(prev_obj_value - curr_obj_value)/prev_obj_value)  < diff_threshold) {
             break
        }
        ## Update the previous objective value
        prev_obj_value <- curr_obj_value
        ## update counter
        counter = counter+1
    }

    W <- W$numpy()
    H <- H$numpy()
    colnames(W) <- paste0("LV",1:k)
    rownames(H) <- paste0("LV",1:k)
    ##W <- min_max_scale(W, min_vals = apply(W,2,min), max_vals = apply(W,2,max)) 
    ##H <- t(min_max_scale(t(H), min_vals = apply(H,1,min), max_vals = apply(H,1,max)))

    if(verbose) {
        print(paste("reconstruction loss:", reconstruction_loss))
        print(paste("pat graph loss:", pat_graph_loss))
        print(paste("gene graph loss:", gene_graph_loss))
        print(paste("total_loss:",curr_obj_value))
        print(counter)
    }

        
    return(list(W =W, H= H, svd_W = svd_W, svd_H=svd_H))


    
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


