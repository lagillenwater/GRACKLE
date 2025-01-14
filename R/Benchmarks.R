## Functions for benchmarking algorithms. 


#' runNMF
#' 
#' runNMF runs the version of NMF from the NMF package. 
#'
#' @import NMF
#' @param expression_matrix Matrix of gene expression data (samples as rows, genes as columns)
#' @param k A numeric variable representing the number of latent variables (aka, rank). Default is 5.
#' @param method A string representing the nmf algorithm to run. See nmf for options. Default is 'lee'.
#' @param seed Specifies the starting point or seeding method. See the NMF package for details.  Default is 42.
#' @return A list containing the following components:
#'
#' \item{W}{The matrix \code{W} obtained from the NMF calculation.}
#' @export
runNMF <- function(expression_matrix,k = 5, method = 'lee', seed = 42) {
 
  
  # Run NMF with parallel computing
  res <- NMF::nmf(expression_matrix, rank = k, seed = seed)
  
    
    W <- basis(res)
    H <- coef(res)
    
  return(list(W = W,H = H))
} 

#' runGRNMF
#' 
#' runGRNMF implements a version of graph regularized NMF
#'
#' @importFrom FNN get.knn
#' @param expression_matrix Matrix of gene expression data (samples as rows, genes as columns)
#' @param k A numeric variable representing the number of latent variables (aka, rank). Default is 5.
#' @param nearest_neighbors A numeric value representing the number of nearest neighbors. (Default is 5. )
#' @param seed Specifies the starting point or seeding method. See the NMF package for details.  (Default is 42.)
#' @param max_iter Numeric for the maximum number of iterations. (Default is 200)
#' @param alpha Numeric value representing the effect of graph regularization on the model. (Default is .1)
#' @return
#' @export

runGRNMF <- function(expression_matrix, k = 5,  seed = 42, max_iter = 200, alpha = .1){
  set.seed(seed)
  # Create an affinity graph using k-nearest neighbors
  knn_graph <- get.knn(t(expression_matrix), k = k)$nn.index
  W_graph <- matrix(0, nrow = ncol(expression_matrix), ncol = ncol(expression_matrix))
  for (i in 1:nrow(knn_graph)) {
    W_graph[i, knn_graph[i, ]] <- 1
    W_graph[knn_graph[i, ], i] <- 1
  }
  
  # Initialize W and H matrices
    W <- matrix(runif(nrow(expression_matrix) * k), nrow = nrow(expression_matrix), ncol = k)
    H <- matrix(runif(k * ncol(expression_matrix)), nrow = k, ncol = ncol(expression_matrix))
  
  
  # GNMF update function
  for (i in 1:max_iter) {
      
      W <- W * (expression_matrix %*% t(H)) / (W %*% H %*% t(H))
      W <- apply(W,2,function(x) scale(x,center = F))      
      H <- H * (t(W) %*% expression_matrix) / (t(W) %*% W %*% H + alpha * H %*% W_graph)
       H <- t(apply(H,1,function(x) scale(x,center = F))) 

  }
  
    # Extract the basis and coefficient matrices
  return(list(W =  W, H = H))
  
} 

#' runPLIER
#' 
#' runPLIER implements a version of PLIER for comparison. 
#'
#' @import PLIER
#' @param expression_matrix Matrix of gene expression data (samples as rows, genes as columns)
#' @param k A numeric variable representing the number of latent variables (aka, rank). Default is 5.
#' @return
#' @export
runPLIER <- function(expression_matrix, k) {
  data(canonicalPathways)
  res <- PLIER(data = expression_matrix,priorMat = canonicalPathways,k,rseed = 42)
  W= res$B
  H = res$Z
  return(list(W =  W, H = H))
}
