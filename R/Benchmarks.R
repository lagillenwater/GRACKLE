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
#' \item{H}{The matrix \code{H} obtained from the NMF calculation.}
#' \item{W}{The matrix \code{W} obtained from the NMF calculation.}
#' @export
runNMF <- function(expression_matrix,k = 5, method = 'lee', seed = 42) {
  res <- nmf(expression_matrix,rank = k,method = method, seed = seed)
  W <- basis(res)
  H <- coef(res)
  return(list(W = W,H = H))
} 