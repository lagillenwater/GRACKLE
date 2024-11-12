# Similarity Functions
#' similarity_calc
#' 
#'  This function calculates the similarity matrix from an input matrix
#'
#' @useDynLib GRACKLE
#' @import Rcpp
#' @param input_matrix Input matrix of adjacency matrix or 
#' @export
similarity_calc <- function(input_matrix) {
    if (!is.matrix(input_matrix)) {
        stop("Error: Input is not a matrix.")
          }
    
    distance <- as.matrix(dist(input_matrix, method = "euclidean"))
 ##   similarity <- 1 / (1 + as.matrix(distance))
    similarity <- exp(-distance)
#    similarity[similarity < threshold] <- 0
    colnames(similarity) <- rownames(input_matrix)
    rownames(similarity) <- rownames(input_matrix)
    return(similarity)
}

