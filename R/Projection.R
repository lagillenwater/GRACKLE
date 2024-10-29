## Projecting data into the factorized latent space

#' project_W
#' 
#'  Function for predicting categorical labels in test data
#'
#' @useDynLib GRACKLE
#' @import nnls
#' @param test_expression Data.frame of test expression data
#' @param H_train Matrix of trained H matrix
#' @return W_test Matrix of project W values for the test data
#' @export
project_W <- function(test_data,H_train,k) {

    k <- nrow(H_train)
    ## Initialize W_test
    W_test <- matrix(0, nrow=nrow(test_data), ncol=k)

    ## convert test_data to a matrix
    test_data <- as.matrix(test_data)
  
    ## Solve for W_test using NNLS
    for (i in 1:nrow(test_data)) {
        nnls_result <- nnls(t(H_train), test_data[i, ])
        W_test[i, ] <- coef(nnls_result)
    }

    rownames(W_test) <- rownames(test_data)
    colnames(W_test) <- rownames(H_train)
  
  return(W_test)
}
