## Evaluation functions

#' categoricalPredication
#' 
#'  Function for predicting categorical labels in test data
#'
#' @import glmnet
#' @param W_test Matrix of projected W from test expression
#' @param test_metadata Matrix of test metadata
#' @return nonzero Named list of nonzero coefficients by subgroup prediction
#' @export
categoricalPrediction <- function(W_test,test_metadata){
    
    test_metadata <- as.data.frame(test_metadata)
        
    ## Fit penalized logistic regression models for each outcome
    models <- lapply(colnames(test_metadata), function(outcome) {
        y <- as.matrix(test_metadata[[outcome]])
        cv.glmnet(W_test, y,  alpha = 1, intercept = F, lower.limit = 0, stadardize = T)  # alpha = 1 for Lasso, alpha = 0 for Ridgee
    })

    
    ## coefficients
    coefficients <- lapply(models, function(model) {
        coef(model, s = model$lambda.min)
        
    })

    ## nonzero coefficients
    nonzero <- lapply(coefficients, function(x) {
        rownames(x)[x[,1] != 0]
        
    })
    names(nonzero) <- names(test_metadata)

    ## top coefficient
    top <- lapply(coefficients, function(x) {
        rownames(x)[x[,1] == max(x)]
    })

    names(top) <- names(test_metadata)

    return(list(nonzero = nonzero, top = top))
    
    
}


#' geneLoadingsEvaluation
#' 
#'  geneLoadingsEvaluation is a function for identifying the top loading genes on each latent variable and aligning them with input data
#'
#' @param H_train Matrix of trained H matrix
#' @param module_size Size of consecutive modules for DREAM simulations.
#' @export

geneLoadingsEvaluation <- function(H_train, module_size = 10) {
    n_genes <- ncol(H_train)
    lower_bound <- seq(1,n_genes,module_size)
    upper_bound <- lower_bound + module_size -1

    module_scores <- lapply(1:length(lower_bound), function(i) {
        rowSums(H_train[,lower_bound[i]:upper_bound[i]])
    })
    
    return(module_scores)

}
