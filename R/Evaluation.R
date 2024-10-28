## Evaluation functions

#' categoricalPredication
#' 
#'  Function for predicting categorical labels in test data
#'
#' @import glmnet
#' @param W_test Matrix of projected W from test expression
#' @param test_metadata Matrix of test metadata
#' @export
categoricalPrediction <- function(W_test,test_metadata){

    ?glmnet()
    ## W_test <- as.data.frame(W_test)
    test_metadata <- as.data.frame(test_metadata)

    ## Fit penalized logistic regression models for each outcome
    models <- lapply(colnames(test_metadata), function(outcome) {
        y <- as.matrix(test_metadata[[outcome]])
        glmnet(W_test, y, family = "binomial", alpha = 1)  # alpha = 1 for Lasso, alpha = 0 for Ridge
        })
    ## Fit logistic regression models for each outcome
    ## models <- lapply(names(test_metadata), function(outcome) {
    ##     glm(as.formula(paste(outcome, "~ .")), data = cbind(test_metadata[outcome], W_test), family = binomial)
    ## })
    
    ## ## Predict probabilities for each outcome
    ## predicted_probabilities <- lapply(models, function(model) {
    ##     predict(model, newx = test_metadata, type = "response", s = "lambda.min")
    ## })

    ## ## Combine predicted probabilities into a data frame
    ## predicted_probabilities_df <- do.call(cbind, predicted_probabilities)
    ## colnames(predicted_probabilities_df) <- names(test_metadata)

    ## ## Print predicted probabilities
    ## print(predicted_probabilities_df)
    return(models)
    
    
}
