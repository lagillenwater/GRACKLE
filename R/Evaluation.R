## Evaluation functions
## Projecting data into the factorized latent space

#' project_W
#' 
#'  Function for predicting categorical labels in test data
#'
#' @import nnls
#' @param test_expression Data.frame of test expression data
#' @param H_train Matrix of trained H matrix
#' @param k Numeric value representing the number of groups
#' @return W_test Matrix of project W values for the test data
#' @export
project_W <- function(test_data,H_train,k) {
  
  k <- nrow(H_train)
  ## Initialize W_test
  W_test <- matrix(0, nrow=nrow(test_data), ncol=k)
  
  row_names <- rownames(test_data)
  
  ## convert test_data to a matrix
  #    test_data <- as.matrix(test_data)
  test_data <- apply(test_data,2,scale)
  
  ## Solve for W_test using NNLS
  for (i in 1:nrow(test_data)) {
    nnls_result <- nnls(t(H_train), test_data[i, ])
    W_test[i, ] <- coef(nnls_result)
  }
  
  rownames(W_test) <- row_names
  colnames(W_test) <- rownames(H_train)
  
  
  return(W_test)
}


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
  
    ## name LV's
    n_lvs <- ncol(W_test)
    colnames(W_test) <- paste0("LV", 1:n_lvs)
  
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
#' @param clusters Igraph object that is a list of cluster (aka, module) assignments
#' @param aligned_clusters The cluster labels aligned to the metadata group assignments
#' @return top A list of the module score. 
#' @export

geneLoadingsEvaluation <- function(H_train, clusters, aligned_clusters ) {
    # name LV's

    n_lvs <- nrow(H_train)
    rownames(H_train) <- paste0("LV", 1:n_lvs)
  
    ## calculate the module scores
    module_scores <- lapply(aligned_clusters, function(i) {
        rowSums(H_train[,clusters[[i]]])
    })
    
   
    module_scores <- lapply(module_scores, function(x) {
      names(x) <- paste0("LV", 1:n_lvs)
      return(x)
    })

    ## top modules
    top <- lapply(module_scores, function(x) {
        names(x)[x == max(x)]
    })
    names(top) <- paste0("subgroup",1:length(top))
    
    return(list(module_scores = module_scores, top = top))

}



#' sampleLoadingsEvaluation
#' 
#'  sampleLoadingsEvaluation is a function for identifying the top loading genes on each latent variable and aligning them with input data
#'
#' @param W_test Matrix of projected W from test data
#' @param module_size Size of consecutive modules for DREAM simulations.
#' @return top A list of the module score. 
#' @export

sampleLoadingsEvaluation <- function(W_test, test_metadata) {
  # name LV's
  ## name LV's
  n_lvs <- ncol(W_test)
  colnames(W_test) <- paste0("LV", 1:n_lvs)
  
  test_metadata <- as.data.frame(test_metadata)
  
  subgroups <- names(test_metadata)
  
    ## calculate the module scores
  subgroup_scores <- lapply(subgroups, function(i) {
    subgroup_samples <- rownames(test_metadata)[test_metadata[,i] == 1]
    num_samples = length(subgroup_samples)
    if(num_samples == 0) {
      res <- as.data.frame(t(W_test[1,]))
      res[1,] <-0
      return(res)
    } else if(num_samples ==1) {
      W_test[subgroup_samples,]
    } else {
      colSums(W_test[rownames(W_test) %in% subgroup_samples,])
    }
  })
  
  ## top coefficient
  top <- lapply(subgroup_scores, function(x) {
    names(x)[x == max(x)]
  })
  
  top <- lapply(top, function(x) {
    if(length(x) > 1) { NA} else {x}
  })
  
  names(top) <- names(test_metadata)
  
  return(list(subgroup_scores = subgroup_scores, top = top))
  
}


#' evaluationWrapper
#' 
#'  evaluationWrapper is function for finding the alignment of results
#'
#' @param test_expression Data.frame of test expression data
#' @param test_metadata Data.frame of the test metadata
#' @param H_train Matrix of trained H matrix
#' @param k Numeric value representing the number of groups
#' @param clusters igraph communities object with cluster membership
#' @param aligned_clusters clusters aligned with metadata assignments
#' @return 
#' @export

evaluationWrapper <- function(test_expression,  test_metadata,H_train,k,clusters,  aligned_clusters) {

  W_test <-  project_W(test_expression, H_train)
  
  ## evaluate sample loadings
  top_sample_LVs <- sampleLoadingsEvaluation(W_test, test_metadata)
  
  ## evaluate gene loadings
  top_loadings <- geneLoadingsEvaluation(H_train, clusters = clusters, aligned_clusters = aligned_clusters)
  #   
  ## correspondence between selected W LV's and top loading gene modules
  score <- mean(unlist(lapply(1:5, function(x) {
    if(!(is.na(top_sample_LVs$top[x][[1]])) | is.na(top_loadings$top[x][[1]])) {
      if(identical(top_sample_LVs$top[x][[1]], top_loadings$top[x][[1]])) { 1} else {0}
    } })))
  
  return(score)
  
  
  
}