## Functions  for runnning the training and testing experiments

#' split_data
#' 
#'  split_data is a function for splitting the gene expression and patient metadata matrices into training and testing data.
#'
#' @param expression Data.frame of simulated gene expression data
#' @param metadata Data.frame of patient metadata
#' @return train_expression Data.frame of training expression data
#' @return test_expression Data.frame of testing expression data
#' @return train_metadata Data.frame of training metadata
#' @return test_metadata Data.frame of testing metadata
#' @export
split_data <- function(expression,metadata,  training_size = .7) {

    set.seed(NULL)
    train_index <- sample(seq_len(nrow(metadata)), size = training_size * nrow(metadata))

    train_expression <- expression[train_index, ]
    test_expression <- expression[-train_index, ]


    train_metadata <- metadata[train_index,]
    test_metadata <- metadata[-train_index,]

    return(list(train_expression = train_expression, test_expression = test_expression, train_metadata = train_metadata, test_metadata = test_metadata))
    
}
