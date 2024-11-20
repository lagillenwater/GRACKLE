## These are functions for  aggregating the DREAM4 simulated data sets for testing GRACKLE.

#' aggregateExpression
#' 
#'  This function reads in the simulated expression profiles from the DREAM4 challenge. Depending on the perturbation parameter, it will read in a single data type or all available. 
#'
#' @importFrom dplyr bind_rows
#' @param training_directory Path to the input directory of training data
#' @param perturbation The in silico perturbation condition. Default is NULL, in which case all perturbation conditions are combined. Other options include "wildtype", "knockdown", "knockouts", and "multifactorial"
#' @return A list containing the following components:
#' \item{expression}{Aggregate expression dataset}
#' \item{metadata}{Aggregate metadata labels}
#' @export
aggregateExpression <- function(training_directory, perturbation = NULL) {
    ## find the files in the directory
    dirs <- list.dirs(training_directory, full.names = TRUE, recursive = FALSE)
    expression <- data.frame()
    metadata <- data.frame()
    
    for(dir in dirs) {
        files <- list.files(dir)
        
        
        if(is.null(perturbation)){
            files <- files[grepl("knockdown|knockout|wildtype|multifactorial", files)]
        } else {
            files <- files[grepl(perturbation, files)]
        }

        for(file in files) {
            file_name <- paste0(dir,"/",file)
            if(grepl("100", file_name)) {
              subtype <- sub(".*_size100_(\\d+)_.*", "\\1", file)
            } else if(grepl("10", file_name)) {
              subtype <- sub(".*_size10_(\\d+)_.*", "\\1", file)
            }
            single_exp <- read.delim(file_name)
            names(single_exp) <- paste0(subtype, "|",names(single_exp))
            row_labels <- paste0(subtype,"|", 1:nrow(single_exp))
            rownames(single_exp) <- row_labels
            expression <- bind_rows(expression, single_exp)
            

            single_metadata <- as.data.frame(rep(subtype,nrow(single_exp)))
            rownames(single_metadata) <- row_labels
            metadata <- bind_rows(metadata,single_metadata)
            
            
        }
    }

    #metadata <- as.data.frame(metadata)
    names(metadata) <- "subgroup"
   # data_rownames <- paste0("S",1:nrow(expression))

    metadata <- model.matrix(~subgroup-1, metadata)
    
  #  rownames(expression) <- data_rownames
   # rownames(metadata) <- data_rownames


    return(list(expression = expression, metadata = metadata))
    
    
}

#' imputeExpression
#' 
#'  This function imputes missing values in the gene expression matrix randomly. 
#'
#' @param expression_matrix Data.frame of simulated gene expression data. 
#' @return expression_matrix Imputed matrix
#' @export
imputeExpression <- function(expression_matrix) {
    exp_mean <- mean(unlist(expression_matrix), na.rm = T)
    exp_sd <- sd(unlist(expression_matrix), na.rm = T)
    expression_matrix[is.na(expression_matrix)] <- abs(rnorm(sum(is.na(expression_matrix)), exp_mean, exp_sd))
    return(expression_matrix)
}

#' noisyImputation
#' 
#'  noisyImputation is a function for imputing data with noise while maintaining the original correlation substructure
#'
#' @param expression_matrix Data.frame of simulated gene expression data. 
#' @param scaling_factor A numeric value to scale the imputed matrix by. The default is .1.
#' @return expression_matrix Imputed matrix
#' @export
noisyImputation <- function(expression_matrix, scaling_factor = .1) {
  n_col <- ncol(expression_matrix)
  n_row <- nrow(expression_matrix)
  noise <- matrix(rnorm(n_row * n_col ), ncol = n_col, nrow = n_row)
  scaled_noise <- noise * scaling_factor
  
  expression_matrix[is.na(expression_matrix)] <- 0
  noisy_gene_expression <- expression_matrix + scaled_noise
  
  preserved_noisy_gene_expression <- preserveCorrelation(expression_matrix,noisy_gene_expression)
  
  return(noisy_gene_expression)
  
  
  }

#' preserveCorrelation
#' 
#'  preserveCorrelation is a function that uses pca to ensure that noisy data maintains the original correlation structure
#'
#' @param original A matrix of the original gene expression values. NA values set to 0.  
#' @param noisy A matrix of the noisy gene expression values, including NA values. 
#' @return expression_matrix Imputed matrix
#' @export
# Function to preserve correlation structure using PCA
preserveCorrelation <- function(original, noisy) {
  pca <- prcomp(original, scale. = TRUE)
  scores <- pca$x
  loadings <- pca$rotation
  reconstructed <- scores %*% t(loadings)
  return(reconstructed)
}

#' aggregateNetworks
#' 
#'  This function reads in the gene regulatory networks used to create the simulated expression data. . 
#'
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph as_adjacency_matrix
#' @importFrom dplyr bind_rows
#' @param gold_standards_directory Directory where the gold standard networks exist. Gold standard networks are edgelists with the format (G1	G2	1)
#' @return network Aggregated network
#' @export
aggregateNetworks <- function(gold_standards_directory) {
    
    files <- list.files(gold_standards_directory)
    files <- files[!(grepl("Bonus",files))]

    network <- data.frame()
    for(file in files) {
        file_name <- paste0(gold_standards_directory,"/",file)
        if(grepl("100", file_name)) {
          subtype <- sub(".*_size100_(\\d+)_.*", "\\1", file)
        } else if(grepl("10", file_name)) {
          subtype <- sub(".*_size10_(\\d+)_.*", "\\1", file)
        }
        single_net <- read.delim(file_name, header = F)
        single_net$V1 <- paste0(single_net$V1,".",subtype)
        single_net$V2 <- paste0(single_net$V2,".",subtype)
        names(single_net) <- c("from", "to", "weight")
        single_graph <- graph_from_data_frame(single_net, directed = TRUE)

        single_adjacency <- as.data.frame(as_adjacency_matrix(single_graph,sparse =F, attr = "weight"))

        network <- bind_rows(network,single_adjacency)
        
    }

    
    network[is.na(network)] <- 0
    return(network)
}

#' permuteNetwork
#' 
#'  permuteNetwork is a function for permuting the network to test as a negative control
#'
#' @param network Matrix for the adjacency to permute
#' @return modified_matrix Matrix of the permuted network
#' @export

permuteNetwork <- function(network, replace_percentage = .1) {
  
  network <- as.matrix(network)
  permuted_network <-network[sample(nrow(network)), sample(ncol(network))]
  
  num_elements <- length(network)
  num_replace <- round(replace_percentage * num_elements)
  
  # Get indices to replace
  replace_indices <- sample(1:num_elements, num_replace)
  
  # Flatten matrices for easy indexing
  original_flat <- as.vector(network)
  permuted_flat <- as.vector(permuted_network)
  
  # Replace values
  original_flat[replace_indices] <- permuted_flat[replace_indices]
  
  # Reshape back to matrix
  modified_matrix <- matrix(original_flat, nrow = nrow(network), ncol = ncol(network))
  return(modified_matrix)
  
}