## These are functions for  aggregating the DREAM4 simulated data sets for testing GRACKLE.

#' aggregateExpression
#' 
#'  This function reads in the simulated expression profiles from the DREAM4 challenge. Depending on the perturbation parameter, it will read in a single data type or all available. 
#'
#' @import dplyr
#' @param training_directory Path to the input directory of training data
#' @param perturbation The in silico perturbation condition. Default is NULL, in which case all perturbation conditions are combined. Other options include "wildtype", "knockdown", "knockouts", and "multifactorial"
#' @return expression Aggregate expression dataset
#' @return metadata Aggregate metadata labels
#' @export
aggregateExpression <- function(training_directory, perturbation = NULL) {
    ## find the files in the directory
    dirs <- list.dirs(training_directory, full.names = TRUE, recursive = FALSE)
    expression <- data.frame()
    metadata <- character()
    
    for(dir in dirs) {
        files <- list.files(dir)
        
        
        if(is.null(perturbation)){
            files <- files[grepl("knockdown|knockout|wildtype|multifactorial", files)]
        } else {
            files <- files[grepl(perturbation, files)]
        }

        for(file in files) {
            file_name <- paste0(dir,"/",file)
            subtype <- sub(".*_size10_(\\d+)_.*", "\\1", file)
            single_exp <- read.delim(file_name)
            names(single_exp) <- paste0(names(single_exp),".",subtype)
            expression <- bind_rows(expression, single_exp)
            

            single_metadata <- rep(subtype,nrow(single_exp))
            metadata <- c(metadata,single_metadata)
            
        }
    }

    metadata <- as.data.frame(metadata)
    names(metadata) <- "subgroup"
    data_rownames <- paste0("S",1:nrow(expression))

    metadata <- model.matrix(~subgroup-1, metadata)
    
    rownames(expression) <- data_rownames
    rownames(metadata) <- data_rownames


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


#' aggregateNetworks
#' 
#'  This function reads in the gene regulatory networks used to create the simulated expression data. . 
#'
#' @import igraph
#' @import dplyr
#' @param gold_standards_directory Directory where the gold standard networks exist. Gold standard networks are edgelists with the format (G1	G2	1)
#' @return network Aggregated network
#' @export
aggregateNetworks <- function(gold_standards_directory) {
    
    files <- list.files(gold_standards_directory)
    files <- files[!(grepl("Bonus",files))]

    network <- data.frame()
    for(file in files) {
        file_name <- paste0(gold_standards_directory,"/",file)
        subtype <- sub(".*_size10_(\\d+)_.*", "\\1", file)
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

