## These are functions for  aggregating the DREAM4 simulated data sets for testing GRACKLE.

#' AggregateExpression
#' 
#'  This function reads in the simulated expression profiles from the DREAM4 challenge. Depending on the perturbation parameter, it will read in a single data type or all available. 
#' 
#' @param training_directory Path to the input directory of training data
#' @param perturbation The in silico perturbation condition. Default is NULL, in which case all perturbation conditions are combined. Other options include "wildtype", "knockdown", "knockouts", and "multifactorial"
#' @return expression Aggregate expression dataset
#' @return metadata Aggregate metadata labels
#' @export
AggregateExpression <- function(training_directory, perturbation = NULL) {
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
            single_exp <- read.delim(file_name)
            expression <- rbind(expression, single_exp)

            single_metadata <- sub(".*_size10_(\\d+)_.*", "\\1", file)
            metadata <- rbind(metadata,single_metadata)
            
        }
    }

    names(metadata) <- "subgroup"
    data_rownames <- paste0("S",1:nrow(expression))

    rownames(expression) <- data_rownames
    rownames(metadata) <- data_rownames

    return(list(expression = expression, metadata = metadata))
    
    
}

