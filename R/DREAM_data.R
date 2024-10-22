## These are functions for  aggregating the DREAM4 simulated data sets for testing GRACKLE.

AggregateExpression <- function(training_directory, perturbation = NULL) {
    ## find the files in the directory
    dirs <- list.dirs(training_directory, full.names = TRUE, recursive = FALSE)

    for(dir in dirs) {
        files <- list.files(dir)
        dir_exp <- list()
        
        if(is.null(perturbation)){
            files <- files[grepl("knockdown|knockout|wildtype|multifactorial", files)]
        } else {
            files <- files[grepl(perturbation, files)]
        }

        for(file in files) {
            file_name <- paste0(dir,"/",file)
            single_exp <- read.delim(file_name)
            print(head(single_exp))
        }
    }
    
}

