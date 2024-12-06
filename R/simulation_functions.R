## Simulation functions

#' randomNetwork
#' 
#' @description randomNetwork applies igraph functions for creating a random network
#' @importFrom igraph barabasi.game
#' @param n Numeric value indicating the size of the graph. (Default is 100)
#' @param edge_quantiles Numeric vector of length 5 with reported quantile values for graph edges
#' @return random graph with the defined number of edges
#' @export
#' 
randomNetwork <- function(n = 100,edge_quantiles) {
  g <- barabasi.game(n = n, m = 2, directed = TRUE)
  V(g)$name = as.character(1:n)
  
  sm <- randomQuantileVector(ecount(g), edge_quantiles) 
  E(g)$op <- sm
  
  return(g)
}


#' modulePermutations
#' 
#'  modulePermutations is a function for creating permuted graphs based on permuting specific modules within the graph
#'
#' @param g igraph graph object
#' @param membership module membership based on graph clustering.
#' @param module_ids Numeric vector for the membership ids to permute.
#' @param increase_percentage Numeric value to scale the edge weights by. Default is 1.15
#' @return permuted graph with the original degree distribution.
#' @export
modulePermutations <- function(g,membership, module_ids, increase_percentage = 1.15){
  
  module_vertices <- which(membership %in% module_ids)

  # Get the edges connected to these vertices
  edges_of_interest <- E(g)[.from(module_vertices) | .to(module_vertices)]
  
  # Extract the weights of these edges
  edge_weights <- E(g)$weight[edges_of_interest]
  
  # increase the edge weights by a defined factor
  scaled_weights <- edge_weights * increase_percentage
  
   E(g)$weight[edges_of_interest] <- scaled_weights
  
  return(g)
}



#' randomQantileVector
#' 
#' @description randomQunatileVector calculates a random vector based on an input quantile vector for quantiles (0,.25,.5,.75,1)
#' @param n Numeric value for the length of random vector
#' @param values Numeric vector of length 5 with reported quantile values
#' @return random quantile vector of length n
#' @export
#' 
randomQuantileVector <- function(n,  values) {
  quantiles <- c(0, 0.25, 0.5, 0.75, 1)
  # Generate uniform random numbers
  u <- runif(n)
  # Interpolate to get the corresponding values
  interpolated_values <- approx(quantiles, values, u)$y
  return(interpolated_values)
}

#' simualateExpression
#' 
#' @description simulateExpression is a function for simulating gene expression data from a GRN. It applies the functions outlined in the sgresR package vignette.
#' @import sgnesR 
#' @param g igraph graph object
#' @param iterations numeric indicating how many iterations to perform. (Default = 10)
#' @param max_expression Numeric value for the maximum_gene expression. (Default is 2000)
#' @param num_samples Numeric values for the number of samples to simulate. (Default is 5)
#' @export
simulateExpression <- function(g,  iterations = 10, max_expression = 2000, num_samples = 5) {
 
  res <- lapply(1:iterations, function(x) {
    
    # Specifying global reaction parameters.
    rp<-new("rsgns.param",time=0,stop_time=100,readout_interval= 100)
  
    # Specifying the reaction rate constant vector for following reactions: (1) Translation rate, (2) RNA degradation rate, (3) Protein degradation rate, (4) Protein binding rate, (5) unbinding rate, (6) transcription rate.
    
    rc <- c(0.002, 0.005, 0.005, 0.005, 0.01, 0.02)
   
    # count nodes
    n_nodes <- vcount(g)
   
      # Assigning initial values to the RNAs and protein products to each node randomly based on breast gene expression counts.
     
      V(g)$Ppop <- (sample(max_expression,n_nodes, rep=TRUE))
      V(g)$Rpop <- (sample(max_expression, n_nodes, rep=TRUE))
      
    #Declaring input data object
    rsg <- new("rsgns.data",network=g, rconst=rc)
   
    #Call the R function for SGN simulator
    capture.output(xx <- rsgns.rn(rsg, rp, timeseries=FALSE, sample=num_samples))
    
    return(xx$expression)
  })
    
    sum_matrix <- Reduce("+", res)

    average_matrix <- sum_matrix / iterations

    return(average_matrix)
    
}

#' makeDirected
#' 
#' @description makedDirected is a function for  calculating the pairwise correlations coefficients between genes and assigning direction to the network based on correlation. 
#' @importFrom parallel detectCores  mclapply
#' @importFrom ppcor pcor
#' @importFrom dplyr filter %>% select arrange mutate
#' @param expression Data Frame of gene expression data
#' @param network Data Frame of the network edgelist with column names 'from', 'to', 'weight'.
#' @return Dataframe of network edgelist with direction (1 = activation, -1 = suppression)
#' @export
makeDirected <- function(expression,network) {
    
  ## arrange network by targets
  network <- network %>%
    arrange(to)
  ## identify the target genes
  targets <- unique(network$to)
  
  ## find number of cores for parallelization
  num_cores <- detectCores() - 1
  correlations <- mclapply( targets, function(i) {
    ## filter network by targets
    tmp_network <- network %>%
      filter(to == i)
    tmp_expression <- expression %>%
      dplyr::select(c(tmp_network$from,i))
        
    ## calculate the partial correlations
    partial_correlations <- suppressWarnings(pcor(tmp_expression)$estimate)
    colnames(partial_correlations) <- names(tmp_expression)
    rownames(partial_correlations) <- names(tmp_expression)
    
    ## 
    res <- partial_correlations %>%
      as.data.frame  %>%
      dplyr::select(i)

    if(!(i %in% tmp_network$from)) {
      res <- res %>%
        filter(rownames(.) != i)
    }
    
    return(res)
    
  },mc.cores = num_cores)
  
  network$correlation <- c(unlist(correlations))
  network <- network %>%
    mutate(weight = ifelse(correlation > 0, probability,-probability))
  
  return(network)
}

