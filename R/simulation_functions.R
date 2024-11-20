## Simulation functions

#' permuteNetwork
#' 
#'  permuteNetwork is a function for permuting the network while maintaining the the degree distribution
#'
#' @importFrom igraph rewire
#' @param g igraph graph object
#' @param edge_quantiles Numeric vector of length 5 with reported quantile values for graph edges
#' @return permuted graph with the original degree distribution.
#' @export
permuteNetwork <- function(g, edge_quantiles){
  g <- rewire(g, with = keeping_degseq(niter = 100))
  #g <- rewire(g, each_edge(p = 0.1, loops = FALSE))
  sm <- randomQuantileVector(ecount(g), edge_quantiles) 
  E(g)$op <- sm
  
  return(g)
}

#' randomNetwork
#' 
#' @description randomNetwork applies igraph functions for creating a random network
#' @importFrom igraph barabasi.game
#' @param n Numeric value indicating the size of the graph. (Default is 100)
#' @return random graph with the defined number of edges
#' @export
#' 
randomNetwork <- function(n = 100, plot = FALSE) {
  g <- barabasi.game(n = n, m = 3, directed = TRUE)
  V(g)$name = as.character(1:n)
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
#' @param g igraph graph object
#' @param iterations numeric indicating how many iterations to perform. (Default = 10)
#' @param max_expression Numeric value for the maximum_gene expression. (Default is 2000)
#' @param num_samples Numeric values for the number of samples to simulate. (Default is 5)
#' @export
simulateExpression <- function(g, edge_quantiles, iterations = 10, max_expression = 2000, num_samples = 5) {
 
      
    # Specifying global reaction parameters.
    rp<-new("rsgns.param",time=0,stop_time=1000,readout_interval=1000)
  
    # Specifying the reaction rate constant vector for following reactions: (1) Translation rate, (2) RNA degradation rate, (3) Protein degradation rate, (4) Protein binding rate, (5) unbinding rate, (6) transcription rate.
    rc <- c(0.002, 0.005, 0.005, 0.005, 0.01, 0.02)
    
    # count nodes
    n_nodes <- vcount(g)
    
   
      # Assigning initial values to the RNAs and protein products to each node randomly based on breast gene expression counts.
     
      V(g)$Ppop <- (sample(max_expression,n_nodes, rep=TRUE))
      V(g)$Rpop <- (sample(max_expression, n_nodes, rep=TRUE))
      
    #Declaring input data object
    rsg <- new("rsgns.data",network=g, rconst=rc)
  res <- lapply(1:iterations, function(x) {
    #Call the R function for SGN simulator
    capture.output(xx <- rsgns.rn(rsg, rp, timeseries=FALSE, sample=num_samples))
    
    return(xx$expression)
  })
    
    sum_matrix <- Reduce("+", res)

    average_matrix <- sum_matrix / iterations

    return(average_matrix)
    
}


#' networkSimilarity
#' 
#' @description networkSimilarity is a function that calculates the similarity between permuted 
#' @param g igraph graph object
#' @param edge_quantiles Numeric vector of length 5 with reported quantile values for graph edges
#' @param iterations numeric indicating how many iterations to perform. (Default = 10)
#' @param max_expression Numeric value for the maximum_gene expression. (Default is 2000)
#' @param num_samples Numeric values for the number of samples to simulate. (Default is 5)
#' @export