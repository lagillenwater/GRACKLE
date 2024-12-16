## Simulation functions

#' randomNetwork
#' 
#' @description randomNetwork applies igraph functions for creating a random network
#' @importFrom igraph degree.sequence.game
#' @param prior_graph Igraph object of the directed, weighted prior network to model  the simulated network on. 
#' @param connected Boolean to determine if the random graph should be directed or not. (Default is TRUE)
#' @param num_nodes Numerical vector representing the size of the graph (Default is 400)
#' @return random graph with similar in and out degree
#' @export
#' 
randomNetwork <- function(prior_graph, connected = TRUE, num_nodes = 400){
    set.seed(42)
    in_degrees <- degree(prior_graph,mode = "in")
    out_degrees <- degree(prior_graph,mode = "out")
    g_length <- 0
    
    while(g_length < num_nodes) {
        g <- degree.sequence.game(out_degrees, in_degrees, method = "simple")
        sm <- sample(E(prior_graph)$weight, ecount(g), rep = FALSE)
        E(g)$weight <- sm
        
        if(connected){  
            components <- decompose(g)
            
            g <- components[[which.max(sapply(components, vcount))]]
        }
        g_length <- length(V(g))
    }
  print(g)
  return(g)
}


#' modulePermutations
#' 
#'  modulePermutations is a function for creating permuted graphs based on permuting specific modules within the graph
#'
#' @param g igraph graph object
#' @param membership module membership based on graph clustering.
#' @param module_ids Numeric vector for the membership ids to permute.
#' @param increase_constant Numeric value to scale the edge weights by. Default is 1.15
#' @return permuted graph with the original degree distribution.
#' @export
modulePermutations <- function(g,membership, module_ids, increase_constant = 1.15){
  
  module_vertices <- which(membership %in% module_ids)
  non_module_vertices <- which(!(membership %in% module_ids))

  # Get the edges connected to these vertices
  edges_of_interest <- E(g)[.from(module_vertices) | .to(module_vertices)]
  other_edges <- E(g)[.from(non_module_vertices) | .to(non_module_vertices)]
  
  # Extract the weights of these edges
  edge_weights <- E(g)$weight[edges_of_interest]
  other_edge_weights <- E(g)$weight[other_edges]
  
  # increase the edge weights by a defined factor
  scaled_weights <- edge_weights * increase_constant
  other_scaled_weights <- other_edge_weights * .8
  
  E(g)$weight[edges_of_interest] <- scaled_weights
  E(g)$weight[other_edges] <- other_scaled_weights
  
  return(g)
}



#' simualateExpression
#' 
#' @description simulateExpression is a function for simulating gene expression data from a GRN. It applies the functions outlined in the sgresR package vignette.
#' @import sgnesR 
#' @importFrom igraph V degree
#' @param g igraph graph object
#' @param max_expression Numeric value for the maximum_gene expression. (Default is 2000)
#' @param iterations Numeric value for how many iterations to perform and find average (Default is 3)
#' @param rc Numeric vector of length 5 representing recaction and decay rates.
#' @param select_nodes Numeric vector of specific nodes to perturb
#' @param perturbation_constant Numberic value to increase initial values of perturbed nodes by. (Default is 1.5)
#' @param noise_constant  Numeric value representing how much of the rest of the matrix to add noise to. 
#' @export
simulateExpression <- function(g, max_expression = 2000, iterations = 3, rc = NULL, select_nodes = NULL, perturbation_constant = 1.5, noise_constant = .3) {
  
  res <- lapply(1:iterations, function(x) {
    # Specifying global reaction parameters.
    time <- sample(600:800,1)
    rp<-new("rsgns.param",time=0,stop_time=time,readout_interval= time)
    
    # Specifying the reaction rate constant vector for following reactions: (1) Translation rate, (2) RNA degradation rate, (3) Protein degradation rate, (4) Protein binding rate, (5) unbinding rate, (6) transcription rate.
  
    if(is.null(rc)) {
      rc <- c(0.002, 0.005, 0.005, 0.005, 0.01, 0.02)
    } 
    # count nodes
    n_nodes <- vcount(g)
  
    # Assigning initial values to the RNAs and protein products to each node randomly based on breast gene expression counts.
    V(g)$Ppop <- sample(max_expression,n_nodes, rep=TRUE)
    V(g)$Rpop <- sample(max_expression, n_nodes, rep=TRUE)
    
    if(!is.null(select_nodes)){
      
      
      # select_degrees <- degree(g, v = select_nodes)
      # high_degree_nodes <- select_nodes[select_degrees > 2]
      # # <- sample(high_degree_nodes,1,replace = F)
      nodes <- select_nodes
      
      
      V(g)$Ppop[nodes] <- V(g)$Ppop[nodes] * perturbation_constant
      V(g)$Rpop[nodes] <- V(g)$Rpop[nodes] * perturbation_constant
    }
      
    if(noise_constant > 0){
      noise_nodes <- n_nodes* noise_constant
      other_nodes <- sample(V(g)[-select_nodes],noise_nodes) 
  
      V(g)$Ppop[other_nodes] <- V(g)$Ppop[other_nodes]*perturbation_constant
      V(g)$Rpop[other_nodes] <- V(g)$Rpop[other_nodes]*perturbation_constant
  
    }
  
    #Declaring input data object
    rsg <- new("rsgns.data",network=g, rconst=rc)
  
    #Call the R function for SGN simulator
    xx <- rsgns.rn(rsg, rp, timeseries=FALSE, sample=1)
  
    return(xx$expression)
  })
  
  sum_matrix <- Reduce("+", res)
  
  average_matrix <- sum_matrix / iterations
  
  return(average_matrix)

}

#' parallelSimulateExpression
#' 
#' @description parallelizes the simulate expression function to simulate the gene expression for multiple samples in the same subgroup
#' @importFrom parallel mclapply detectCores
#' @param g igraph graph object
#' @param max_expression Numeric value for the maximum_gene expression. (Default is 2000)
#' @param num_samples Numeric values for the number of samples to simulate. (Default is 5)
#' @param iterations numeric value to pass to simulation function. (Default is 3)
#' @param select_nodes Numeric vector of specific nodes to perturb
#' @param perturbation_constant Numberic value to increase initial values of perturbed nodes by. (Default is 1.5)
#' @param noise_constant  Numeric value representing how much of the rest of the matrix to add noise to. 
#' @return data frame of simulated gene expression data
#' @export
parallelSimulateExpression <- function(g,max_expression = 2000, num_samples= 5, iterations = 3, rc = NULL, select_nodes = NULL, perturbation_constant = 1.5,noise_constant = .3) {
  num_cores <- detectCores() - 1
  res <- mclapply(1:num_samples, function(x) {
    simulateExpression(g,max_expression, iterations, rc, select_nodes , perturbation_constant, noise_constant)
  })
  
  exp <- do.call(cbind,res)

  colnames(exp) <- paste0("sample", 1:ncol(exp))
  exp <- as.data.frame(t(exp))

  return(exp)
    
}



#' simulateMetadata
#' 
#' @description generates simulated metadata based on cluster assignments
#' @param group_size Numeric vector for the size of each group
#' @param group_labels Character vector of group labels
#' @param row_names Character vector of rownames to assign the metadata
#' @export
simulateMetadata <- function(group_size, group_labels, row_names) {
  
  
  metadata <- as.data.frame(rep(group_labels, times = group_size))
  names(metadata) <- "subgroup"
  rownames(metadata) <- row_names

  metadata <- metadata %>%
    rownames_to_column("sample") %>%  
    mutate(values = 1) %>%
    pivot_wider(names_from = subgroup,  values_from = values, values_fill = 0) %>%
    column_to_rownames("sample")

  return(metadata)
}

#' plotSimulatedExpression
#' 
#' @description plots heatmap of simulated expression
#' @importFrom RColorBrewer brewer.pal
#' @param expression simulated gene expression data. sample X genes
#' @param metadata simulated metadata. sample by column
#' @return
#' @export
plotSimulatedExpression <- function(expression,metadata){
  subgroups <- names(metadata)
  colors <- brewer.pal(length(subgroups), "Set1")
  group_colors <- setNames(colors, subgroups)
  ha <- rowAnnotation(subgroup = metadata %>% pivot_longer(cols = everything()) %>% filter(value == 1) %>% .$name, col = list(subgroup = group_colors), show_legend = FALSE)
  h1 <- Heatmap(as.matrix(expression), right_annotation = ha, show_row_names = FALSE, show_column_names = FALSE, show_heatmap_legend = FALSE )
  return(h1)

}


