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
   # set.seed(42)
    in_degrees <- degree(prior_graph,mode = "in")
    out_degrees <- degree(prior_graph,mode = "out")
    g_diff <- 100
    scale_factor <- 1
    g <- degree.sequence.game(out_degrees, in_degrees, method = "simple")
    sm <- sample(E(prior_graph)$weight, ecount(g), rep = FALSE)
    E(g)$weight <- sm
        
    while(abs(g_diff) > 1) {
        
        ## random_edges <- sample(E(g), num_nodes * scale_factor)
        ## g_subgraph <- subgraph.edges(g, eids = random_edges)

        random_nodes<- sample(V(g), num_nodes * scale_factor)
        g_subgraph <- induced_subgraph(g, vids = random_nodes)
        
        print(length(V(g_subgraph)))
        
        if(connected){  
          components <- decompose(g_subgraph)
                  g_subgraph <- components[[which.max(sapply(components, vcount))]]
        }

        g_length <- length(V(g_subgraph))
        

        g_diff <- g_length - num_nodes
        ## print(g_diff)
        
        if(g_diff < 0) {
            scale_factor <- scale_factor + .1
        } else {
            scale_factor <- scale_factor - .1
        }

        ## print(scale_factor)
        
    }
  return(g_subgraph)
}


#' old_randomNetwork
#' #'
#' #' @description randomNetwork applies igraph functions for creating a random network
#' #' @importFrom igraph degree.sequence.game
#' #' @param prior_graph Igraph object of the directed, weighted prior network to model  the simulated network on.
#' #' @param connected Boolean to determine if the random graph should be directed or not. (Default is TRUE)
#' #' @param num_nodes Numerical vector representing the size of the graph (Default is 400)
#' #' @return random graph with similar in and out degree
#' #' @export
#' #'
 old_randomNetwork <- function(prior_graph, connected = TRUE, num_nodes = 400){
    set.seed(42)
     in_degrees <- degree(prior_graph,mode = "in")
     out_degrees <- degree(prior_graph,mode = "out")
     g_length <- 0
     scale_factor <- 5

     while(g_length <  num_nodes) {
         g <- degree.sequence.game(out_degrees, in_degrees, method = "simple")
         sm <- sample(E(prior_graph)$weight, ecount(g), rep = FALSE)
         E(g)$weight <- sm

         random_edges <- sample(E(g), num_nodes * scale_factor)
         g <- subgraph.edges(g, eids = random_edges)

         if(connected){
           components <- decompose(g)

           g <- components[[which.max(sapply(components, vcount))]]
         }
         g_length <- length(V(g))
         scale_factor <- scale_factor + .1
    #     (paste0(scale_factor, " : ", g_length))
     }
   print(g)
   return(g)
 }

#' adjustTransitivity
#' 
#' @description adjustTransitivity iteratively rewires the graph to achieve the desired transitivity.
#' @importFrom igraph degree.sequence.game rewire
#' @param graph Igraph object 
#' @param target_transitivity target_transitivity measure for the simulated graph. (Default is .5)
#' @param max_iter the maximum number of iterations to try to achieve the target transitivity. 
#' @return random graph with similar in and out degree
#' @export
#' 
adjustTransitivity <- function(graph, target_transitivity = .5, max_iter = 1000) {
  current_transitivity <- transitivity(graph, type = "global")

  while(current_transitivity < target_transitivity) {
    
    # Randomly rewire edges to adjust transitivity
    graph <- increaseTransitivity(graph)
    
    # Ensure the graph remains connected
    # if (!is_connected(graph)) {
    #   components <- clusters(graph)
    #   for (i in 2:length(components$csize)) {
    #     graph <- add_edges(graph, edges = c(sample(which(components$membership == 1),1), sample(which(components$membership == i),1)))
    #   }
    # }
    current_transitivity <- transitivity(graph, type = "global")
  }
  
  return(graph)
}

#' increaseTransitivity
#' 
#' @description Function for increasing the transitivity of the graph based on 
#' @importFrom igraph neighbors add_edges
#' @param graph Igraph object 
#' @return random graph with similar in and out degree
#' @export
#' 
increaseTransitivity <- function(graph) {
  for (v in sample(V(graph),2)) {
    neighbors_v <- neighbors(graph, v)
    n_neighbors <- length(neighbors_v)
    if(n_neighbors > 2) {
      for (i in sample(1:(n_neighbors - 1),1)) {
        for (j in sample((i + 1):n_neighbors,1)) {
          if (!are_adjacent(graph, neighbors_v[i], neighbors_v[j])) {
            result <- try(graph <- add_edges(graph, c(neighbors_v[i], neighbors_v[j])))
            
            if (inherits(result, "try-error")) {
              # Skip the iteration if an error occurs
              next
            }
            
          
          }
        }
      }
    }
  }
  return(graph)
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

#' networkNoise
#' 
#' networkNoise is a function for increasing the noise of the graph 
#'
#' @importFrom igraph modularity cluster_louvain rewire
#' @param graph igraph graph object
#' @param goal_modularity The goal modularity to acheive in the permuted graph
#' @param membership Vertex membership of igraph object
#' @param large_clusters list of large_cluster ids
#' @return permuted graph
#' @export
networkNoise<- function(graph, goal_modularity, membership, large_clusters) {
  
  current_modularity <- modularity(cluster_louvain(as.undirected(graph),weights = NA))
  
  weights <- E(g)$weight 
  
  while(current_modularity > goal_modularity) {

      old_edge <- sample(E(graph),1)
      graph <- delete_edges(graph, old_edge)
    new_edge <- sample(V(graph), 2)
    graph <- add_edges(graph, new_edge)
    E(graph)[length(E(graph))]$weight  <- sample(weights,1)
    
    # for(i in large_clusters) {
    #   # print(i)
    #   community_nodes <- which(membership == i)
    #   community_edges <- E(graph)[.inc(V(graph)[community_nodes])]
    #   rand_edge <- sample(1:length(community_edges),1)
    #   rand_nodes <- c(ends(graph, community_edges[rand_edge]))
    #   community_edge <- get.edge.ids(graph, rand_nodes)
    #   graph <- delete_edges(graph, community_edge)
    #   
    #   graph <- add_edges(graph,c(rand_nodes[1], sample(which(membership != i),1)))
    #               
    # 
    # }
    current_modularity <- modularity(cluster_louvain(as.undirected(graph),weights = NA))
   #  print(current_modularity)
  }
  return(graph)
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
    time <- sample(400:600,1)
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
  res <- mclapply(1:num_samples, function(x) {
    simulateExpression(g,max_expression, iterations, rc, select_nodes , perturbation_constant, noise_constant)
  }, mc.cores = 2)
  
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


