

library(tidyverse)
library(igraph)
library(devtools)
library(parallel)
library(optparse)
load_all()
library(reticulate)
use_virtualenv("/mnt/grackle_env")
library(tensorflow)



                                        # Define the list of options
option_list <- list(
    make_option(c("-n", "--num_nodes"), type = "numeric",  help = "Number of nodes int he graph", metavar = "NUMBER"),
    make_option(c("-i", "--iterations"), type = "numeric", help = "Iterations to run", metavar = "NUMBER"),
    make_option(c("-g", "--num_per_group"), type = "numeric",  help = "Number of samples per group", metavar = "NUMBER"),
    make_option(c("-k", "--num_groups"), type = "numeric",  help = "Number of groups", metavar = "NUMBER"),
    make_option(c("-c", "--num_cores"), type = "numeric",  help = "Number of cores", metavar = "NUMBER")
    
    
)


## Parse the command line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)



breast_g <- get(load("../GRACKLE_data/data/breast_igraph_prob_1_cor_0_05.RData"))

num_nodes = opt$num_nodes
num_groups = opt$num_groups
num_per_group = opt$num_per_group

### noise test
noise_sequence <- c(seq(.3,.7, .2))
t_sequence <- c(0,seq(.1,.3,.1))

num_nodes = 400
num_groups = 5
num_per_group = 20


rand_graph_file = paste0("../GRACKLE_results/small_simulation_A/random_graph_", num_nodes, "nodes.RData")


if(file.exists(rand_graph_file)) {
    print("loading random graph")
    load(rand_graph_file)
} else {
    print('generating random graph')
                                       # generate random breast_g# generate random nework 
    g <- old_randomNetwork(prior_graph = breast_g, num_nodes = num_nodes)
    save(g,file = rand_graph_file)
}



for(x in t_sequence) {

    cat("\n")
    print(paste("t_sequence", x))

    if(x == 0) {
        new_g = g
    } else {
                new_g <- adjustTransitivity(g,target_transitivity = x)
    }
    
    network_t = transitivity(as_undirected(new_g), weights = NA)
    print(paste("network transivity: ", network_t))
    
    ## generate list of permuted networks with distinctly permuted modules
    set.seed(42)
    clusters <- cluster_louvain(as.undirected(new_g),weights = NA)
    membership <- clusters$membership
    membership_counts <- table(as.factor(membership))
    print(membership_counts)
    large_clusters <- names(membership_counts[order(membership_counts, decreasing = T)][1:5])

# Assign colors to each community
V(g)$color <- membership(clusters)

clusters_to_highlight <- lapply(large_clusters, function(x) clusters[[x]])
# Plot the graph with module overlay
pdf("./results/networks/network_modules.pdf")
plot(
  clusters, 
  g, 
  mark.groups = clusters_to_highlight,
  vertex.size = 3,  # Adjust node size as needed
  vertex.label = NA,  # Hide vertex labels for clarity
  edge.width = 0.3,  # Adjust edge width as needed
  edge.arrow.size = .15,
  layout = layout_with_lgl
)
dev.off()

## plot a single cluster
cluster_id <- large_clusters[1]

vertices_in_cluster <- which(membership(clusters) == cluster_id)

# Create a subgraph with only the vertices in the selected cluster
subgraph <- induced_subgraph(g, vids = vertices_in_cluster)


# Plot the subgraph
pdf("./results/networks/single_module.pdf")
plot(
  subgraph, 
  vertex.size = 10,  # Adjust node size as needed
  vertex.label = NA,  # Hide vertex labels for clarity
  edge.width = 0.5,
  edge.arrow.size = .5,# Adjust edge width as needed,
  layout = layout_with_lgl
)
dev.off()


 
## simulate expression using the base network
base_exp <- parallelSimulateExpression(g,max_expression = 300,num_samples = 100,iterations = 3)

## compare expression value distribution to actual data
load("./data/Breast/correlation_filtered_breast_expression.RData")
quantile(as.matrix(correlations_filtered_breast_expression))
quantile(as.matrix(base_exp))

## simulate metadata
metadata <- simulateMetadata(group_size = rep(20,5), group_labels = paste0("subgroup", 1:5), row_names = rownames(base_exp))

pdf("./results/simulations/plots/full_network_heatmap.pdf")
plotSimulatedExpression(base_exp,metadata)
dev.off()
#permuted_networks <- lapply(large_clusters, function(x) modulePermutations(g,membership,x, increase_constant = 2 ))

# simulate gene expression for the permuted networks
sim_exp <- lapply(large_clusters[1:5], function(x) parallelSimulateExpression(g, iterations = 3, max_expression = 100,num_samples = 20,select_nodes = which(membership == x),perturbation_constant = 4,noise_constant = .7))

## combine expression data
sim_exp <- lapply(sim_exp, as.data.frame)

combined_exp <- do.call(rbind, sim_exp)

rownames(combined_exp) <- paste0("sample", 1:nrow(combined_exp))


pdf("./results/simulations/plots/4_permuted_network_heatmap.pdf")
Heatmap(as.matrix(combined_exp), right_annotation = ha, show_row_names = FALSE, show_column_names = FALSE, show_heatmap_legend = FALSE )
dev.off()

  ## split the data into train and test
  g_adjacency <- as_adjacency_matrix(g)  
  
  set.seed(42)
  
  dat <- split_data(combined_exp, metadata , training_size = .7)
  
  pat_sim <- as.matrix(dat$train_metadata) %*% t(dat$train_metadata)
    
  
  min_vals <- apply(dat$train_expression,2, min)
  max_vals <- apply(dat$train_expression,2, max)

    ## min-max scale the input matrix
  dat$train_expression <-as.matrix( min_max_scale(dat$train_expression,min_vals,max_vals))
  dat$test_expression <- as.matrix(min_max_scale(dat$test_expression,min_vals,max_vals))
  
  
  i_seq <-seq(0,1,.5) 
  j_seq <-seq(0,1,.5)
  
  
  grid_search <- as.data.frame(expand.grid(i_seq,j_seq))
  names(grid_search) <- c("lambda_1", "lambda_2")
  grid_search$score <- 0
  total <- ncol(grid_search) * nrow(grid_search)  
  
  
    
for(i in 1:nrow(grid_search)){
      
      print(i/nrow(grid_search))
      ## run GRACKLE NMF

      grackle <- GRACKLE(
        Y = dat$train_expression,
        net_similarity = as.matrix(g_adjacency),
        patient_similarity = pat_sim,
        diff_threshold = 1e-6,
        lambda_1 = grid_search$lambda_1[i],
        lambda_2 =  grid_search$lambda_2[i],
        k = 5, 
        verbose = F,
        beta = .0)
      
      
      ## correspondence between selected W LV's and top loading gene modules
      grid_search$score[i] <- evaluationWrapper(test_expression = dat$test_expression,
                                                test_metadata = dat$test_metadata,
                                                H_train = grackle$H,
                                                k = 5,
                                                clusters = clusters,
                                                aligned_clusters = large_clusters)
}
    
  p1 <- ggplot(grid_search, aes(x = lambda_1, y= lambda_2, fill = score)) +
    geom_tile() +
    scale_fill_continuous(limits = c(0,1)) + 
  theme_classic()    

  ggsave("./results/simulations/scores/heatmap_grid_search_4.pdf",plot = p1)
    
    
    
    nmf <- runNMF(dat$train_expression,5, "lee", seed = 42)
    evaluationWrapper(test_expression = dat$test_expression,
                      test_metadata = dat$test_metadata,
                      H_train = nmf$H,
                      k = 5,
                      clusters = clusters,
                      aligned_clusters = large_clusters)
    
    grnmf <- runGRNMF(expression_matrix = dat$train_expression, k = 5, seed = 42, max_iter = 200, alpha = .1)
    evaluationWrapper(test_expression = dat$test_expression,
                      test_metadata = dat$test_metadata,
                      H_train = nmf$H,
                      k = 5,
                      clusters = clusters,
                      aligned_clusters = large_clusters)
    
    
    