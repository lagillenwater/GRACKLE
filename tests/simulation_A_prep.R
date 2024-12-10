## set seed
set.seed(42)
setwd("~/OneDrive - The University of Colorado Denver/Projects/GRACKLE/")
library(tidyverse)
library(igraph)
library(devtools)
library(ggplot2)
library(ComplexHeatmap)
load_all()

load("./data/Breast/directed_breast_igraph.RData")

## in and out degree of breast network

degrees <- degree(breast_g)

pdf("./results/networks/degree_scatter.pdf")
qplot(y  = seq_along(degrees), x= degrees, geom = "point", ylab = "" ) + theme_classic()
dev.off()
  
# generate random breast_g# generate random nework 
g <- randomNetwork(prior_graph = breast_g)

## plot network
pdf("./results/networks/network.pdf")
plot(
  g, 
  vertex.size = 3,  # Adjust node size as needed
  vertex.label = NA,  # Hide vertex labels for clarity
  edge.width = 0.3,  # Adjust edge width as needed
  edge.arrow.size = .15,
  layout = layout_with_lgl
)
dev.off()


## permute modules
# generate list of permuted networks with distinctly permuted modules
clusters <- cluster_walktrap(g,weights = NA)
membership <- clusters$membership
membership_counts <- table(as.factor(membership))
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

### noise test
noise_sequence <- c(seq(0,.8,.1))


set.seed(42)

noisy_expression <- lapply(noise_sequence, function(z){
  print(z)
  # simulate gene expression for the permuted networks
  sim_exp <- lapply(large_clusters, function(x) parallelSimulateExpression(g, iterations = 3, max_expression = 100,num_samples = 20,select_nodes = which(membership == x),perturbation_constant = 2.8,noise_constant = z))

  ## combine expression data
  sim_exp <- lapply(sim_exp, as.data.frame)
  combined_exp <- do.call(rbind, sim_exp)
  rownames(combined_exp) <- paste0("sample", 1:nrow(combined_exp))
  
  metadata <- simulateMetadata(group_size = rep(20,5), group_labels = paste0("subgroup", 1:5), row_names = paste0("sample", 1:nrow(combined_exp)))

  
   return(combined_exp)
})

names(noisy_expression) <- noise_sequence

lapply(noisy_expression,dim)

metadata <- simulateMetadata(group_size = rep(20,5), group_labels = paste0("subgroup", 1:5), row_names = paste0("sample", 1:nrow(noisy_expression[[1]])))

lapply(noise_sequence, function(x){
  pdf(paste0("./results/simulations/plots/permuted_network_heatmap_",x,"_noise.pdf"))
  print(plotSimulatedExpression(noisy_expression[[as.character(x)]],metadata))
  dev.off()
})

save(noisy_expression, file = "./results/simulations/data/noisy_simulations.RData")
