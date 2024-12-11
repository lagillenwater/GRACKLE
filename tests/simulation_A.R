## set seed
set.seed(42)
#setwd("~/OneDrive - The University of Colorado Denver/Projects/GRACKLE/")
library(tidyverse)
library(igraph)
library(devtools)
library(ggplot2)
library(ComplexHeatmap)
library(parallel)
load_all()

load("./data/Breast/directed_breast_igraph.RData")

degrees <- degree(breast_g)

pdf("./results/networks/degree_scatter.pdf")
qplot(y  = seq_along(degrees), x= degrees, geom = "point", ylab = "" ) + theme_classic()
dev.off()

# generate random breast_g# generate random nework 
g <- randomNetwork(prior_graph = breast_g, num_nodes = 400)

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


### noise test
noise_sequence <- c(seq(0,.8,.1))

noisy_expression <- mclapply(noise_sequence, function(z){
  print(z)
  # simulate gene expression for the permuted networks
  sim_exp <- mclapply(large_clusters, function(x) parallelSimulateExpression(g, iterations = 3, max_expression = 100,num_samples = 20,select_nodes = which(membership == x),perturbation_constant = 2.8,noise_constant = z), mc.cores = length(large_clusters))
  
  ## combine expression data
  sim_exp <- lapply(sim_exp, as.data.frame)
  combined_exp <- do.call(rbind, sim_exp)
  rownames(combined_exp) <- paste0("sample", 1:nrow(combined_exp))
  
  metadata <- simulateMetadata(group_size = rep(20,5), group_labels = paste0("subgroup", 1:5), row_names = paste0("sample", 1:nrow(combined_exp)))
  
  
  return(combined_exp)
}, mc.cores = length(noise_sequence))

names(noisy_expression) <- noise_sequence

metadata <- simulateMetadata(group_size = rep(20,5), group_labels = paste0("subgroup", 1:5), row_names = paste0("sample", 1:nrow(noisy_expression[[1]])))

lapply(noise_sequence, function(x){
  pdf(paste0("./results/simulations/plots/permuted_network_heatmap_",x,"_noise.pdf"))
  print(plotSimulatedExpression(noisy_expression[[as.character(x)]],metadata))
  dev.off()
})

save(noisy_expression, file = "./results/simulations/data/noisy_simulations.RData")



noise_sequence <- c(seq(0,.8,.1))

metadata <- simulateMetadata(group_size = rep(20,5), group_labels = paste0("subgroup", 1:5), row_names = rownames(noisy_expression[[1]]))

## split the data into train and test
g_adjacency <- as_adjacency_matrix(g)

noise_simulation_results <- mclapply(noisy_expression, function(x) {

  results <- mclapply(1:10, function(y) {

    dat <- split_data(x, metadata , training_size = .7)

    pat_sim <- as.matrix(dat$train_metadata) %*% t(dat$train_metadata)

    min_vals <- apply(dat$train_expression,2, min)
    max_vals <- apply(dat$train_expression,2, max)

    ## min-max scale the input matrix
    dat$train_expression <-as.matrix( min_max_scale(dat$train_expression,min_vals,max_vals))
    dat$test_expression <- as.matrix(min_max_scale(dat$test_expression,min_vals,max_vals))


    i_seq <-seq(0,1,.1)
    j_seq <-seq(0,1,.1)

    grid_search <- as.data.frame(expand.grid(i_seq,j_seq))
    names(grid_search) <- c("lambda_1", "lambda_2")
    grid_search$score <- 0
    total <- ncol(grid_search) * nrow(grid_search)

    print(y)

    for(i in 1:nrow(grid_search)){

       #   print(i/nrow(grid_search))

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





        nmf <- runNMF(dat$train_expression,5, "lee", seed = 42)
        nmf_res <- evaluationWrapper(test_expression = dat$test_expression,
                          test_metadata = dat$test_metadata,
                        H_train = nmf$H,
                        k = 5,
                        clusters = clusters,
                        aligned_clusters = large_clusters)

      grnmf <- runGRNMF(expression_matrix = dat$train_expression, k = 5, seed = 42, max_iter = 200, alpha = .1)
      grnmf_res <- evaluationWrapper(test_expression = dat$test_expression,
                        test_metadata = dat$test_metadata,
                        H_train = grnmf$H,
                        k = 5,
                        clusters = clusters,
                        aligned_clusters = large_clusters)


      return(list(grid_search = grid_search, nmf_res = nmf_res, grnmf_res = grnmf_res))
}, mc.cores = 10)



      return(results)


}, mc.cores = 3)



save(noise_simulation_results, file = "./results/simulations/data/NMF_model_results.RData")


