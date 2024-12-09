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
in_degrees <- degree(breast_g,mode = "in")
out_degrees <- degree(breast_g,mode = "out")
degrees <- degree(breast_g)

pdf("./results/networks/degree_scatter.pdf")
qplot(y  = seq_along(degrees), x= degrees, geom = "point", ylab = "" ) + theme_classic()
dev.off()
  
# generate random breast_g# generate random nework 
g <- degree.sequence.game(out_degrees, in_degrees, method = "simple")
sm <- sample(E(breast_g)$weight, ecount(g), rep = FALSE)
E(g)$weight <- sm
g

is_connected(g)
components <- decompose(g)
  
g <- components[[which.max(sapply(components, vcount))]]
is_connected(g)
g

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
base_exp <- parallelSimulateExpression(g,max_expression = 300,num_samples = 50,iterations = 3)

## compare expression value distribution to actual data
load("./data/Breast/correlation_filtered_breast_expression.RData")
quantile(as.matrix(correlations_filtered_breast_expression))
quantile(as.matrix(base_exp))

## simulate metadata
metadata <- simulateMetadata(group_size = rep(10,5), group_labels = paste0("subgroup", 1:5), row_names = rownames(base_exp))

group_colors <- c("subgroup1" = "black",  "subgroup2" ="darkorange", "subgroup3" ="darkgreen","subgroup4" ="red","subgroup5" ="blue")
ha <- rowAnnotation(subgroup = metadata %>% pivot_longer(cols = everything()) %>% filter(value == 1) %>% .$name, col = list(subgroup = group_colors), show_legend = FALSE)
pdf("./results/simulations/plots/full_network_heatmap.pdf")
Heatmap(as.matrix(base_exp), right_annotation = ha, show_row_names = FALSE, show_column_names = FALSE, show_heatmap_legend = FALSE )
dev.off()

#permuted_networks <- lapply(large_clusters, function(x) modulePermutations(g,membership,x, increase_constant = 2 ))

# simulate gene expression for the permuted networks
sim_exp <- lapply(large_clusters, function(x) parallelSimulateExpression(g, iterations = 3, max_expression = 200,num_samples = 3,select_nodes = which(membership == x),perturbation_constant = 3 ))

## combine expression data
sim_exp <- lapply(sim_exp, as.data.frame)

combined_exp <- do.call(rbind, sim_exp)

rownames(combined_exp) <- paste0("sample", 1:nrow(combined_exp))

## simulate metadata
metadata <- simulateMetadata(group_size = rep(3,3), group_labels = paste0("subgroup", 1:3), row_names = rownames(combined_exp))

group_colors <- c("subgroup1" = "black",  "subgroup2" ="darkorange", "subgroup3" ="darkgreen","subgroup4" ="red","subgroup5" ="blue")
ha <- rowAnnotation(subgroup = metadata %>% pivot_longer(cols = everything()) %>% filter(value == 1) %>% .$name, col = list(subgroup = group_colors), show_legend = FALSE)
 Heatmap(as.matrix(combined_exp), right_annotation = ha, show_row_names = FALSE, show_column_names = FALSE, show_heatmap_legend = FALSE )

pdf("heatmap.pdf")
draw(p1)
dev.off()


## explore the increase percentage as it relates to the data separation
num_cores <- detectCores() - 1
  
mclapply(seq(.5,4,.5), function(i) {
    

  
   
  
    save(combined_exp, file = paste0("./results/simulations/data/expression/", i,"_simulated_expression_100_360.RData"))
  },mc.cores = num_cores )
  
  
  
  
  ## TODO check the quantiles of the simulated expression data against observed expression data
  
  # lapply(sim_exp, quantile)
  
  load("./results/simulations/data/expression/3_simulated_expression_100_360.RData")
 
  
  ## split the data into train and test
  dat <- split_data(combined_exp, metadata , training_size = .7)
  
  set.seed(42)
  
  
 
  
  pat_sim <- as.matrix(dat$train_metadata) %*% t(dat$train_metadata)
  
  g_adjacency <- as_adjacency_matrix(g)
  
  net_sim <- similarityCalc(as.matrix(as_adjacency_matrix(g)))
  # pat_sim[pat_sim < 1] <- 0
  min_vals <- apply(dat$train_expression,2, min)
  max_vals <- apply(dat$train_expression,2, max)
  ## min-max scale the input matrix
  Y <-as.matrix( min_max_scale(dat$train_expression,min_vals,max_vals))
  
  i_seq <-seq(0,1,.1) 
  j_seq <-seq(0,1,.1)
  
  
  grid_search <- as.data.frame(expand.grid(i_seq,j_seq))
  names(grid_search) <- c("lambda_1", "lambda_2")
  grid_search$score <- 0
  total <- ncol(grid_search) * nrow(grid_search)  
  
  
    
for(i in 1:nrow(grid_search)){
      
      print(grid_search[i,])
      ## run GRACKLE NMF
      g_res <- GRACKLE(
        Y = Y,
        net_similarity = as.matrix(g_adjacency),
        patient_similarity = pat_sim,
        diff_threshold = 1e-6,
        lambda_1 = grid_search$lambda_1[i],
        lambda_2 = grid_search$lambda_2[i],
        k = 5, 
        verbose = F,
        beta = .0)
      
      ## project the test data into W_test
      W_test <- project_W(min_max_scale(dat$test_expression,min_vals,max_vals) ,g_res$H)
      
      ## evaluate sample loadings
      top_sample_LVs <- sampleLoadingsEvaluation(W_test,dat$test_metadata)
      #print(top_sample_LVs$subgroup_scores)
        
      ## evaluate gene loadings
      top_loadings <- geneLoadingsEvaluation(g_res$H, clusters = clusters, aligned_clusters = large_clusters)
     # print(top_loadings$module_scores)
      
      ## correspondence between selected W LV's and top loading gene modules
      grid_search$score[i] <- mean(unlist(lapply(1:5, function(x) {
        if(!(is.na(top_sample_LVs$top[x][[1]])) | is.na(top_loadings$top[x][[1]])) {
          if(identical(top_sample_LVs$top[x][[1]], top_loadings$top[x][[1]])) { 1} else {0}
        } })))
}
    
  p1 <- ggplot(grid_search, aes(x = lambda_1, y= lambda_2, fill = score)) +
    geom_tile() +
    scale_fill_continuous(limits = c(0,1))
  p1
  ggsave("./results/simulations/scores/heatmap_grid_search_25_360.pdf",plot = p1)
    
    
    
   l_res <- nmf(as.matrix(Y),5,'lee', seed = 42)
     
   W_test <-  project_W(min_max_scale(dat$test_expression,min_vals,max_vals), coef(l_res))
     ## evaluate sample loadings
     top_sample_LVs <- sampleLoadingsEvaluation(W_test,dat$test_metadata)
  #   
     ## evaluate gene loadings
     top_loadings <- geneLoadingsEvaluation(g_res$H, clusters = clusters, aligned_clusters = large_clusters)
  #   
   ## correspondence between selected W LV's and top loading gene modules
     mean(unlist(lapply(1:5, function(x) {
       if(!(is.na(top_sample_LVs$top[x][[1]])) | is.na(top_loadings$top[x][[1]])) {
         if(identical(top_sample_LVs$top[x][[1]], top_loadings$top[x][[1]])) { 1} else {0}
       } })))

  #   
  # 
  # toplot <- data.frame(noise = seq(.1,.9,.2),  nmf = as.numeric(res[2,]))
  # library(tidyverse)
  # toplot <- toplot %>%
  #   pivot_longer(cols = -noise)
  # 
  # ggplot(toplot, aes(x = noise, y = value, color = name)) +
  #   ylim(c(0,1))+
  #   geom_line(linewidth = 2) +
  #   theme_classic()
  # 
