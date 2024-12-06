  ### Simulation study A
  
## set seed
set.seed(42)
setwd("~/OneDrive - The University of Colorado Denver/Projects/GRACKLE/")
library(tidyverse)
library(igraph)
library(ComplexHeatmap)
library(parallel)

# # Assign random edge weights based on the distribution of in and out degree over regulatory relationships in GRAND Breast network
breast_expression <- read.csv("./data/Breast_expression.csv")
breast_network <- read.csv("./data/Breast_network.csv")

## pivot longer
breast_network_long <- breast_network %>%
  pivot_longer(cols = -X)

# ## Remove interactions below a threshold of 1 This threshold is based on the evidence of interaction, regardless of directionality.
breast_network_long <- breast_network_long %>%
 filter(value >= 1)

# change ensembl ID's to HGNC ID's
network_symbol_map <- ensemblToHGNC(breast_network_long$name)
expression_symbol_map <- ensemblToHGNC(breast_expression$X)

breast_network_long <- breast_network_long %>%
    left_join(network_symbol_map, by = c("name" = "ensembl_gene_id"), relationship = "many-to-many")

breast_network_long <- breast_network_long %>%
    filter(!is.na(symbol))

breast_network_long <- breast_network_long %>%
    dplyr::select(X,symbol,value) %>%
    rename(from = X,to = symbol,probability = value)

breast_expression <- breast_expression %>%
    left_join(expression_symbol_map, by = c("X" = "ensembl_gene_id"), relationship = "many-to-many")

breast_expression <- breast_expression %>%
    filter(!is.na(symbol)) %>%
    select(-X)

## filter expression and network
unique_genes <- unique(union(breast_network_long$from,breast_network_long$to))
filtered_breast_expression <- breast_expression %>%
  filter(symbol %in% unique_genes) %>%
  column_to_rownames('symbol') %>%
  t() %>%
  as.data.frame()

filtered_breast_network_long <- breast_network_long %>%
  filter(from %in% names(filtered_breast_expression) & to %in% names(filtered_breast_expression))


#Calculate the correlation coefficients between pairs of genes to determine
# directed_breast_network <- makeDirected(filtered_breast_expression,filtered_breast_network_long)
# save(directed_breast_network, file = "directed_breast_network.RData")
load("directed_breast_network.RData")

directed_breast_network <- directed_breast_network %>% filter(abs(correlation) > .4)

## remove self loops
## directed_breast_network <- directed_breast_network %>% filter(from != to)
 
## find percent positive and negative edges
positive_edges_percent <- sum(directed_breast_network$weight > 0)/nrow(directed_breast_network)

# ## Find network properties
# breast_g <- graph_from_data_frame(directed_breast_network)
# save(breast_g, file = "directed_breast_igraph.RData")
load("directed_breast_igraph.RData")
## in and out degree of breast network
in_degrees <- degree(breast_g,mode = "in")
out_degrees <- degree(breast_g,mode = "out")
degrees <- degree(breast_g)
  
  
# generate random breast_g# generate random nework 
g <- degree.sequence.game(out_degrees, in_degrees, method = "simple")
  
sm <- sample(E(breast_g)$weight, ecount(g), rep = FALSE)
E(g)$weight <- sm
g
is_connected(g)
components <- decompose(g)
  
g <- components[[which.max(sapply(components, vcount))]]
is_connected(g)
## permute modules
# generate list of permuted networks with distinctly permuted modules
  
clusters <- cluster_walktrap(g)
membership <- cut_at(clusters,no = 100)
membership_counts <- table(as.factor(membership))
large_clusters <- names(membership_counts[order(membership_counts, decreasing = T)][1:5])
  
  metadata <- as.data.frame(c(rep(1, 20), rep(2, 20), rep(3, 20), rep(4, 20), rep(5,20)))
  names(metadata) <- "label"
  rownames(metadata) <- rownames(combined_exp)
  
  metadata <- metadata %>%
    rownames_to_column("sample") %>%  
    mutate(values = 1) %>%
    pivot_wider(names_from = label, names_prefix = "group", values_from = values, values_fill = 0) %>%
    column_to_rownames("sample")
  
  save(metadata, file = "./results/simulations/data/metadata/simulated_metadata_100_360.RData")
  
  ## explore the increase percentage as it relates to the data separation
  num_cores <- detectCores() - 1
  
  mclapply(seq(.5,4,.5), function(i) {
    
    permuted_networks <- lapply(large_clusters, function(x) modulePermutations(g,membership,x, increase_percentage =i ))
  
    # simulate gene expression for the permuted networks
    sim_exp <- lapply(permuted_networks, function(x) simulateExpression(x, iterations = 3, max_expression = 3000, num_samples = 20))
  
    ## combine expression data
    sim_exp <- lapply(sim_exp, as.data.frame)
    sim_exp <- lapply(sim_exp, function(x) {
      x %>% rownames_to_column("gene")
    })
    combined_exp <- sim_exp %>% reduce(full_join, by = "gene")
  
    combined_exp <- combined_exp %>%
      column_to_rownames("gene")
    
  
 
      names(combined_exp) <- paste0("sample", 1:ncol(combined_exp))
  
    combined_exp <- as.data.frame(t(combined_exp))
  
    save(combined_exp, file = paste0("./results/simulations/data/expression/", i,"_simulated_expression_100_360.RData"))
  },mc.cores = num_cores )
  
  
  
  
  ## TODO check the quantiles of the simulated expression data against observed expression data
  # breast_exp_quantiles <- quantile(as.matrix(breast_expression))
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
