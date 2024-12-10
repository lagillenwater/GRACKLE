## set seed
set.seed(42)
setwd("~/OneDrive - The University of Colorado Denver/Projects/GRACKLE/")
library(tidyverse)
library(igraph)
library(devtools)
library(ggplot2)
library(ComplexHeatmap)
library(parallel)
load_all()

load("./data/Breast/directed_breast_igraph.RData")


# generate random breast_g# generate random nework 
g <- randomNetwork(prior_graph = breast_g)


## permute modules
# generate list of permuted networks with distinctly permuted modules
clusters <- cluster_walktrap(g,weights = NA)
membership <- clusters$membership
membership_counts <- table(as.factor(membership))
large_clusters <- names(membership_counts[order(membership_counts, decreasing = T)][1:5])
## simulate metadata
metadata <- simulateMetadata(group_size = rep(20,5), group_labels = paste0("subgroup", 1:5), row_names = rownames(base_exp))

noise_sequence <- c(seq(0,.8,.1))

load("./results/simulations/data/noisy_simulations.RData")

## split the data into train and test
g_adjacency <- as_adjacency_matrix(g)  

noise_simulation_results <- lapply(noisy_expression, function(x) {
  
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
}, mc.cores = detectCores() - 4)
      
      
      
      return(results)
      
      
})



save(noise_simulation_results, file = "./results/simulations/data/NMF_model_results.RData")


