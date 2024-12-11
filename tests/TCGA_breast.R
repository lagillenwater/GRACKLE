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

load("./data/Breast/TCGA/Breast_metadata.RData")
breast_subtype_metadata <- metadata %>%
  as.data.frame() %>%
  dplyr::select(paper_BRCA_Subtype_PAM50) %>%
  drop_na(paper_BRCA_Subtype_PAM50) %>%
  rownames_to_column("ID") %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = paper_BRCA_Subtype_PAM50, values_from = value) %>%
  column_to_rownames("ID")

breast_subtype_metadata[is.na(breast_subtype_metadata)] <- 0


load("./data/Breast/TCGA/directed_breast_g_with_PAM50.RData")

## split the data into train and test
breast_g_adjacency <- as_adjacency_matrix(directed_breast_g_with_PAM50)

load("./data/Breast/TCGA/Breast_filtered_gene_expression_with_PAM50.RData")

expression <- as.data.frame(t(expression_data))


## results <- mclapply(1:10, function(y) {
    
dat <- split_data(expression, breast_subtype_metadata , training_size = .7)
    
   pat_sim <- as.matrix(dat$train_metadata) %*% t(dat$train_metadata)
    
    min_vals <- apply(dat$train_expression,2, min)
    max_vals <- apply(dat$train_expression,2, max)
    
    ## min-max scale the input matrix
    dat$train_expression <-as.matrix( min_max_scale(dat$train_expression,min_vals,max_vals))
    dat$test_expression <- as.matrix(min_max_scale(dat$test_expression,min_vals,max_vals))
    
    
    ## i_seq <-seq(0,1,.1)
    ## j_seq <-seq(0,1,.1)
    
    ## grid_search <- as.data.frame(expand.grid(i_seq,j_seq))
    ## names(grid_search) <- c("lambda_1", "lambda_2")
    ## grid_search$score <- 0
#    total <- ncol(grid_search) * nrow(grid_search)
    

    
#    for(i in 1:nrow(grid_search)){
## scores <- mclapply(  1:nrow(grid_search), function(i) {
      #   print(i/nrow(grid_search))
      
      ## run GRACKLE NMF
      
      grackle <- GRACKLE(
        Y = dat$train_expression,
        net_similarity = as.matrix(breast_g_adjacency),
        patient_similarity = pat_sim,
        diff_threshold = 1e-6,
        lambda_1 =0, # grid_search$lambda_1[i],
        lambda_2 =0, #  grid_search$lambda_2[i],
        k = 5,
        verbose = T,
        beta = 0)
      
      
      ## correspondence between selected W LV's and top loading gene modules
      score <- evaluationWrapper(test_expression = dat$test_expression,
                                                test_metadata = dat$test_metadata,
                                                H_train = grackle$H,
                                                k = 5,
                                                clusters = clusters,
                                 aligned_clusters = large_clusters)
##        return(score)
#    }
    ## }, mc.cores = 10)   

    
##    grid_search$score <- unlist(scores)
  print(score)  
    
##     nmf <- runNMF(dat$train_expression,5, "lee", seed = 42)
##     nmf_res <- evaluationWrapper(test_expression = dat$test_expression,
##                                  test_metadata = dat$test_metadata,
##                                  H_train = nmf$H,
##                                  k = 5,
##                                  clusters = clusters,
##                                  aligned_clusters = large_clusters)
    
##     grnmf <- runGRNMF(expression_matrix = dat$train_expression, k = 5, seed = 42, max_iter = 200, alpha = .1)
##     grnmf_res <- evaluationWrapper(test_expression = dat$test_expression,
##                                    test_metadata = dat$test_metadata,
##                                    H_train = grnmf$H,
##                                    k = 5,
##                                    clusters = clusters,
##                                    aligned_clusters = large_clusters)
    
    
##     return(list(grid_search = grid_search, nmf_res = nmf_res, grnmf_res = grnmf_res))
##   }, mc.cores = 3)
  
  
  

## save(results, file = "./results/simulations/data/TCGA_results.RData")


