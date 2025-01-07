
#setwd("~/OneDrive - The University of Colorado Denver/Projects/GRACKLE/")
library(tidyverse)
library(igraph)
library(devtools)
library(ggplot2)
library(parallel)
library(aricode)
load_all()

load("../GRACKLE_data/data/Breast/TCGA/Breast_metadata.RData")

breast_subtype_metadata <- metadata %>%
  as.data.frame() %>%
  dplyr::select(paper_BRCA_Subtype_PAM50) %>%
  drop_na(paper_BRCA_Subtype_PAM50) %>%
  rownames_to_column("ID") %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = paper_BRCA_Subtype_PAM50, values_from = value) %>%
  column_to_rownames("ID")

breast_subtype_metadata[is.na(breast_subtype_metadata)] <- 0


load("../GRACKLE_data/data/Breast/TCGA/directed_breast_g_with_PAM50.RData")

## split the data into train and test
breast_g_adjacency <- as_adjacency_matrix(directed_breast_g_with_PAM50)

load("../GRACKLE_data/data/Breast/TCGA/Breast_filtered_gene_expression_with_PAM50.RData")


expression <- as.data.frame(t(expression_data))


expression <- expression %>%
  filter(rownames(.) %in% rownames(breast_subtype_metadata))

dim(breast_g_adjacency)
dim(expression)

results <- list()
for ( x in 1:100) {
  
    print(x)

    dat <- split_data(expression, breast_subtype_metadata , training_size = .7)
  
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
  
  k = 5
  
  pam50_clusters <- dat$test_metadata %>%
    rownames_to_column("sample") %>%
    pivot_longer(cols = -sample) %>%
    filter(value ==1) %>%
    select(-value)
  
  #    for(i in 1:nrow(grid_search)){

    scores <- mclapply(  1:nrow(grid_search), function(i) {
        print(i/nrow(grid_search))
        ## run GRACKLE NMF
        grackle <- GRACKLE(
            Y = dat$train_expression,
            net_similarity = as.matrix(breast_g_adjacency),
            patient_similarity = pat_sim,
            diff_threshold = 1e-5,
            lambda_1 = grid_search$lambda_1[i],
            lambda_2 =grid_search$lambda_2[i],
            k = k,
            verbose = F,
            beta = 0)
        ## correspondence between selected W LV's and top loading gene modules
        W_test <- project_W(dat$test_expression,grackle$H,k) 
        top <- apply(W_test,1, function(x) which(x == max(x)))
                                        # top <- lapply(top, function(x) ifelse(length(x) ==1, x, NA))
        both_clusters <- cbind(pam50_clusters,top = unlist(top))
        both_clusters <- both_clusters %>%
            filter(!is.na(top))
        score = ARI(as.factor(both_clusters$name), as.factor(both_clusters$top))
        return(score)
  },mc.cores = 1) 
     
  grid_search$score <- unlist(scores)

      
                  
 
    
     nmf <- runNMF(dat$train_expression,k, "lee", seed = 42)
     ## correspondence between selected W LV's and top loading gene modules
     W_test <- project_W(dat$test_expression,nmf$H,k) 
     top <- apply(W_test,1, function(x) which(x == max(x)))
     both_clusters <- cbind(pam50_clusters, top = unlist(top))
     
     nmf_res = ARI(as.factor(both_clusters$name), as.factor(both_clusters$top))
     nmf_res
     
     
     
     
      grnmf <- runGRNMF(expression_matrix = dat$train_expression, k = 5, seed = 42, max_iter = 200, alpha = 1)
      W_test <- project_W(dat$test_expression,grnmf$H,k) 
      top <- apply(W_test,1, function(x) which(x == max(x)))
      both_clusters <- cbind(pam50_clusters, top = unlist(top))
      
      grnmf_res = ARI(as.factor(both_clusters$name), as.factor(both_clusters$top))
      grnmf_res
      
      return(list(grid_search = grid_search, nmf_res = nmf_res, grnmf_res = grnmf_res))
}
      
      save(results, file = "../GRACKLE_results/TCGA_BRCA/TCGA_breast_with_PAM50_results.RData")
    
    
