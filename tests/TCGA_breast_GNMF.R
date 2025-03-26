#setwd("~/OneDrive - The University of Colorado Denver/Projects/GRACKLE/")
library(tidyverse)
library(igraph)
library(devtools)
library(ggplot2)
library(parallel)
library(aricode)
library(SummarizedExperiment)
## library(reticulate)
## use_virtualenv("/mnt/grackle_env")
## library(tensorflow)
library(devtools)
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


## load the graph and expression data

load( file = "../GRACKLE_data/data/Breast/TCGA/directed_breast_g_with_PAM50.RData")
load( file = "../GRACKLE_data/data/Breast/TCGA/Breast_filtered_gene_expression_with_PAM50.RData")

breast_g_adjacency <- as_adjacency_matrix(directed_breast_g_with_PAM50)

expression <- as.data.frame(t(expression_data))
tokeep <- colnames(expression)[colnames(expression) %in% colnames(breast_g_adjacency)]
expression <- expression %>%
    filter(rownames(.) %in% rownames(breast_subtype_metadata)) %>%
    select(tokeep)


results <- list()
for ( x in 1:50) {
  
    print(paste("itertion",x))

  dat <- split_data(expression, breast_subtype_metadata , training_size = .7, seed = x)
  pat_sim <- as.matrix(dat$train_metadata) %*% t(dat$train_metadata)
  min_vals <- apply(dat$train_expression,2, min) + 1e-10
  max_vals <- apply(dat$train_expression,2, max) 
  ## min-max scale the input matrix
  dat$train_expression <-as.matrix( min_max_scale(dat$train_expression,min_vals,max_vals))
  dat$test_expression <- as.matrix(min_max_scale(dat$test_expression,min_vals,max_vals))

    k = ncol(breast_subtype_metadata)


    
  pam50_clusters <- dat$test_metadata %>%
    rownames_to_column("sample") %>%
      pivot_longer(cols = -sample)%>%
    filter(value ==1) %>%
    select(-value)

                    
        grnmf <- runGRNMF(expression_matrix = dat$train_expression, k = 5, seed = 42, max_iter = 100, alpha = a)
        W_test <- project_W(dat$test_expression,grnmf$H,k)
        top <- as.data.frame(apply(W_test,1, function(x) which(x == max(x)))) %>%
            rownames_to_column("sample")
        names(top)[2] <-  "top"
        both_clusters <- pam50_clusters %>%
            left_join(top, by = "sample")
        both_clusters <- both_clusters %>%
            filter(!is.na("top"))
        tmp_res= ARI(as.factor(both_clusters$name), as.factor(both_clusters$top))
        print(paste("GRNMF score for alpha", a, " = ", tmp_res))
        grnmf_res[grnmf_res$alpha == a, "score"]  = tmp_res
    }
        
    results[[x]] <- list(grnmf_res = grnmf_res)
    save(results, file = paste0("TCGA_breast_with_PAM50_results_GNMF_", Sys.Date(), ".RData"))

}



