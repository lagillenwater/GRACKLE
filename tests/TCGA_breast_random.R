#setwd("~/OneDrive - The University of Colorado Denver/Projects/GRACKLE/")
library(tidyverse)
library(igraph)
library(devtools)
library(ggplot2)
library(parallel)
library(aricode)
library(SummarizedExperiment)
library(reticulate)
use_virtualenv("/mnt/grackle_env")
library(tensorflow)
library(devtools)
load_all()

load("../GRACKLE_data/data/Breast/TCGA/Breast_metadata.RData")


##### PAM50
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


dim(breast_g_adjacency)
dim(expression)

load("TCGA_breast_with_PAM50_results.RData")

res <-list()
iterations <- 50
for(x in 1:iterations){

    dat <- split_data(expression, breast_subtype_metadata , training_size = .7, seed = x)
    min_vals <- apply(dat$train_expression,2, min) + 1e-10
    max_vals <- apply(dat$train_expression,2, max) 
    ## min-max scale the input matrix
    dat$train_expression <-as.matrix( min_max_scale(dat$train_expression,min_vals,max_vals))
    dat$test_expression <- as.matrix(min_max_scale(dat$test_expression,min_vals,max_vals))
    k = ncol(breast_subtype_metadata)

    num_rows_to_permute <- nrow(dat$train_metadata)
    set.seed(x) # For reproducibility
    rows_to_permute <- sample(rownames(dat$train_metadata), num_rows_to_permute)
    dat$train_metadata[rows_to_permute, ] <- dat$train_metadata[sample(rows_to_permute), ]
    permuted_pat_sim <- as.matrix(dat$train_metadata) %*% t(dat$train_metadata)

    pam50_clusters <- dat$test_metadata %>%
        rownames_to_column("sample") %>%
        pivot_longer(cols = -sample)%>%
        filter(value ==1) %>%
        select(-value)

    tmp <-results[[x]][[1]]
    l_1 <- tmp[which(tmp$score == max(tmp$score)),"lambda_1"][1]
    l_2 <- tmp[which(tmp$score == max(tmp$score)),"lambda_2"][1]
        
   grackle<- GRACKLE(
       Y = dat$train_expression,
       net_similarity = as.matrix(breast_g_adjacency),
       patient_similarity = permuted_pat_sim,
       diff_threshold = 1e-4,
       lambda_1 = l_1,
       lambda_2 =l_2,
       k = k,
       verbose = F,
       beta = 0,
       error_terms = F,
       iterations = 100)

    z = 100
    ## correspondence between selected W LV's and top loading gene modules
    while(any(is.na(grackle$H))) {
        z = z-5
        print(z)
        grackle <- GRACKLE(
            Y = dat$train_expression,
            net_similarity = as.matrix(breast_g_adjacency),
            patient_similarity = permuted_pat_sim,
            diff_threshold = 1e-4,
            lambda_1 = l_1,
            lambda_2 =l_2,
            k = k,
            verbose = F,
            beta = 0,
            error_terms = F,
            iterations = z)
    }
    
    W_test <- project_W(dat$test_expression,as.matrix(grackle$H),k)
    top <- as.data.frame(apply(W_test,1, function(x) which(x == max(x)))) %>%
        rownames_to_column("sample")
    names(top)[2] <-  "top"
    both_clusters <- pam50_clusters %>%
        left_join(top, by = "sample")
    both_clusters <- both_clusters %>%
        filter(!is.na("top"))
    grackle_res = ARI(as.factor(both_clusters$name), as.factor(both_clusters$top))

    res[[x]] <- grackle_res
}

    save(res, file = "TCGA_breast_random_patients_to_PAM50_results.RData")




##### Methylation
names(metadata) <- make.names(names(as.data.frame(metadata)))
breast_subtype_metadata <- metadata %>%
    as.data.frame() %>%
    select( paper_DNA.Methylation.Clusters) %>%
    drop_na( "paper_DNA.Methylation.Clusters" ) %>%
  rownames_to_column("ID") %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = "paper_DNA.Methylation.Clusters", values_from = value) %>%
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

dim(breast_g_adjacency)
dim(expression)

load("TCGA_breast_Methylation_to_PAM50_results.RData")

res <-list()
iterations <- 50
for(x in 1:iterations){
    print(x)
    dat <- split_data(expression, breast_subtype_metadata , training_size = .7, seed = x)
    min_vals <- apply(dat$train_expression,2, min) + 1e-10
    max_vals <- apply(dat$train_expression,2, max) 
    ## min-max scale the input matrix
    dat$train_expression <-as.matrix( min_max_scale(dat$train_expression,min_vals,max_vals))
    dat$test_expression <- as.matrix(min_max_scale(dat$test_expression,min_vals,max_vals))
    k = 5

    num_rows_to_permute <- nrow(dat$train_metadata)
    set.seed(x) # For reproducibility
    rows_to_permute <- sample(rownames(dat$train_metadata), num_rows_to_permute)
    dat$train_metadata[rows_to_permute, ] <- dat$train_metadata[sample(rows_to_permute), ]
    permuted_pat_sim <- as.matrix(dat$train_metadata) %*% t(dat$train_metadata)

    pam50_clusters <- dat$test_metadata %>%

    rownames_to_column("sample") %>%
        pivot_longer(cols = -sample)%>%
        filter(value ==1) %>%
        select(-value)

    tmp <-results[[x]][[1]]
    l_1 <- tmp[which(tmp$score == max(tmp$score)),"lambda_1"][1]
    l_2 <- tmp[which(tmp$score == max(tmp$score)),"lambda_2"][1]
        
   grackle<- GRACKLE(
       Y = dat$train_expression,
       net_similarity = as.matrix(breast_g_adjacency),
       patient_similarity = permuted_pat_sim,
       diff_threshold = 1e-4,
       lambda_1 = l_1,
       lambda_2 =l_2,
       k = k,
       verbose = F,
       beta = 0,
       error_terms = F,
       iterations = 100)

    z = 100
    ## correspondence between selected W LV's and top loading gene modules
    while(any(is.na(grackle$H))) {
        z = z-5
        print(z)
        grackle <- GRACKLE(
            Y = dat$train_expression,
            net_similarity = as.matrix(breast_g_adjacency),
            patient_similarity = permuted_pat_sim,
            diff_threshold = 1e-4,
            lambda_1 = l_1,
            lambda_2 =l_2,
            k = k,
            verbose = F,
            beta = 0,
            error_terms = F,
            iterations = z)
    }
    
    W_test <- project_W(dat$test_expression,as.matrix(grackle$H),k)
    top <- as.data.frame(apply(W_test,1, function(x) which(x == max(x)))) %>%
        rownames_to_column("sample")
    names(top)[2] <-  "top"
    both_clusters <- pam50_clusters %>%
        left_join(top, by = "sample")
    both_clusters <- both_clusters %>%
        filter(!is.na("top"))
    grackle_res = ARI(as.factor(both_clusters$name), as.factor(both_clusters$top))

    res[[x]] <- grackle_res
}

    save(res, file = "TCGA_breast_random_patients_Methylation_to_PAM50_results.RData")






##### Independent

breast_subtype_metadata <- metadata %>%
    as.data.frame() %>%
    select( paper_DNA.Methylation.Clusters, paper_lncRNA.Clusters, paper_Protein.Clusters, paper_CNV.Clusters) %>%
    rownames_to_column("ID") %>%
    mutate(value = 1) %>% pivot_wider(names_from = c(paper_DNA.Methylation.Clusters, paper_lncRNA.Clusters, paper_Protein.Clusters, paper_CNV.Clusters), values_from = value) %>%
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

dim(breast_g_adjacency)
dim(expression)

load("TCGA_breast_Independent_to_PAM50_results.RData")

res <-list()
iterations <- 50
for(x in 1:iterations){

    print(x)
    dat <- split_data(expression, breast_subtype_metadata , training_size = .7, seed = x)
    min_vals <- apply(dat$train_expression,2, min) + 1e-10
    max_vals <- apply(dat$train_expression,2, max) 
    ## min-max scale the input matrix
    dat$train_expression <-as.matrix( min_max_scale(dat$train_expression,min_vals,max_vals))
    dat$test_expression <- as.matrix(min_max_scale(dat$test_expression,min_vals,max_vals))
    k = 5

    num_rows_to_permute <- nrow(dat$train_metadata)
    set.seed(x) # For reproducibility
    rows_to_permute <- sample(rownames(dat$train_metadata), num_rows_to_permute)
    dat$train_metadata[rows_to_permute, ] <- dat$train_metadata[sample(rows_to_permute), ]
    permuted_pat_sim <- as.matrix(dat$train_metadata) %*% t(dat$train_metadata)

    pam50_clusters <- dat$test_metadata %>%
        rownames_to_column("sample") %>%
        pivot_longer(cols = -sample)%>%
        filter(value ==1) %>%
        select(-value)

    tmp <-results[[x]][[1]]
    l_1 <- tmp[which(tmp$score == max(tmp$score)),"lambda_1"][1]
    l_2 <- tmp[which(tmp$score == max(tmp$score)),"lambda_2"][1]
        
   grackle<- GRACKLE(
       Y = dat$train_expression,
       net_similarity = as.matrix(breast_g_adjacency),
       patient_similarity = permuted_pat_sim,
       diff_threshold = 1e-4,
       lambda_1 = l_1,
       lambda_2 =l_2,
       k = k,
       verbose = F,
       beta = 0,
       error_terms = F,
       iterations = 100)

    z = 100
    ## correspondence between selected W LV's and top loading gene modules
    while(any(is.na(grackle$H))) {
        z = z-5
        print(z)
        grackle <- GRACKLE(
            Y = dat$train_expression,
            net_similarity = as.matrix(breast_g_adjacency),
            patient_similarity = permuted_pat_sim,
            diff_threshold = 1e-4,
            lambda_1 = l_1,
            lambda_2 =l_2,
            k = k,
            verbose = F,
            beta = 0,
            error_terms = F,
            iterations = z)
    }
    
    W_test <- project_W(dat$test_expression,as.matrix(grackle$H),k)
    top <- as.data.frame(apply(W_test,1, function(x) which(x == max(x)))) %>%
        rownames_to_column("sample")
    names(top)[2] <-  "top"
    both_clusters <- pam50_clusters %>%
        left_join(top, by = "sample")
    both_clusters <- both_clusters %>%
        filter(!is.na("top"))
    grackle_res = ARI(as.factor(both_clusters$name), as.factor(both_clusters$top))

    res[[x]] <- grackle_res
}

    save(res, file = "TCGA_breast_random_patients_Independent_to_PAM50_results.RData")
