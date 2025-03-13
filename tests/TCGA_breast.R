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
   i_seq <-seq(0,100,10)
   j_seq <-seq(0,20,2)
   grid_search <- as.data.frame(expand.grid(i_seq,j_seq))
   names(grid_search) <- c("lambda_1", "lambda_2")
   grid_search$score <- 0
    total <- ncol(grid_search) * nrow(grid_search)
    k = ncol(breast_subtype_metadata)

    dim(dat$train_expression)
    dim(dat$test_expression)    
    
  pam50_clusters <- dat$test_metadata %>%
    rownames_to_column("sample") %>%
      pivot_longer(cols = -sample)%>%
    filter(value ==1) %>%
    select(-value)

    scores <- lapply( 1:nrow(grid_search), function(i) {
        ## run GRACKLE NMF

        grackle <- GRACKLE(
            Y = dat$train_expression,
            net_similarity = as.matrix(breast_g_adjacency),
            patient_similarity = pat_sim,
            diff_threshold = 1e-4,
            lambda_1 = grid_search$lambda_1[i],
            lambda_2 =grid_search$lambda_2[i],
            k = k,
            verbose = F,
            beta = 0,
            error_terms = F,
            iterations = 100)
        z = 100
        ## correspondence between selected W LV's and top loading gene modules
        while(any(is.na(grackle$H))) {
            z = z-5
        grackle <- GRACKLE(
            Y = dat$train_expression,
            net_similarity = as.matrix(breast_g_adjacency),
            patient_similarity = pat_sim,
            diff_threshold = 1e-4,
            lambda_1 = grid_search$lambda_1[i],
            lambda_2 =grid_search$lambda_2[i],
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
        score = ARI(as.factor(both_clusters$name), as.factor(both_clusters$top))

        
            return(score)
   })

    grid_search$score <- unlist(scores)

    ## nmf <- runNMF(dat$train_expression,k, "lee", seed = 42)
     ## ## correspondence between selected W LV's and top loading gene modules
     ## W_test <- project_W(dat$test_expression,nmf$H,k) 
     ## top <- apply(W_test,1, function(x) which(x == max(x)))
     ## both_clusters <- cbind(pam50_clusters, top = unlist(top))
     
     ## nmf_res = ARI(as.factor(both_clusters$name), as.factor(both_clusters$top))
     ## nmf_res

    print(paste("avg GRACKLE score", round(mean(grid_search$score),3)))
    print(paste("top GRACKLE score", max(grid_search %>% filter(scores == max(grid_search$score)) %>% .$score)))
    print(grid_search %>% filter(scores == max(grid_search$score)) )
    
    nmf_res <- grid_search[1,"score"]

    
    
    grnmf <- runGRNMF(expression_matrix = dat$train_expression, k = 5, seed = 42, max_iter = 100, alpha = 1)
    W_test <- project_W(dat$test_expression,grnmf$H,k)
    top <- as.data.frame(apply(W_test,1, function(x) which(x == max(x)))) %>%
        rownames_to_column("sample")
    names(top)[2] <-  "top"
    both_clusters <- pam50_clusters %>%
        left_join(top, by = "sample")
    both_clusters <- both_clusters %>%
        filter(!is.na("top"))
    grnmf_res = ARI(as.factor(both_clusters$name), as.factor(both_clusters$top))

    print(paste("NMF score", nmf_res))
    print(paste("GRNMF score", grnmf_res))
    
    results[[x]] <- list(grid_search = grid_search, nmf_res = nmf_res, grnmf_res = grnmf_res)



}

save(results, file = "TCGA_breast_with_PAM50_results.RData")



#### Full data
pat_sim <- tcrossprod(as.matrix(breast_subtype_metadata))
min_vals <- apply(expression,2, min) + 1e-10
max_vals <- apply(expression,2, max) 
## min-max scale the input matrix
train_expression  <-as.matrix( min_max_scale(expression,min_vals,max_vals))

k = ncol(breast_subtype_metadata)

    
pam50_clusters <- breast_subtype_metadata %>%
    rownames_to_column("sample") %>%
    pivot_longer(cols = -sample)%>%
    filter(value ==1) %>%
    select(-value)




grackle <- GRACKLE(
    Y = train_expression,
    net_similarity = as.matrix(breast_g_adjacency),
    patient_similarity = pat_sim,
    diff_threshold = 1e-4,
    lambda_1 = 50,
    lambda_2 =16,
    k = k,
    verbose = F,
    beta = 0,
    error_terms = F,
    iterations = 100)



rownames(grackle$W) <- rownames(train_expression)
top <- as.data.frame(apply(grackle$W,1, function(x) which(x == max(x)))) %>%
    rownames_to_column("sample")
names(top)[2] <-  "top"

both_clusters <- pam50_clusters %>%
    left_join(top, by = "sample")

both_clusters <- both_clusters %>%
    filter(!is.na("top"))

table(both_clusters$name, both_clusters$top)

pam50 <- read.delim("pam50_annotation.txt")

pam50LVs <- grackle$H
rownames(pam50LVs) <- c("Basal", "none", "LumA_Normal", "Her2", "LumB")
colnames(pam50LVs) <- colnames(train_expression)
    
write.csv(as.data.frame(t(pam50LVs)), file = "PAM50_gene_LVs.csv")

score = ARI(as.factor(both_clusters$name), as.factor(both_clusters$top))
