## set seed
set.seed(42)
#setwd("~/OneDrive - The University of Colorado Denver/Projects/GRACKLE/")
library(tidyverse)
library(igraph)
library(devtools)
library(ggplot2)
library(ComplexHeatmap)
library(parallel)
library(aricode)
library(umap)
library(gridExtra)
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


expression <- expression %>%
  filter(rownames(.) %in% rownames(breast_subtype_metadata))
## results <- mclapply(1:10, function(y) {
    
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


    
#    for(i in 1:nrow(grid_search)){
## scores <- mclapply(  1:nrow(grid_search), function(i) {
      #   print(i/nrow(grid_search))
      
      ## run GRACKLE NMF
    i = 61
    grid_search[i,]
      
    grackle <- GRACKLE(
        Y = dat$train_expression,
        net_similarity = as.matrix(breast_g_adjacency),
        patient_similarity = pat_sim,
        diff_threshold = 1e-5,
        lambda_1 = grid_search$lambda_1[i],
        lambda_2 =grid_search$lambda_2[i],
        k = 5,
        verbose = T,
        beta = 0)
      
      
      ## correspondence between selected W LV's and top loading gene modules
      W_test <- project_W(dat$test_expression,grackle$H,k) 
      kmeans_res <- kmeans(W_test, centers = 5)
      tcga_clusters <- as.data.frame(kmeans_res$cluster) %>%
        rownames_to_column("sample")
      names(tcga_clusters)[2] <- "GRACKLE"
      
      pam50_clusters <- dat$test_metadata %>%
        rownames_to_column("sample") %>%
        pivot_longer(cols = -sample) %>%
        filter(value ==1)
  
      both_clusters <- pam50_clusters %>%
        left_join(tcga_clusters, by= "sample")
      
      score = ARI(as.factor(both_clusters$name), as.factor(both_clusters$GRACKLE))
      
      toplot <- umap(W_test)$layout
      colnames(toplot) <- c("umap1", "umap2")
      toplot <- toplot %>%
        as.data.frame() %>%
        rownames_to_column("sample") %>%
        left_join(both_clusters, by = "sample")
      
      p1<- ggplot(toplot, aes(x = umap1, y = umap2, color = name)) +
        geom_point()
      p2<- ggplot(toplot, aes(x = umap1, y = umap2, color = GRACKLE)) +
        geom_point()
      
      ggsave(filename = paste0("./results/TCGA/grackle_", grid_search$lambda_1[i],"_",  grid_search$lambda_2[i],".pdf"), grid.arrange(p1,p2))
      
                  
##        return(score)
#    }
    ## }, mc.cores = 10)   

    
##    grid_search$score <- unlist(scores)
 
    
     nmf <- runNMF(dat$train_expression,5, "lee", seed = 42)
     ## correspondence between selected W LV's and top loading gene modules
     W_test <- project_W(dat$test_expression,nmf$H,k) 
     kmeans_res <- kmeans(W_test, centers = 5)
     tcga_clusters <- as.data.frame(kmeans_res$cluster) %>%
       rownames_to_column("sample")
     names(tcga_clusters)[2] <- "kmeans"
     
     pam50_clusters <- dat$test_metadata %>%
       rownames_to_column("sample") %>%
       pivot_longer(cols = -sample) %>%
       filter(value ==1)
     
     both_clusters <- pam50_clusters %>%
       left_join(tcga_clusters, by= "sample")
     
     nmf_res = ARI(as.factor(both_clusters$name), as.factor(both_clusters$kmeans))
    
      toplot <- umap(W_test)$layout
     colnames(toplot) <- c("umap1", "umap2")
     toplot <- toplot %>%
       as.data.frame() %>%
       rownames_to_column("sample") %>%
       left_join(both_clusters, by = "sample")
     
     p1<- ggplot(toplot, aes(x = umap1, y = umap2, color = name)) +
       geom_point()
     p2<- ggplot(toplot, aes(x = umap1, y = umap2, color = kmeans)) +
       geom_point()
     
     ggsave(filename = "./results/TCGA/nmf.pdf", grid.arrange(p1,p2))    

     
     
     
      grnmf <- runGRNMF(expression_matrix = dat$train_expression, k = 5, seed = 42, max_iter = 200, alpha = .1)
    
      W_test <- project_W(dat$test_expression,grnmf$H,k) 
      kmeans_res <- kmeans(W_test, centers = 5)
      tcga_clusters <- as.data.frame(kmeans_res$cluster) %>%
        rownames_to_column("sample") 
     names(tcga_clusters)[2] <- "kmeans"
     
     pam50_clusters <- dat$test_metadata %>%
       rownames_to_column("sample") %>%
       pivot_longer(cols = -sample) %>%
       filter(value ==1)
     
     both_clusters <- pam50_clusters %>%
       left_join(tcga_clusters, by= "sample")
     
     grnmf_res = ARI(as.factor(both_clusters$name), as.factor(both_clusters$kmeans))
    
     toplot <- umap(W_test)$layout
     colnames(toplot) <- c("umap1", "umap2")
     toplot <- toplot %>%
       as.data.frame() %>%
       rownames_to_column("sample") %>%
       left_join(both_clusters, by = "sample")
     
     p1<- ggplot(toplot, aes(x = umap1, y = umap2, color = name)) +
       geom_point()
     p2<- ggplot(toplot, aes(x = umap1, y = umap2, color = kmeans)) +
       geom_point()
     
     ggsave(filename = "./results/TCGA/grnmf.pdf", grid.arrange(p1,p2))
    
     return(list(grid_search = grid_search, nmf_res = nmf_res, grnmf_res = grnmf_res))
##   }, mc.cores = 3)
  
  
  

## save(results, file = "./results/simulations/data/TCGA_results.RData")


