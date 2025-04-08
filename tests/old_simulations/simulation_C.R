
## load hallmark gene sets
library(msigdbr)
pathway_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")

## Filter for unique gene sets
library(tidyverse)
library(dplyr)
library(stringr)
library(tidyr)
library(devtools)
load_all()
library(tensorflow)
install_tensorflow(envname = "/mnt/grackle")

path_counts <- table(pathway_gene_sets$gs_name)
## create a gene by gene set binary matrix

gene_pathway_matrix <- pathway_gene_sets %>%
    mutate(present = 1) %>%
    pivot_wider(id_cols = c(gene_symbol), names_from = gs_name, values_from = present, values_fill = 0, values_fn = max) %>%
    column_to_rownames("gene_symbol")
    

## create gene by gene matrix
gene_sim <- tcrossprod(as.matrix(gene_pathway_matrix))
diag(gene_sim) <- 1

# gene_sim <- similarityCalc(as.matrix(reactome_matrix))

## create simulated data
num_groups <-5
num_per_group <-20
n_row <- (num_groups+1)*num_per_group

pathways <- sample(colnames(gene_pathway_matrix), num_groups)
gene_groups <- lapply(1:num_groups, function(x) {
    tokeep <- which(gene_pathway_matrix[[pathways[x]]] == 1)
    tmp <- gene_pathway_matrix[tokeep,]
    return(rownames(tmp))
})
gene_groups[[num_groups+1]] <- rownames(gene_pathway_matrix)[!(rownames(gene_pathway_matrix) %in% unlist(gene_groups))]

pairs <- combn(seq_along(gene_groups), 2, simplify = FALSE)
pairwise_intersections <- lapply(pairs, function(pair) {
    intersect(gene_groups[[pair[1]]], gene_groups[[pair[2]]])
})
shared_genes <- unlist(pairwise_intersections)

sim_exp <- matrix(nrow =n_row , ncol = ncol(gene_sim), runif(n_row*ncol(gene_sim), min = 0, max = 25))
colnames(sim_exp) <- colnames(gene_sim)
rownames(sim_exp) <- paste0("sample", 1:n_row)

## simulate up/down regulation of pathways
for(i in 1:num_groups) {
     end <- i * num_per_group
    start <- end - num_per_group + 1
     genes <- gene_groups[[i]][!(gene_groups[[i]] %in% shared_genes)]
     sim_exp[start:end, genes] <-     matrix(nrow =num_per_group ,ncol = length(genes),runif(num_per_group*length(genes), min =25 , max = 30 ))
}    


gene_dist_matrix <- dist(t(sim_exp))
gene_hc <- hclust(gene_dist_matrix)
hc_gene_groups <- cutree(gene_hc, k =6)
print(table(hc_gene_groups))
print(unlist(lapply(gene_groups,length)))

for(x in unique(hc_gene_groups)) {
    tmp <- names(hc_gene_groups)[hc_gene_groups==x]
    overlap = lapply(gene_groups, function(y) {
        sum(tmp %in% y)
    })
    print(unlist(overlap))
}

dist_matrix <- dist(sim_exp)
hc <- hclust(dist_matrix)
print(table(cutree(hc, k = num_groups+1)))


metadata <- simulateMetadata(group_size = rep(num_per_group,(num_groups+1)), group_labels = paste0("subgroup", 1:(num_groups+1)), row_names = rownames(sim_exp))


iterations <- 5
for(d in  1:iterations) {

    print(paste("Iteration", d))
    dat <- split_data(sim_exp, metadata , training_size = .7, seed = d)
    pat_sim <- as.matrix(dat$train_metadata) %*% t(dat$train_metadata)
    min_vals <- apply(dat$train_expression,2, min)
    max_vals <- apply(dat$train_expression,2, max)
    ## min-max scale the input matrix
    dat$train_expression <-as.matrix( min_max_scale(dat$train_expression,min_vals,max_vals))
    dat$test_expression <- as.matrix(min_max_scale(dat$test_expression,min_vals,max_vals))
    i_seq <-seq(0,1,.2)
    j_seq <-seq(0,1,.2)
    print("GRACKLE grid search ...") 
    grid_search <- as.data.frame(expand.grid(i_seq,j_seq))
    names(grid_search) <- c("lambda_1", "lambda_2")
    grid_search$score <- 0

    scores <- lapply( 1: nrow(grid_search),function(i) {
        ##     ## run GRACKLE
        ##        print(i)

        grackle <- GRACKLE(
            Y = dat$train_expression,
            net_similarity = gene_sim,
            patient_similarity = pat_sim,
            diff_threshold = 1e-4,
            lambda_1 = grid_search$lambda_1[i],
            lambda_2 =  grid_search$lambda_2[i],
            k = num_groups+1,
            verbose = F,
            error_terms = F,
            beta = 0,
            iterations = 20)
        
        colnames(grackle$H) <- colnames(gene_sim)
        ## correspondence between selected W LV's and top loading gene modules

        score <- evaluationWrapper(test_expression = dat$test_expression,
                                   test_metadata = dat$test_metadata,
                                   H_train = grackle$H,
                                   k = num_groups+1,
                                   clusters = gene_groups,
                                   aligned_clusters= 1:(num_groups+1))

        scores[[i]] <- score
     ##   print(score)
        
    })

    grid_search$score <- unlist(scores)
        
    print(paste("avg GRACKLE score", round(mean(grid_search$score),3)))
    print(paste("top GRACKLE score", max(grid_search %>% filter(scores == max(grid_search$score)) %>% .$score)))
    print(grid_search %>% filter(scores == max(grid_search$score)))
    print(paste("NMF score", round(grid_search[1,"score"],3)))

    grnmf <- runGRNMF(expression_matrix = dat$train_expression, k = (num_groups+1), seed = 42, max_ = 200, alpha = 1)
    colnames(grnmf$H) <- colnames(gene_sim)
    grnmf_res <- evaluationWrapper(test_expression = dat$test_expression,
                                   test_metadata = dat$test_metadata,
                                   H_train = grnmf$H,
                                   k = (num_groups + 1),
                                   clusters = gene_groups,
                                   aligned_clusters = 1:(num_groups+1))
    print(paste("GRNMF score", grnmf_res))
}



