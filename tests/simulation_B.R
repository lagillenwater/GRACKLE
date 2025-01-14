library(devtools)
load_all()

## ## load gene expression and network to model
## network <- read.csv("../GRACKLE_data/data/mammary_epithelium_top", sep = "\t", header = F)

## from <- entrezToHGNC(network$V1)
## to <- entrezToHGNC(network$V2)

## network <- network %>%
##     left_join(from, by = c(  "V1" = "entrez_id")) %>%
##     left_join(to, by = c(  "V2" = "entrez_id")) 
    

## network <- network %>%
##     filter(!(is.na(symbol.x) | is.na(symbol.y))) %>%
##     select(symbol.x, symbol.y,V3)

## names(network) <- c("from", "to", "prob")

## subnetwork <- network %>%
##     filter(prob > .9) %>%
##     mutate(prob = 1)

## subnetwork[is.na(subnetwork)] <- 0
## save(subnetwork, file = "../GRACKLE_data/breast_subnetwork.RData")

library(igraph)

load("../GRACKLE_data/breast_subnetwork.RData")

g <- graph_from_data_frame(subnetwork)

components <- decompose(g)
g <- components[[1]]

## generate list of permuted networks with distinctly permuted modules
set.seed(42)
clusters <- cluster_louvain(as.undirected(g), resolution = 1)
membership <- clusters$membership
membership_counts <- table(as.factor(membership))
gene_clusters <- 1:5

adj <- as.matrix(as_adjacency_matrix(g))

num_groups <-5
num_per_group <-100
n_row <- num_groups*num_per_group
sim_exp <- matrix(nrow =n_row , ncol = ncol(adj), runif(n_row*ncol(adj), min = 0, max = 20))
colnames(sim_exp) <- colnames(adj)
rownames(sim_exp) <- paste0("sample", 1:n_row)

## simulate up/down regulation of pathways
for(i in 1:num_groups) {
    end <- i * num_per_group
    start <- end - num_per_group + 1

    ## not working to have up and down regulation 1/13/25

        sim_exp[start:end, membership == gene_clusters[i]] <-     matrix(nrow =num_per_group ,ncol = sum(membership == gene_clusters[i]),runif(num_per_group*sum(membership == gene_clusters[i]), min =5 , max = 20))




}

dist_matrix <- dist(sim_exp)
hc <- hclust(dist_matrix)
#table(cutree(hc, k = 5))

gene_dist_matrix <- dist(t(sim_exp))
gene_hc <- hclust(gene_dist_matrix)
table(cutree(gene_hc, k = 24))
membership_counts

metadata <- simulateMetadata(group_size = rep(num_per_group,num_groups), group_labels = paste0("subgroup", 1:num_groups), row_names = rownames(sim_exp))


iterations <- 5
for( d in 1:iterations) {

    cat("\n")
    print(paste("Iteration: ", d))

    dat <- split_data(sim_exp, metadata , training_size = .7, seed = d)
    pat_sim <- as.matrix(dat$train_metadata) %*% t(dat$train_metadata)
    min_vals <- apply(dat$train_expression,2, min)
    max_vals <- apply(dat$train_expression,2, max)
    ## min-max scale the input matrix
    dat$train_expression <-as.matrix( min_max_scale(dat$train_expression,min_vals,max_vals))
    dat$test_expression <- as.matrix(min_max_scale(dat$test_expression,min_vals,max_vals))

    i_seq <-seq(0,1,.25)
    j_seq <-seq(0,1,.25)

    print("GRACKLE grid search ...") 
    grid_search <- as.data.frame(expand.grid(i_seq,j_seq))
    names(grid_search) <- c("lambda_1", "lambda_2")
    grid_search$score <- 0

##    scores <- list()
    scores <- mclapply(1:nrow(grid_search), function(i) {
        ##     ## run GRACKLE
        grackle <- GRACKLE(
            Y = dat$train_expression,
            net_similarity = adj,
            patient_similarity = pat_sim,
            diff_threshold = 1e-4,
            lambda_1 = grid_search$lambda_1[i],

            lambda_2 =  grid_search$lambda_2[i],
            k = num_groups,
            verbose = F,
            error_terms = F,
            beta = 0,
            iterations = 50)
        colnames(grackle$H) <- colnames(adj)
        ## correspondence between selected W LV's and top loading gene modules
        score <- evaluationWrapper(test_expression = dat$test_expression,
                                   test_metadata = dat$test_metadata,
                                   H_train = grackle$H,
                                   k = num_groups,
                                   clusters = clusters,
                                   aligned_clusters= 1:num_groups)
        return(score)
    }, mc.cores = 6)

    grid_search$score <- unlist(scores)
        
    print(paste("avg GRACKLE score", round(mean(grid_search$score),3)))
    print(paste("top GRACKLE score", max(grid_search %>% filter(scores == max(grid_search$score)) %>% .$score)))
    print(grid_search %>% filter(scores == max(grid_search$score)))
    print(paste("NMF score", round(grid_search[1,"score"],3)))

    ## nmf <- runNMF(dat$train_expression,num_groups, "lee", seed = 42)
    ## nmf_res <- evaluationWrapper(test_expression = dat$test_expression,
    ##                              test_metadata = dat$test_metadata,
    ##                              H_train = nmf$H,
    ##                              k = num_groups,
    ##                              clusters = clusters,
    ##                              aligned_clusters = large_clusters)


    grnmf <- runGRNMF(expression_matrix = dat$train_expression, k = num_groups, seed = 42, max = 50, alpha = 1)
    colnames(grnmf$H) <- colnames(adj)
    grnmf_res <- evaluationWrapper(test_expression = dat$test_expression,
                                   test_metadata = dat$test_metadata,
                                   H_train = grnmf$H,
                                   k = num_groups,
                                   clusters = clusters,
                                   aligned_clusters = 1:5)
    print(paste("GRNMF score", grnmf_res))
}
