

## Libraries

## For graph clustering
library(igraph)
## GRACKLE package
library(devtools)
load_all()
## statistical testing
library(aricode)
library(cluster)


### load in the ppi network
load("./GRACKLE_data/STRING_long.RData")
## transform weight to 0-1 scale
names(string_long)[3] <- "weight"  ## move to string processing file
string_long$weight <- string_long$weight/1000
## Convert to igraph
new_g <- graph_from_edgelist(as.matrix(string_long[,1:2]), directed = F)
E(new_g)$weight <- string_long$weight

## create the simularity matrices for the algorithm
## net_similarity <- as_adjacency_matrix(new_g, attr = "weight", sparse = F)
##save(net_similarity, file = "./GRACKLE_data/STRING_adjacency.RData")

## load the full adjacency
load("./GRACKLE_data/STRING_adjacency.RData")
## set seed for reproducibility
set.seed(42)
## cluster the graph 
clusters <- cluster_louvain(new_g, resolution = 5)
membership <- clusters$membership
names(membership) <- V(new_g)$name
membership_counts <- table(membership)

### Filter for testing
##modules_tokeep <- sample(names(membership_counts)[membership_counts > 100 & membership_counts < 250],15)
counter = 1
post_ari <- 0
post_gene_ari <- 0
sample_ari <- .7
gene_ari <- .7
threshold <- .05
sample_n <- 1000
group_n <- 5

metadata <- simulateMetadata(group_size = sample_n/group_n,  group_labels = paste0("subgroup", 1:group_n), row_names = paste0("s",1:sample_n))
metadata_labels <- metadata %>% pivot_longer(cols = everything()) %>% filter(value == 1) %>% select(name)

while((abs(post_ari- sample_ari) > threshold) | (abs(post_gene_ari -gene_ari) > threshold)) {
    set.seed(counter)
    modules_tokeep <- sample(names(membership_counts)[membership_counts > 50 & membership_counts < 300], 30)
    genes_tokeep<- names(membership)[which(membership %in% modules_tokeep)]
    net_similarity_filtered <- net_similarity[genes_tokeep,genes_tokeep]
    print(dim(net_similarity_filtered))
    ## simualate the metadata
    ##### simualte the gene expression data
    ## create a random background distribution
    genes_n <- ncol(net_similarity_filtered)
    sim_exp <- matrix(nrow = sample_n, ncol = genes_n, runif(sample_n * genes_n, min = 0, max = 50))
    rownames(sim_exp) <- rownames(metadata)
    colnames(sim_exp) <- colnames(net_similarity_filtered)
    ## Cluster the simulated expression matrix
    set.seed(42)
    prior_kmeans <- kmeans(sim_exp,centers = group_n, iter.max = 100, nstart = 3)$cluster
    ## find the ARI between the clusters and labels
    prior_ari = ARI(metadata_labels$name, prior_kmeans)
    ## same for eghe gene clusters
    set.seed(42)
    pre_gene_kmeans = kmeans(t(sim_exp), centers = 15, iter.max = 100, nstart = 5)$cluster
    ## pre_gene_kmeans = kmeans(t(sim_exp[,colnames(sim_exp) %in% names(membership)[membership %in% modules_tokeep]]), centers = 15, iter.max = 100, nstart = 5)$cluster
    sample_gene_membership <- membership[which(membership %in% modules_tokeep)]
    pre_gene_ari = ARI(sample_gene_membership, pre_gene_kmeans)
    ## simulate coordinated up and down regulation of genes in pathways within groups
    pathway_n <- 2
    ## number of pathways to perturb in each group
    all_pathways <- list()
    post_sim_exp <- sim_exp
    for(i in 1:group_n) {
        group <- paste0("subgroup",i)
        samples <- metadata %>%
            filter(!!sym(group) == 1) %>%
            rownames()
        ##    pathways <- sample(names(membership_counts)[membership_counts>6],pathway_n)
        set.seed(i)
        pathways <- sample(modules_tokeep,pathway_n)
        ##        print(pathways)
        all_pathways[[i]] <- pathways
        for(j in pathways) {
            if(sample(1:2, 1) ==1) {
                tmp_exp <- matrix(nrow = length(samples), ncol = membership_counts[j], runif(length(samples) * membership_counts[j], min = 20,max = 40))
            } else {
                tmp_exp <- matrix(nrow = length(samples), ncol = membership_counts[j], runif(length(samples) * membership_counts[j], min = 0,max = 20))
            }
        }
        pathway_genes <- names(membership)[membership == j]
        post_sim_exp[samples,pathway_genes] <- tmp_exp
    }
    ## Cluster the simulated expression matrix
    set.seed(42)
    post_kmeans <- kmeans(post_sim_exp,centers = group_n, iter.max = 100, nstart = 5)$cluster
    ## find the ARI between the clusters and labels
    post_ari = ARI(metadata_labels$name, post_kmeans)
    ## cluster the genes in the simulated data
    ## find the ARI between the gene clusters and modules
    all_pathways <- unlist(all_pathways)
    unique_pathways <- unique(all_pathways)
    set.seed(42)
    post_gene_kmeans = kmeans(t(post_sim_exp[,colnames(sim_exp) %in% names(membership)[membership %in% unique_pathways]]), centers = length(unique_pathways), iter.max = 100, nstart = 5)$cluster
    sample_gene_membership <- membership[which(membership %in% unique_pathways)]
    post_gene_ari = ARI(sample_gene_membership, post_gene_kmeans)
    print(prior_ari)
    print(post_ari)
    print(pre_gene_ari)
    print(post_gene_ari)
    counter = counter+1
}

simulated_data <- list('sim_exp' = post_sim_exp, 'unique_pathways' = unique_pathways, 'net_similarity_filtered' = net_similarity_filtered, 'metadata' = metadata, 'membership' = membership, 'membership_counts' = membership_counts)
save(simulated_data, file = "simulated_data_n_1000_ari_7.Rdata")

