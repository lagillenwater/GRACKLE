

library(igraph)

## load GRN
breast_g <- get(load("../GRACKLE_data/data/breast_igraph_prob_1_cor_0_05.RData"))


## load breast metadata
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



breast_subtype_metadata <- as.matrix(breast_subtype_metadata)

pat_sim <- breast_subtype_metadata %*% t(breast_subtype_metadata)

pat_g <- graph_from_adjacency_matrix(pat_sim)

## calculate the entorpy
entropy  <- function(g) {
    degree_dist <- degree_distribution(g)
    entropy <- -sum(degree_dist * log(degree_dist, base =2), na.rm = T)
    return(entropy)
}

pat_degree_dist <- degree_distribution(pat_g)
grn_degree_dist <- degree_distribution(breast_g)

pat_entropy <- -sum(pat_degree_dist * log(pat_degree_dist, base =2), na.rm = T)
grn_entropy <- -sum(grn_degree_dist * log(grn_degree_dist, base =2), na.rm = T)

entropy(pat_g)



load("./GRACKLE_data/STRING_long.RData")
17 ## transform weight to 0-1 scale
18 names(string_long)[3] <- "weight"  ## move to string processing file
19 string_long$weight <- string_long$weight/1000
20 ## Convert to igraph
21 new_g <- graph_from_edgelist(as.matrix(string_long[,1:2]), directed = F)
22 E(new_g)$weight <- string_long$weight



## calcuate the modularity of graphs

transitivity(breast_g)
transitivity(pat_g)

grn_clusters <- cluster_louvain(as.undirected(breast_g),weights = NA)
pat_clusters <- cluster_louvain(as.undirected(pat_g), weights = NA)

modularity(grn_clusters)
modularity(pat_clusters)
