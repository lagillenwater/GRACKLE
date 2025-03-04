########
## This is a script for simulation studies to evaluate the GRACKLE algorithm.
#######

## Libraries

## For graph clustering
library(igraph)

## GRACKLE package
library(devtools)
load_all()

## statistical testing
library(aricode)

### load in the ppi network
load("./GRACKLE_data/STRING_long.RData")
ls()

## covert to igraph object
names(string_long)[3] <- "weight"  ## move to string processing file
new_g <- graph_from_data_frame(string_long, directed = F)

## set seed for reproducibility
set.seed(42)

## cluster the graph 
clusters <- cluster_louvain(new_g, resolution = 3)
membership <- clusters$membership
membership_counts <- table(membership)

## set the sample size, number of groups, number of genes
sample_n <- 100
group_n <- 5
genes_n <- length(V(new_g))


## simualate the metadata
metadata <- simulateMetadata(group_size = sample_n/group_n,  group_labels = paste0("subgroup", 1:group_n), row_names = paste0("s",1:sample_n))
metadata_labels <- metadata %>% pivot_longer(cols = everything()) %>% filter(value == 1) %>% select(name)



##### simualte the gene expression data

## NB: write a function for this process

## create a random background distribution
sim_exp <- matrix(nrow = sample_n, ncol = genes_n, runif(sample_n * genes_n, min = 0, max = 40))
rownames(sim_exp) <- rownames(metadata)

## Cluster the simulated expression matrix
set.seed(42)
prior_kmeans <- kmeans(sim_exp,centers = group_n)$cluster

## find the ARI between the clusters and labels
prior_ari = ARI(metadata_labels$name, prior_kmeans)

## simulate coordinated up and down regulation of genes in pathways within groups
pathway_n <- 5 ## number of pathways to perturb in each group

for(i in 1:group_n) {
    group <- paste0("subgroup",i)
    samples <- metadata %>%
        filter(!!sym(group) == 1) %>%
        rownames()
    pathways <- sample(names(membership_counts),pathway_n)
    for(j in pathways) {
        if(sample(1:2, 1) ==1) {
            tmp_exp <- matrix(nrow = length(samples), ncol = membership_counts[j], runif(length(samples) * membership_counts[j], min = 30,max = 50))
        } else {
            tmp_exp <- matrix(nrow = length(samples), ncol = membership_counts[j], runif(length(samples) * membership_counts[j], min = 0,max = 20))
        }
    }
    pathway_genes <- which(membership == j)
    sim_exp[samples,pathway_genes] <- tmp_exp
}


            
## Cluster the simulated expression matrix
set.seed(42)
post_kmeans <- kmeans(sim_exp,centers = group_n)$cluster

## find the ARI between the clusters and labels
post_ari = ARI(metadata_labels$name, post_kmeans)
                              

## create the simularity matrices for the algorithm
patient_distance <- dist(metadata,diag = F)
patient_similarity <- as.matrix(1/(1+as.matrix(patient_distance)))

net_similarity <- as.








