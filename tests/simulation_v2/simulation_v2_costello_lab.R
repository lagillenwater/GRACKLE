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
## for running tensorflow
library(reticulate)
use_virtualenv("./grackle_env", required = TRUE)
library(tensorflow)


### load in the ppi network
load("./GRACKLE_data/STRING_long.RData")

## process weight
names(string_long)[3] <- "weight"  ## move to string processing file
string_long$weight <- string_long$weight/1000


## filter to a smaller graph object
sting_long <- string_long[1:1000,]

## converte to igrpah
new_g <- graph_from_edgelist(as.matrix(string_long[,1:2]), directed = F)
E(new_g)$weight <- string_long$weight

## create the simularity matrices for the algorithm
net_similarity <- as_adjacency_matrix(new_g, attr = "weight", sparse = F)
diag(net_similarity) <- 0

## set seed for reproducibility
set.seed(42)

## cluster the graph 
clusters <- cluster_louvain(new_g, resolution = 5)
membership <- clusters$membership
membership_counts <- table(membership)

## set the sample size, number of groups, number of genes
sample_n <- 1000
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
pathway_n <- 8## number of pathways to perturb in each group
for(i in 1:group_n) {
    group <- paste0("subgroup",i)
    samples <- metadata %>%
        filter(!!sym(group) == 1) %>%
        rownames()
    pathways <- sample(names(membership_counts)[membership_counts>6],pathway_n)
    for(j in pathways) {
        if(sample(1:2, 1) ==1) {
            tmp_exp <- matrix(nrow = length(samples), ncol = membership_counts[j], runif(length(samples) * membership_counts[j], min = 35,max = 50))
        } else {
            tmp_exp <- matrix(nrow = length(samples), ncol = membership_counts[j], runif(length(samples) * membership_counts[j], min = 0,max = 15))
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
print(prior_ari)
print(post_ari)


## run grackle
dat <- split_data(sim_exp, metadata , training_size = .7, seed = 10)

metadata_distance <- dist(dat$train_metadata)
metadata_similarity <- as.matrix(1/(as.matrix(metadata_distance) + 1))
patient_similarity <-  metadata_similarity
diag(patient_similarity) <- 0




min_vals <- apply(dat$train_expression,2, min)
max_vals <- apply(dat$train_expression,2, max)
## min-max scale the input matrix
dat$train_expression <-as.matrix( min_max_scale(dat$train_expression,min_vals,max_vals))
dat$test_expression <- as.matrix(min_max_scale(dat$test_expression,min_vals,max_vals))

grackle <- GRACKLE(
    Y = dat$train_expression,
    net_similarity = net_similarity,
    patient_similarity = patient_similarity,
    diff_threshold = 1e-6,
    lambda_1 = 100,
    lambda_2 = 0,
    k = group_n,
    iterations = 1000)


set.seed(42)
train_kmeans = kmeans(grackle$W,centers = group_n)$cluster
sample_labels <- dat$train_metadata %>% pivot_longer(cols = everything()) %>% filter(value == 1) %>% select(name)
ARI(train_kmeans, sample_labels$name)

W_test <- project_W(dat$test_expression,grackle$H,group_n)
test_kmeans <- kmeans(W_test,centers = group_n)$cluster
sample_labels <- dat$test_metadata %>% pivot_longer(cols = everything()) %>% filter(value == 1) %>% select(name)
ARI(test_kmeans, sample_labels$name)













