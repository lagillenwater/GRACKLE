### Simulation study A

## set seed
set.seed(42)
setwd("~/OneDrive - The University of Colorado Denver/Projects/GRACKLE/")
library(tidyverse)
library(igraph)
# test ~/OneDrive - The University of Colorado Denver/Projects/GRACKLE/data/# test <- read.csv('./data/Breast_Mammary_Tissue_TF_BindTFDF.csv')
# tfnet <- read.delim("./data/TF-RISet.tsv")
# 
# trrust <- read.delim("./data/trrust_rawdata.human.tsv", header = F)
# trrust <- trrust %>%
#   mutate(weight = ifelse(V3 == "Repression",-1,
#                          ifelse(V3 == "Activation",1,0))) %>%
#   rename(from = V1, to = V2) %>%
#   select(c(from,to,weight)) %>%
#   filter(weight != 0)
# 

# Assign random edge weights based on the distribution of in and out degree over regulatory relationships in GRAND Breast network
breast_expression <- read.csv("./data/Breast_expression.csv")
breast_network <- read.csv("./data/Breast_network.csv")

# ## pivot longer
breast_network_long <- breast_network %>%
  pivot_longer(cols = -X)

## Remove interactions below a threshold of 1 This threshold is based on the evidence of interaction, regardless of directionality.
breast_network_long <- breast_network_long %>%
  filter(value >= 1)

# change ensembl ID's to HGNC ID's
network_symbol_map <- ensemblToHGNC(breast_network_long$name)
expression_symbol_map <- ensemblToHGNC(breast_expression$X)

breast_network_long <- breast_network_long %>%
  left_join(network_symbol_map, by = c("name" = "ensembl_gene_id"), relationship = "many-to-many")

breast_network_long <- breast_network_long %>%
  filter(!is.na(symbol))

breast_network_long <- breast_network_long %>%
  select(X,symbol,value) %>%
  rename(from = X,to = symbol,weight = value)


breast_expression <- breast_expression %>%
  left_join(expression_symbol_map, by = c("X" = "ensembl_gene_id"), relationship = "many-to-many")

breast_expression <- breast_expression %>%
  filter(!is.na(symbol)) %>%
  select(-X)

## filter expression and network 
unique_genes <- unique(union(breast_network_long$from,breast_network_long$to))
filtered_breast_expression <- breast_expression %>%
  filter(symbol %in% unique_genes) %>%
  column_to_rownames('symbol') %>%
  t() %>%
  as.data.frame()

filtered_breast_network_long <- breast_network_long %>%
  filter(from %in% rownames(filtered_breast_expression) & to %in% rownames(filtered_breast_expression))


###### TODO:threshold the directionality based on significance/rho value

## Calculate the correlation coefficients between pairs of genes to determine
# directed_breast_network <- makeDirected(filtered_breast_expression,filtered_breast_network_long)
# save(directed_breast_network, file = "directed_breast_network.RData")
load("directed_breast_network.RData")
directed_breast_network <- directed_breast_network %>% filter(abs(correlation) > .85)

## remove self loops 
directed_breast_network <- directed_breast_network %>% filter(from != to)



## Find network properties
breast_g <- graph_from_data_frame(directed_breast_network)
save(breast_g, file = "directed_breast_igraph.RData")
load("directed_breast_igraph.RData")
## in and out degree of breast network
in_degrees <- degree(breast_g,mode = "in")
out_degrees <- degree(breast_g,mode = "out")
degrees <- degree(breast_g)


# generate random breast_g# generate random nework 
g <- degree.sequence.game(out_degrees, in_degrees, method = "simple")
sum(directed_breast_network$direction == 1)/nrow(directed_breast_network)
sm <- sample(c(1,-1), ecount(g), rep=TRUE, p=c(.51,.49))
E(g)$op <- sm
g
is_connected(g)
components <- decompose(g)

g <- components[[which.max(sapply(components, vcount))]]
is_connected(g)
## permute modules
# generate list of permuted networks with distinctly permuted modules

###### TODO: identify modules and then extract subgraph with distinct and overlapping modules

clusters <- cluster_walktrap(g)
membership <- cut_at(clusters,no = 100)
membership_counts <- table(as.factor(membership))
large_clusters <- membership_counts[order(membership_counts, decreasing = T)][1:5]

permuted_networks <- lapply(large_clusters, function(x) modulePermutations(g,membership,x))

# simulate gene expression for the permuted networks
sim_exp <- lapply(permuted_networks, function(x) simulateExpression(x, iterations = 3, max_expression = 3000, num_samples = 5))

## combine expression data
sim_exp <- lapply(sim_exp, as.data.frame)
sim_exp <- lapply(sim_exp, function(x) {
  x %>% rownames_to_column("gene")
})
combined_exp <- sim_exp %>% reduce(full_join, by = "gene")

combined_exp <- combined_exp %>%
  column_to_rownames("gene")
  

names(combined_exp) <- paste0("sample", 1:ncol(combined_exp))

combined_exp <- as.data.frame(t(combined_exp))

## TODO check the quantiles of the simulated expression data against observed expression data
# breast_exp_quantiles <- quantile(as.matrix(breast_expression))
# lapply(sim_exp, quantile)



metadata <- as.data.frame(c(rep(1, 5), rep(2, 5), rep(3, 5), rep(4, 5), rep(5, 5)))
names(metadata) <- "label"
rownames(metadata) <- rownames(combined_exp)

## split the data into train and test
dat <- split_data(combined_exp, metadata , training_size = .7)

set.seed(42)


  
  pat_sim <- similarityCalc(as.matrix(dat$train_metadata))
  # pat_sim[pat_sim < 1] <- 0
  
  ## min-max scale the input matrix
  Y <- apply(dat$train_expression, 2, min_max_scale)
  
  i_seq <-seq(0,1,.1) 
  j_seq <-seq(0,1,.1)
  
  grid_search <- as.data.frame(expand.grid(i_seq,j_seq))
  names(grid_search) <- c("lambda_1", "lambda_2")
  grid_search$score <- 0
  
  for(i in 1:nrow(grid_search)){
    ## run GRACKLE NMF
    g_res <- GRACKLE(
      Y = Y,
      net_similarity = net_sim,
      patient_similarity = pat_sim,
      diff_threshold = 1e-6,
      lambda_1 = grid_search$lambda_1[i],
      lambda_2 = grid_search$lambda_2[i],
      k = 5, 
      verbose = F,
      beta = 0)
    
    ## project the test data into W_test
    W_test <- project_W(apply(dat$test_expression,2,min_max_scale), g_res$H)
    
    ## evaluate sample loadings
    top_sample_LVs <- sampleLoadingsEvaluation(W_test,dat$test_metadata)
    
    ## evaluate gene loadings
    top_loadings <- geneLoadingsEvaluation(g_res$H, 10)
    
    ## correspondence between selected W LV's and top loading gene modules
    grid_search$score[i] <- sum(unlist(lapply(1:5, function(x) {
      identical(top_sample_LVs$top[x], top_loadings$top[x])
    })))/5
  }
  
  print( ggplot(grid_search, aes(x = lambda_1, y= lambda_2, fill = score)) +
           geom_tile() +
           ggtitle(x))
  
  
  
  
  l_res <- nmf(Y,5,'lee', seed = 42)
  
  W_test <-  project_W(apply(dat$test_expression,2,min_max_scale), coef(l_res))
  ## evaluate sample loadings
  top_sample_LVs <- sampleLoadingsEvaluation(W_test,dat$test_metadata)
  
  ## evaluate gene loadings
  top_loadings <- geneLoadingsEvaluation(coef(l_res), 10)
  
  ## correspondence between selected W LV's and top loading gene modules
  l_corr <- sum(unlist(lapply(1:5, function(x) {
    identical(top_sample_LVs$top[x], top_loadings$top[x])
  })))/5
  
  return(list(g = g_corr, l = l_corr))
})  


toplot <- data.frame(noise = seq(.1,.9,.2),  nmf = as.numeric(res[2,]))
library(tidyverse)
toplot <- toplot %>%
  pivot_longer(cols = -noise)

ggplot(toplot, aes(x = noise, y = value, color = name)) +
  ylim(c(0,1))+
  geom_line(linewidth = 2) +
  theme_classic()

