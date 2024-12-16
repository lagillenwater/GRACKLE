### Simulation study A

## set seed
set.seed(42)
#setwd("~/OneDrive - The University of Colorado Denver/Projects/GRACKLE/")
library(tidyverse)
library(igraph)
library(devtools)
load_all()

# # Assign random edge weights based on the distribution of in and out degree over regulatory relationships in GRAND Breast network
breast_expression <- read.csv("./data/Breast_expression.csv")
breast_network <- read.csv("./data/Breast_network.csv")

## pivot longer
breast_network_long <- breast_network %>%
  pivot_longer(cols = -X)

# ## Remove interactions below a threshold of 1 This threshold is based on the evidence of interaction, regardless of directionality.
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
  dplyr::select(X,symbol,value) %>%
  rename(from = X,to = symbol,probability = value)

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
  filter(from %in% names(filtered_breast_expression) & to %in% names(filtered_breast_expression))

save(filtered_breast_expression, file = "./data/Breast/network_filtered_breast_expression.RData ")

#Calculate the correlation coefficients between pairs of genes to determine
directed_breast_network <- makeDirected(filtered_breast_expression,filtered_breast_network_long)
save(directed_breast_network, file = "./data/Breast/directed_breast_network.RData")

## Filter by correlation for simulation experiments
load("./data/Breast/directed_breast_network.RData")
directed_breast_network <- directed_breast_network %>% filter(abs(correlation) > .2)

# some specific modifications for sgnesR. Expression simulator only takes 1 and -1 weights. 
directed_breast_network$weight <- ifelse(directed_breast_network$weight > 0,1,-1)
directed_breast_network <- directed_breast_network %>%
  select(from,to,weight)

##
load("./data/Breast/network_filtered_breast_expression.RData")
correlations_filtered_breast_expression <- filtered_breast_expression %>%
  dplyr::select( union(directed_breast_network$from, directed_breast_network$to))
save(correlations_filtered_breast_expression, file = "./data/Breast/correlation_filtered_breast_expression.RData")

# ## Turn into igraph
breast_g <- graph_from_data_frame(directed_breast_network)
save(breast_g, file = "./data/Breast/directed_breast_igraph.RData")



