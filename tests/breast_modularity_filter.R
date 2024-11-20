## Script for generating simulated gene expression and metadata based on a previously inferred network. 
library(sgnesR)
library(igraph)

## read in the breast cancer adjacency matrix from 
breast <- read.csv("./data/Breast.csv",row.names = 1)

## convert ensembl ID's to hgnc ID's
#hgnc_ids <- ensemblToHGNC(names(breast))$symbol

## convert data.frame into directed, weighted igraph object
g <- graph_from_data_frame(breast,directed = T)

