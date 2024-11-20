## Script for generating simulated gene expression and metadata based on a previously inferred network. 
library(sgnesR)
library(igraph)

## read in the breast adjacency matrix and expression data
breast_network <- read.csv("./data/Breast_Network.csv.csv",row.names = 1)

breast_expression <- read.csv("./data/Breast_Expression.csv", row.names = 1)

breast_gtex <- read.delim("./data/gene_reads_v10_breast_mammary_tissue.gct.gz", skip = 2)
rownames(breast_gtex) <- breast_gtex$Name



## convert ensembl ID's to hgnc ID's
#hgnc_ids <- ensemblToHGNC(names(breast))$symbol

## convert data.frame into directed, weighted igraph object
g_breast <- graph_from_data_frame(breast_network,directed = T)

