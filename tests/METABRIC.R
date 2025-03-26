
library(tidyverse)

load("../GRACKLE_data/data/Breast/directed_breast_network.RData")  

expression_data <- read.table("./GRACKLE_data/METABRIC/brca_metabric/data_mrna_illumina_microarray.txt", header = TRUE)

expression_data <- expression_data %>%
    distinct(Hugo_Symbol, .keep_all = TRUE)%>%
    column_to_rownames("Hugo_Symbol") %>%
    select(-Entrez_Gene_Id)

dim(expression_data)



## filter by 
variance <- apply(expression_data,1,var)
tokeep <- variance[order(variance, decreasing = T)][1:5000]
length(tokeep)
expression_data <- expression_data %>%
  filter(rownames(.) %in% names(tokeep))


directed_breast_network <- directed_breast_network %>%
  filter(from %in% rownames(expression_data) & to %in% rownames(expression_data))

library(igraph)

directed_breast_g_with_PAM50 <- graph_from_data_frame(directed_breast_network %>% dplyr::select(from,to,weight), directed = TRUE)
save(directed_breast_g_with_PAM50, file = "./GRACKLE_data/METABRIC/directed_breast_g_with_PAM50.RData")

expression_data <- expression_data %>%
  filter(rownames(.) %in% union(directed_breast_network$from,directed_breast_network$to))


save(expression_data, file = "./GRACKLE_data/METABRIC/Breast_filtered_gene_expression_with_PAM50.RData")



