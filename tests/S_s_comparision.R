library(tidyverse)
library(devtools)
library(aricode)
load_all()


load("../GRACKLE_data/data/Breast/TCGA/Breast_metadata.RData")
names(metadata) <- make.names(names(as.data.frame(metadata)))

breast_subtype_metadata <- metadata %>%
  as.data.frame() %>%
  dplyr::select(paper_BRCA_Subtype_PAM50) %>%
  drop_na(paper_BRCA_Subtype_PAM50) %>%
  rownames_to_column("ID") %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = paper_BRCA_Subtype_PAM50, values_from = value) %>%
    column_to_rownames("ID") 
breast_subtype_metadata[is.na(breast_subtype_metadata)] <- 0






methylation_subtype_metadata <- metadata %>%
    as.data.frame() %>%
    select( paper_DNA.Methylation.Clusters) %>%
    drop_na( "paper_DNA.Methylation.Clusters" ) %>%
  rownames_to_column("ID") %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = "paper_DNA.Methylation.Clusters", values_from = value) %>%
    column_to_rownames("ID") 
methylation_subtype_metadata[is.na(methylation_subtype_metadata)] <- 0


independent_subtype_metadata <- metadata %>%
    as.data.frame() %>%
    select( paper_DNA.Methylation.Clusters, paper_lncRNA.Clusters, paper_Protein.Clusters, paper_CNV.Clusters) %>%
    rownames_to_column("ID") %>%
    mutate(value = 1) %>% pivot_wider(names_from = c(paper_DNA.Methylation.Clusters, paper_lncRNA.Clusters, paper_Protein.Clusters, paper_CNV.Clusters), values_from = value) %>%
    column_to_rownames("ID") 
independent_subtype_metadata[is.na(independent_subtype_metadata)] <- 0



breast_S <- as.matrix(breast_subtype_metadata) %*% t(breast_subtype_metadata)
methylation_S <- as.matrix(methylation_subtype_metadata) %*% t(methylation_subtype_metadata)
independent_S <- as.matrix(independent_subtype_metadata) %*% t(independent_subtype_metadata)


breast_k <- data.frame(breast = kmeans(breast_S, centers = 5)$cluster) %>% rownames_to_column("ID")
methylation_k <- data.frame(methylation = kmeans(methylation_S, centers = 5)$cluster) %>% rownames_to_column("ID")
independent_k <- data.frame(independent = kmeans(independent_S, centers = 5)$cluster) %>% rownames_to_column("ID")

k_clusters <- breast_k %>%
    left_join(methylation_k, by = "ID") %>%
    left_join(independent_k, by = "ID")

ARI(k_clusters$breast, k_clusters$methylation)

ARI(k_clusters$breast, k_clusters$independent)

ARI(k_clusters$independent, k_clusters$methylation)
