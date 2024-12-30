## Script for loading a processing TCGA data for GRACKLE

BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)

# Query for clinical data
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Clinical", 
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")

# Download the data
GDCdownload(query)

# Prepare the data
clinical_data <- GDCprepare(query)


# Extract subtype information
subtype_data <- clinical_data$clinical_patient_brca

save(subtype_data, file = "../GRACKLE_data/data/Breast/TCGA/Breast_subtype.RData")
# View the first few rows of the subtype data
head(subtype_data)


query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts")

# Download the data
GDCdownload(query)

# Prepare the datap p
gene_expression_data <- GDCprepare(query)

save(gene_expression_data,file = "../GRACKLE_data/data/Breast/TCGA/Breast_gene_expression.RData" )

# View the metadata
colData(gene_expression_data)

metadata <- colData(gene_expression_data)
save(metadata, file = "../GRACKLE_data/data/Breast/TCGA/Breast_metadata.RData")
# View the gene expression data
expression_data <- assay(gene_expression_data, "tpm_unstrand")

load("../GRACKLE_data/data/Breast/TCGA/Breast_filtered_gene_expression_with_PAM50.RData")

library(tidyverse)

expression_data <- expression_data %>%
  as.data.frame() %>%
  mutate(ensembl_ids = gsub("\\..*", "", rownames(.)))

hgnc_map <- ensemblToHGNC(expression_data$ensembl_ids)

expression_data <- expression_data %>%
  filter(ensembl_ids %in% hgnc_map$ensembl_gene_id) %>%
  left_join(hgnc_map, by = c("ensembl_ids" = "ensembl_gene_id"),relationship = "many-to-many") %>%
  distinct(symbol, .keep_all = T) %>%
  column_to_rownames("symbol") %>%
  dplyr::select(-ensembl_ids)



load("../GRACKLE_data/data/Breast/directed_breast_network.RData")  




## filter by variance
variance <- apply(expression_data,1,var)
tokeep <- variance[order(variance, decreasing = T)][1:2000]
length(tokeep)
expression_data <- expression_data %>%
  filter(rownames(.) %in% names(tokeep))



directed_breast_network <- directed_breast_network %>%
  filter(from %in% rownames(expression_data) & to %in% rownames(expression_data))

library(igraph)

directed_breast_g_with_PAM50 <- graph_from_data_frame(directed_breast_network %>% dplyr::select(from,to,weight), directed = TRUE)
save(directed_breast_g_with_PAM50, file = "../GRACKLE_data/data/Breast/TCGA/directed_breast_g_with_PAM50.RData")

expression_data <- expression_data %>%
  filter(rownames(.) %in% union(directed_breast_network$from,directed_breast_network$to))


save(expression_data, file = "../GRACKLE_data/data/Breast/TCGA/Breast_filtered_gene_expression_with_PAM50.RData")


## filtering out PAM50 genes
pam50 <- read.delim("./data/Breast/TCGA/pam50_annotation.txt")
directed_breast_network_without_PAM50 <- directed_breast_network %>%
  filter(!(from %in% pam50$GeneName | to %in% pam50$GeneName))

directed_breast_g_without_PAM50 <- graph_from_data_frame(directed_breast_network_without_PAM50 %>% dplyr::select(from,to,weight), directed = TRUE)
save(directed_breast_g_without_PAM50, file = "./data/Breast/TCGA/directed_breast_g_without_PAM50.RData")

expression_data_without_PAM50 <- expression_data %>%
  filter(!(rownames(.) %in% pam50$GeneName))
save(expression_data_without_PAM50, file = "./data/Breast/TCGA/Breast_filtered_gene_expression_without_PAM50.RData")         

