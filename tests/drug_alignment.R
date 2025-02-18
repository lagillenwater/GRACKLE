setwd("~/OneDrive - The University of Colorado Denver/Projects/GRACKLE/")
library(tidyverse)
library(devtools)


drugs <- read.delim("/Users/lucas/OneDrive - The University of Colorado Denver/Projects/GRACKLE/data/repurposing_drugs_20180907.txt", skip = 9)

drugs_filtered <- drugs %>%
  separate_rows(target, sep = "\\|") %>%
  distinct(target,moa)


loadings <- read.csv("/Users/lucas/OneDrive - The University of Colorado Denver/Projects/GRACKLE/results/PAM50_gene_LVs.csv")
load("../GRACKLE_data/data/Breast/directed_breast_network.RData")

var_loadings <- loadings %>%
  rowwise() %>%
  mutate(variance = var(c_across(-X))) %>%
  arrange(-variance)


drugs_filtered <- drugs_filtered %>%
  filter(target %in% loadings$X)

moas <- table(drugs_filtered$moa)
moas <- moas[names(moas)!=""]

moa_scores <- sapply(2:6, function(y) {
  
  lv <- loadings[,c('X',names(loadings)[y])]
  
  scores<-  sapply(names(moas), function(x) {
    gene_set <- drugs_filtered %>% 
      filter(moa == x) %>%
      .$target
    
    lv_scores <- lv %>% 
      filter(X %in% gene_set) 
    
    mean(lv_scores[,2])
        
  })
})

colnames(moa_scores) <- colnames(loadings)[2:6]

var_moa_scores <- moa_scores %>%
  as.data.frame() %>%
  rownames_to_column("moa") %>%
  rowwise() %>%
  mutate(variance = var(c_across(-moa))) %>%
  distinct(variance, .keep_all = T) %>%
  arrange(-variance) %>%
  column_to_rownames("moa") %>%
  select(-variance)

library(ComplexHeatmap)


ht <- Heatmap(as.matrix(var_moa_scores[1:20,]),
              heatmap_legend_param = list(direction = "horizontal"),
              row_names_side = "left",
              row_dend_side = "right",
              row_names_gp = gpar(fontsize = 12),
              width = ncol(var_moa_scores)*unit(5, "mm")
)
pdf("./results/PAM50_drug_heatmap.pdf")
draw(ht, heatmap_legend_side = "bottom")
dev.off()


load("../GRACKLE_data/data/Breast/TCGA/Breast_metadata.RData")

breast_subtype_metadata <- metadata %>%
  as.data.frame() %>%
  dplyr::select(paper_BRCA_Subtype_PAM50) %>%
    drop_na(paper_BRCA_Subtype_PAM50) 






