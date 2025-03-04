
### This is a script for downloading and processing the STRING network for use in the GRACKLE algorithm.


## load the necessary libraries
library(tidyverse)

### load the current version of string
string <- read.delim("./GRACKLE_data/9606.protein.links.v12.0.txt", sep = " ")

## load the protein metadata
info <- read.delim("./GRACKLE_data/9606.protein.aliases.v12.0.txt", sep= "\t")
info_filt <- info %>%
    filter(source== "Ensembl_HGNC_symbol") %>%
    select(-source)


## map ensembl names to HGNC symbols
string_HGNC <- string %>%
    left_join(info_filt, by = c("protein1" = "X.string_protein_id"), relationship = "many-to-many") %>%
    left_join(info_filt, by = c("protein2" = "X.string_protein_id"), relationship = "many-to-many")
names(string_HGNC)[-c(1:3)] <- c("from","to")

string_long <- string_HGNC %>%
    select(from,to,combined_score) %>%
    filter(!is.na(from) & !is.na(to))
    


### save the full data
save(string_long, file = "./GRACKLE_data/STRING_long.RData")


### subsample the adjacency for testing
