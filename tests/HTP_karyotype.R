#setwd("~/OneDrive - The University of Colorado Denver/Projects/GRACKLE/")
library(tidyverse)
library(igraph)
library(devtools)
library(ggplot2)
library(parallel)
library(aricode)
library(SummarizedExperiment)
library(reticulate)
use_virtualenv("/mnt/grackle_env")
library(tensorflow)
library(devtools)
library(FNN)
library(NEMO)
library(SNFtool)
library(cluster)
library(enrichR)
load_all()


comorbidities <- read.delim("./GRACKLE_data/P4C_Comorbidity_020921.tsv", skip= 1)
comorbidities <- comorbidities %>% pivot_wider(id_cols = RecordID, names_from = Condition, values_from = HasCondition)
names(comorbidities) <- make.names(names(comorbidities))
metadata <- read.delim("./GRACKLE_data/P4C_metadata_021921_Costello.txt")
metadata <- comorbidities %>% full_join(metadata, by = "RecordID")
htp_expression <- get(load("./GRACKLE_data/processed_transcriptome.RData"))$expression
##htp_expression <- get(load("processed_transcriptome_standardized.RData"))
##htp_expression <- exp(htp_expression)
variance <- apply(htp_expression,2,var)
tokeep <- variance[order(variance, decreasing = T)][1:5000]
htp_expression <- htp_expression%>%
     select(names(tokeep))

## load("../GRACKLE_data/data/blood_filtered_network.RData")

## directed_blood_network <- filtered_network_long %>%
##     filter(from %in% names(htp_expression) & to %in% names(htp_expression)) %>%
##     filter(probability > 1)

## blood_g <- graph_from_data_frame(directed_blood_network)
## blood_g_adjacency <- as_adjacency_matrix(blood_g)

### load in the ppi network
string <- read.delim("./GRACKLE_data/9606.protein.links.v12.0.txt", sep = " ")
string$score <- string$combined_score/1000
info <- read.delim("./GRACKLE_data/9606.protein.aliases.v12.0.txt", sep= "\t")
info_filt <- info %>%
    filter(source== "Ensembl_HGNC_symbol") %>%
    select(-source)
string_HGNC <- string %>%
    left_join(info_filt, by = c("protein1" = "X.string_protein_id"), relationship = "many-to-many") %>%
    left_join(info_filt, by = c("protein2" = "X.string_protein_id"), relationship = "many-to-many")
names(string_HGNC)[-c(1:4)] <- c("from","to")
string_filt <- string_HGNC %>%
    filter(from %in% names(htp_expression) & to %in% names(htp_expression)) %>%
    select(from,to,score)
string_adj <- string_filt %>%
    pivot_wider(id_cols = "from", names_from = "to", values_from = "score") %>%
    column_to_rownames("from")
htp_expression <- htp_expression %>%
    select(names(string_adj))
string_adj <- string_adj[match(names(htp_expression), rownames(string_adj)),match(names(htp_expression),names(string_adj))]
string_adj[is.na(string_adj)] <- 0


### 2025-02-18
### for testing why the typical pathways were underrepresented when separating based on other co-occurring conditions
metadata <- metadata %>% filter(LabID %in% rownames(htp_expression)) 
metadata <- metadata %>% select( c("Karyotype","Any.autoimmune.skin.condition"  ,"Any.hypothyroidism" , "Any.sleep.apnea","Obesity",     "LabID", "Anxiety", "Depression"))
## metadata <- metadata %>% filter(LabID %in% rownames(htp_expression)) %>% select(c("LabID", "Karyotype"))
## metadata <- metadata %>%
##     mutate(value = 1)%>%
##     pivot_wider(names_from = "Karyotype", values_from = "value")
## #### end test
metadata <- metadata %>% column_to_rownames("LabID")
metadata$Karyotype <- ifelse(metadata$Karyotype == "T21",1,0)
metadata[is.na(metadata)] <- 0
## IMS <- data.frame( IMS_cluster = get(load("all_clust.RData"))$clustering)
## IMS <- IMS %>%
##     rownames_to_column("sample")
## IMS <- IMS %>%
##     rownames_to_column("sample") %>%
htp_expression <- htp_expression[rownames(htp_expression) %in% rownames(metadata),]
htp_expression <-htp_expression[rownames(metadata),]
min_vals <- apply(htp_expression,2, min) 
max_vals <- apply(htp_expression,2, max) 
  ## min-max scale the input matrix
htp_expression <-as.matrix( min_max_scale(htp_expression,min_vals,max_vals))



## K = 25
## patient_knn <- get.knn(metadata, k = K)

## patient_nearest_neighbors <- patient_knn$nn.index
## patient_knn_distance <- patient_knn$nn.dist

## patient_mat <- matrix(.5, nrow = nrow(metadata), ncol = ncol(metadata))


## patient_distance <- (dist2(as.matrix(metadata),as.matrix(metadata)))^(1/2)
## patient_similarity <- affinityMatrix(patient_distance, K =K, sigma = .5)

patient_distance <- dist(metadata,diag = F)

#### comparing euclidean to jaccard distance
## library(vegan)
## tmp_metadata <- metadata + 1
## patient_jaccard <- vegdist(tmp_metadata, method = "jaccard",na.rm = T)
## euc_jac_cor <- diag(cor(as.matrix(patient_distance), as.matrix(patient_jaccard)))


patient_similarity <- as.matrix(1/(1+as.matrix(patient_distance)))


net_similarity <- as.matrix(string_adj)
diag(net_similarity) <- 1


### extracting specific resutls
k=4
Y = htp_expression
net_similarity = net_similarity
patient_similarity = patient_similarity
diff_threshold = 1e-5
lambda_1 = 1
lambda_2 =1
verbose = F
beta = 0
error_terms = F

dbs <- c(  "GO_Biological_Process_2023")

for(k in 10) {
    grackle <- GRACKLE(
        Y = htp_expression,
        net_similarity = net_similarity,
        patient_similarity = patient_similarity,
        diff_threshold = 1e-6,
        lambda_1 = 1,
        lambda_2 =1,
        k = k)
    W_eval <- cbind(metadata,grackle$W)
    condition_scores <- sapply(names(metadata), function(y){
        sapply(1:k, function(x){
            eval(parse(text = paste0("m0 <- wilcox.test(LV",x,"~as.factor(",y,"),W_eval)")))
            ## print(m0)
            p <- m0$p.value
            ##b <- m0$estimate[2] - m0$estimate[1]
            ##p <- ifelse(b>0,p,-1*p)
            return(p)
            ##qnorm(p/2)/sqrt(nrow(W_eval))
        })
    })
    sig_conditions<- apply(condition_scores,1,function(x) {
        sig <- which(abs(x) < .05)
        sum_sig <- ifelse(x[sig] > 0,"+++","---")
        names_sig <- names(sig)
        paste(c(sum_sig,names_sig),collapse  = " ")
    })
    names(sig_conditions) <- paste0("LV",1:k)
    sig_conditions <- sig_conditions[sig_conditions != ""]
    print(paste("for K ", k, " significant:" ,length(sig_conditions)))
#    print(apply(grackle$H, 1, function(x) cor(grackle$H[1,], x)))
    colnames(grackle$H) <- colnames(htp_expression)
    rownames(grackle$H) <- paste0("LV",1:k)
    ## enrichments <- lapply( names(sig_conditions), function(rank){
    ##     gene_loadings <- grackle$H[rank,]
    ##     cutoff <- floor(.05 * ncol(htp_expression))
    ##     input_genes <- names(sort(gene_loadings, decreasing = T)[1:cutoff])
    ##     enriched <- enrichr(input_genes, dbs, background = names(htp_expression))
    ##     unlist(lapply(enriched, function(x) x %>% filter(Adjusted.P.value < .01) %>% select(Term)))
    ## })
    ## names(enrichments) <- names(sig_conditions)
}



inner_scores <- lapply( names(metadata), function(col) {
    tmp <- as.data.frame(table(top,metadata[[col]]) ) %>%
        pivot_wider(id_cols = top, names_from = Var2, values_from = Freq) %>%
        column_to_rownames("top")
    return(tmp)
})
names(inner_scores) <- names(metadata)



save(inner_scores, file = "HTP_k_6_l_1_100_l_2_2.RData")

 
write.csv(as.data.frame(t(grackle$H)), file = "grackle_H_HTP_k_4.csv")
