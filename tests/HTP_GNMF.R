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
library(NEMO)
library(SNFtool)
library(cluster)
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
length(tokeep)
htp_expression <- htp_expression%>%
    select(names(tokeep))



load("../GRACKLE_data/data/blood_filtered_network.RData")

directed_blood_network <- filtered_network_long %>%
    filter(from %in% names(htp_expression) & to %in% names(htp_expression)) %>%
    filter(probability > 1)


blood_g <- graph_from_data_frame(directed_blood_network, directed = T)
blood_g_adjacency <- as_adjacency_matrix(blood_g)


metadata <- metadata %>% filter(LabID %in% rownames(htp_expression)) %>%
    filter(Karyotype == "T21")
metadata <- metadata %>% select( c("Any.autoimmune.skin.condition"  ,"Any.hypothyroidism" , "Any.sleep.apnea","Obesity",     "LabID", "Anxiety", "Depression"))
metadata <- metadata %>% column_to_rownames("LabID")
metadata[is.na(metadata)] <- 0


## IMS <- data.frame( IMS_cluster = get(load("all_clust.RData"))$clustering)
## IMS <- IMS %>%
##     rownames_to_column("sample")
## IMS <- IMS %>%
##     rownames_to_column("sample") %>%

htp_expression <- htp_expression[rownames(htp_expression) %in% rownames(metadata),]
htp_expression <-htp_expression[rownames(metadata),]
htp_expression <- htp_expression %>%
    select(colnames(blood_g_adjacency))

dim(blood_g_adjacency)
dim(htp_expression)
dim(metadata)
identical(colnames(blood_g_adjacency), names(htp_expression))
identical(rownames(htp_expression), rownames(metadata))

min_vals <- apply(htp_expression,2, min)
max_vals <- apply(htp_expression,2, max)
## min-max scale the input matrix
htp_expression <-as.matrix( min_max_scale(htp_expression,min_vals,max_vals))



patient_similarity <- tcrossprod(as.matrix(metadata))
net_similarity <- as.matrix(blood_g_adjacency)



for(k in seq(2,4,1)) {
        res <- list()
    for(d in 1:50){

        print(paste("K",k,"iteration",d))
        dat <- split_data(htp_expression, metadata , training_size = .7, seed = d)
        pat_sim <- as.matrix(dat$train_metadata) %*% t(dat$train_metadata)
        min_vals <- apply(dat$train_expression,2, min) + 1e-10
        max_vals <- apply(dat$train_expression,2, max)
        ## min-max scale the input matrix
        dat$train_expression <-as.matrix( min_max_scale(dat$train_expression,min_vals,max_vals))
        dat$test_expression <- as.matrix(min_max_scale(dat$test_expression,min_vals,max_vals))

        grnmf_res <- data.frame(alpha = c(10^ 0, 10^1, 10^2, 10^3, 10^4), score = 0)
        ##grnmf_res <- data.frame(alpha = seq(.1), score = 0)
        for (a in grnmf_res$alpha) { 
    

            grnmf <- runGRNMF(expression_matrix = dat$train_expression, k = k, seed = 42, max_iter = 100, alpha = a)
            rownames(grnmf$W) <-  rownames(dat$train_expression)
            top <- as.data.frame(apply(grnmf$W,1, function(x) which(x == max(x)))) %>%
                rownames_to_column("sample")
            names(top)[2] <-  "top"
            inner_scores <- lapply( names(dat$train_metadata), function(col) {
                tmp <- data.frame(table(top = top$top,col = dat$train_metadata[[col]]) ) %>%
                    pivot_wider(id_cols = top, names_from = col, values_from = Freq) %>%
                    column_to_rownames("top")

                score <- fisher.test(tmp, alternative = "greater",simulate.p.value = TRUE)$p.value
                ## print(col)
                ## print(score)
                ## print(tmp)
                return(score)
            })
            
            tmp_res <- sum(unlist(inner_scores))
            print(paste("GRNMF score for alpha", a, " = ", tmp_res))
            grnmf_res[grnmf_res$alpha == a, "score"]  = tmp_res
        }
        res[[d]] <- list(grnmf_res = grnmf_res)
        
        save(res, file = paste0("HTP_results_GNMF_k_",k,"_", Sys.Date(), ".RData"))
    }
}





### extracting specific resutls
## k = 6
## l_1 = 100, l_2 = 2
k=4

Y = htp_expression
net_similarity = net_similarity
patient_similarity = patient_similarity
diff_threshold = 1e-4
lambda_1 = 100
lambda_2 =20
k = k
verbose = F
beta = 0
error_terms = F



grackle <- GRACKLE(
    Y = htp_expression,
    net_similarity = net_similarity,
    patient_similarity = patient_similarity,
    diff_threshold = 1e-4,
    lambda_1 = 100,
    lambda_2 =20,
    k = k,
    verbose = F,
    beta = 0,
    error_terms = F,
    iterations = 50)

top <- apply(grackle$W,1, function(x) which(x == max(x)))
inner_scores <- lapply( names(metadata), function(col) {
    tmp <- as.data.frame(table(top,metadata[[col]]) ) %>%
        pivot_wider(id_cols = top, names_from = Var2, values_from = Freq) %>%
        column_to_rownames("top")
    return(tmp)
})
names(inner_scores) <- names(metadata)
cor(grackle$H[1,], grackle$H[2,])

save(inner_scores, file = "HTP_k_6_l_1_100_l_2_2.RData")

colnames(grackle$H) <- colnames(htp_expression)
write.csv(as.data.frame(t(grackle$H)), file = "grackle_H_HTP_k_4.csv")
