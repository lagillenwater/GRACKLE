########
## This is a script for simulation studies to evaluate the GRACKLE algorithm.
#######

## Libraries

## For graph clustering
library(igraph)
## GRACKLE package
library(devtools)
load_all()
## statistical testing
library(aricode)
library(cluster)
## for running tensorflow
library(reticulate)
use_virtualenv("/mnt/grackle_env")
library(tensorflow)


load("simulated_data_n_1000_ari_7.Rdata")


results <- list()
counter <- 1
group_n <- 5
for(d in 1:50) {

    print(d)
    ## run grackle
    dat <- split_data(simulated_data$sim_exp, simulated_data$metadata , training_size = .7, seed = d)
    metadata_distance <- dist(dat$train_metadata)
    patient_similarity <- as.matrix(1/(as.matrix(metadata_distance) + 1))
    ##patient_similarity <-  metadata_similarity
    diag(patient_similarity) <- 0
    min_vals <- apply(dat$train_expression,2, min)
    max_vals <- apply(dat$train_expression,2, max)
    ## min-max scale the input matrix
    dat$train_expression <-as.matrix( min_max_scale(dat$train_expression,min_vals,max_vals))
    dat$test_expression <- as.matrix(min_max_scale(dat$test_expression,min_vals,max_vals))

    ## patient clustering of data
    set.seed(42)
    train_kmeans = kmeans(dat$train_expression, centers = group_n, iter.max = 100, nstart =3)$cluster
    sample_labels <- dat$train_metadata %>% pivot_longer(cols = everything()) %>% filter(value == 1) %>% select(name)
    ##     print(paste("training data ARI:", ARI(train_kmeans, sample_labels$name)))
    train_ari= ARI(train_kmeans, sample_labels$name)

    set.seed(42)
    test_kmeans = kmeans(dat$test_expression, centers = group_n, iter.max = 100, nstart =3)$cluster
    sample_labels <- dat$test_metadata %>% pivot_longer(cols = everything()) %>% filter(value == 1) %>% select(name)
    ## print(paste("testing data ARI:", ARI(test_kmeans, sample_labels$name)))
    test_ari =ARI(test_kmeans, sample_labels$name)

    ## gene clustering
    set.seed(42)
    train_gene_kmeans = kmeans(t(dat$train_expression[,colnames(simulated_data$sim_exp) %in% names(simulated_data$membership)[simulated_data$membership %in% simulated_data$unique_pathways]]), centers = length(simulated_data$unique_pathways))$cluster
    sample_gene_membership <- simulated_data$membership[which(simulated_data$membership %in% simulated_data$unique_pathways)]
    ## print(paste("training gene ARI:", ARI(sample_gene_membership, train_gene_kmeans)))
    train_gene_ari <- ARI(sample_gene_membership, train_gene_kmeans)

    k = 10

    ## Perform SVD
    svd_result <- svd(dat$train_expression)
    ##find variance explianed
    ## singular_values <- svd_result$d
    ## variance_explained <- (singular_values^2)/sum(singular_values^2)
    ## cumulative_variance_explained <- cumsum(variance_explained)
    ## k <- max(which(cumulative_variance_explained < .9))
    ## Initialize W and H matrices using SVD
    svd_W <- abs(svd_result$u[,1:k] %*% diag(svd_result$d[1:k]))
    svd_H <- t(abs(svd_result$v[,1:k]))
    svd_W <- min_max_scale(svd_W, min_vals = apply(svd_W,2,min), max_vals = apply(svd_W,2,max)) 
    svd_H <- t(min_max_scale(t(svd_H), min_vals = apply(svd_H,1,min), max_vals = apply(svd_H,1,max)))

    set.seed(42)
    svd_kmeans = kmeans(svd_W,centers = group_n, iter.max = 100, nstart =5)$cluster
    sample_labels <- dat$train_metadata %>% pivot_longer(cols = everything()) %>% filter(value == 1) %>% select(name)
    ## print(paste("SVD ARI:",ARI(svd_kmeans, sample_labels$name)))
    svd_ari = ARI(svd_kmeans, sample_labels$name)

    set.seed(42)
    svd_gene_kmeans = kmeans(t(svd_H[,which(colnames(simulated_data$sim_exp) %in% names(simulated_data$membership)[simulated_data$membership %in% simulated_data$unique_pathways])]), centers = length(simulated_data$unique_pathways), iter.max = 100, nstart =5)$cluster
    sample_gene_membership <- simulated_data$membership[which(simulated_data$membership %in% simulated_data$unique_pathways)]
    ## print(paste("SVD gene ARI:", ARI(sample_gene_membership, svd_gene_kmeans)))
    svd_gene_ari = ARI(sample_gene_membership, svd_gene_kmeans)

    ## Y = dat$train_expression
    ## net_similarity = simulated_data$net_similarity_filtered
    ## patient_similarity = patient_similarity
    ## diff_threshold = 1e-6
    ## lambda_1 = 1
    ## lambda_2 = 0
    ## k = 15
    ## iterations = 1000
    ## verbose = T
    ## svd_W = svd_W
    ## svd_H = svd_H

    ## NB: split SVD out from grackle so as not to repeat with multiple iterations
    for(i in c(0,.01,.1,.5)) {
        for(j in c(0,.01,.1, .5)) {
            
            ## print(paste(i,j))        
            grackle <- GRACKLE(
                Y = dat$train_expression,
                net_similarity = simulated_data$net_similarity_filtered,
                patient_similarity = patient_similarity,
                diff_threshold = 1e-5,
                lambda_1 =i,
                lambda_2 = j,
                k = k,
                iterations = 500,
                verbose = F,
                svd_W = svd_W,
                svd_H = svd_H)

            ## patient clustering of results
            set.seed(42)
            train_W_kmeans = kmeans(grackle$W,centers = group_n, iter.max = 100, nstart = 5)$cluster
            ## silhouette_score <- mean(silhouette(train_W_kmeans, dist(grackle$W))[,3])
            ## print(paste("silhouette score", silhouette_score))
            sample_labels <- dat$train_metadata %>% pivot_longer(cols = everything()) %>% filter(value == 1) %>% select(name)
            ## print(paste("train W ARI:", ARI(train_W_kmeans, sample_labels$name)))
            train_w_ari = ARI(train_W_kmeans, sample_labels$name)
            
            set.seed(42)
            W_test <- project_W(dat$test_expression,grackle$H)
            test_W_kmeans <- kmeans(W_test,centers = group_n, iter.max = 100, nstart =5)$cluster
            sample_labels <- dat$test_metadata %>% pivot_longer(cols = everything()) %>% filter(value == 1) %>% select(name)
            ## print(paste("test W ARI:", ARI(test_W_kmeans, sample_labels$name)))
            test_w_ari = ARI(test_W_kmeans, sample_labels$name)

            ## gene clusering results
            set.seed(42)
            train_H_kmeans = kmeans(t(grackle$H[,which(colnames(simulated_data$sim_exp) %in% names(simulated_data$membership)[simulated_data$membership %in% simulated_data$unique_pathways])]), centers = length(simulated_data$unique_pathways),iter.max = 100, nstart =5)$cluster
            ## silhouette_score <- mean(silhouette(train_H_kmeans, dist(t(grackle$H[,which(colnames(simulated_data$sim_exp) %in% names(membership)[membership %in% simulated_data$unique_pathways])])))[,3])
            ## print(paste("silhouette score", silhouette_score))
            ## print(paste("H ARI", ARI(sample_gene_membership, train_H_kmeans)))
            train_h_ari= ARI(sample_gene_membership, train_H_kmeans)

            ## set.seed(42)
            ## H_test <- project_H(dat$test_expression,grackle$W)
            ## test_H_kmeans <- kmeans(t(H_test[,which(colnames(simulated_data$sim_exp) %in% names(membership)[membership %in% simulated_data$unique_pathways])]), centers = length(simulated_data$unique_pathways), iter.max = 100, nstart =5)$cluster
            ## print(paste("test H ARI:", ARI(sample_gene_membership, test_H_kmeans)))

            ##cat("\n")
            results[[counter]] <- c(i,j,train_ari, test_ari, train_gene_ari, train_w_ari, test_w_ari, train_h_ari, svd_ari, svd_gene_ari)
            counter = counter+1
        }
    }
}

results <- do.call(rbind,results)
colnames(results) <- c('lambda_1', 'lambda_2',"train_ari", 'test_ari', 'train_gene_ari', 'train_w_ari', 'test_w_ari', 'train_h_ari', 'svd_ari', 'svd_gene_ari')


results_summary <- results %>%
    as.data.frame() %>%
    group_by(lambda_1,lambda_2) %>%
    mutate(mean_train_w_ari = mean(train_w_ari),
           mean_test_w_ari = mean(test_w_ari),
           mean_train_h_ari = mean(train_h_ari)
           ) %>%
    distinct(lambda_1,lambda_2, .keep_all = T) %>%
    select(-c(train_w_ari, test_w_ari,train_h_ari))



save(results_summary, file = "simulation_v2_results_summary_3.RData")
save(results, file = "simulated_v2_results_3.RData")
    
