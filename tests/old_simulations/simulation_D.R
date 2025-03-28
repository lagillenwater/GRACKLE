## generate list of permuted networks with distinctly permuted modules
set.seed(42)
clusters <- cluster_louvain(as.undirected(g),weights = NA, resolution =1)
membership <- clusters$membership
membership_counts <- table(as.factor(membership))
print(membership_counts)
large_clusters <- sample(names(membership_counts),num_groups)
## large_clusters <- names(membership_counts[order(membership_counts, decreasing = F)][1:5])

cat("\n")
network_noise <- round(modularity(clusters), 4)
print(paste("prior network modularity: ", network_noise))

noisy_expression_file <-  paste0("../GRACKLE_results/small_simulation_D/data/noisy_simulations_num_nodes_", num_nodes,"_num_groups_", num_groups,"_per_group_", num_per_group,".RData")

noise_sequence = .5

if(file.exists(noisy_expression_file)) {
    print("loading noisy expression")
    load(noisy_expression_file)
} else {
    set.seed(42)
    print("creating noisy expression")
    ## simulate gene expression for the permuted networks
    sim_exp <- lapply(large_clusters, function(x) {
        group_exp <- mclapply(1:num_per_group, function(y) {
            simulateExpression(g, iterations = 2, max_expression = 100,select_nodes = which(membership == x),perturbation_constant = 2.8,noise_constant =noise_sequence )
        }, mc.cores = 20)
        tmp <- t(do.call(cbind,group_exp))
        return(tmp)
    })
    ## combine expression data
    sim_exp <- lapply(sim_exp, as.data.frame)
    combined_exp <- do.call(rbind, sim_exp)
    rownames(combined_exp) <- paste0("sample", 1:nrow(combined_exp))
    noisy_expression <- combined_exp
    metadata <- simulateMetadata(group_size = rep(num_per_group,num_groups), group_labels = paste0("subgroup", 1:num_groups), row_names = paste0("sample", 1:nrow(noisy_expression)))
    print("saving noisy expression")
    save(noisy_expression, metadata, file = noisy_expression_file)
}


mod =     round(modularity(cluster_louvain(as.undirected(g), weights = NA)), 1)



##split the data into train and test
g_adjacency <- as_adjacency_matrix(g)

results <- list()
        
iterations <- opt$iterations


iterations <- 100            
for(d in 1:iterations) {
    cat("\n")
    print(paste("iteration", d, "out of", iterations))
    dat <- split_data(noisy_expression, metadata , training_size = .7, seed = d)
    ##                print(head(rownames(dat$train_expression)))
    pat_sim <- as.matrix(dat$train_metadata) %*% t(dat$train_metadata)
    min_vals <- apply(dat$train_expression,2, min)
    max_vals <- apply(dat$train_expression,2, max)
    ## min-max scale the input matrix
    dat$train_expression <-as.matrix( min_max_scale(dat$train_expression,min_vals,max_vals))
    dat$test_expression <- as.matrix(min_max_scale(dat$test_expression,min_vals,max_vals))
    i_seq <-c(0,.1,1,10,100)
    j_seq <-c(0,.1,1,10,100)
    
    print("GRACKLE grid search ...") 
    grid_search <- as.data.frame(expand.grid(i_seq,j_seq))
    names(grid_search) <- c("lambda_1", "lambda_2")
    grid_search$score <- 0

    
    ## cl <- makeCluster(detectCores() -1)
    ## clusterEvalQ(cl, .libPaths("/home/lucas/R/x86_64-pc-linux-gnu-library/4.1"))
    ## clusterEvalQ(cl,{

    ##     library(devtools)
    ##     load_all()
    ## })
    ## clusterExport(cl, varlist = c("pat_sim", "g_adjacency", "dat", "num_groups", "grid_search", "large_clusters", "clusters"))       
    
    scores <- lapply(1:nrow(grid_search), function(i) {
        ## scores <- list()
        ## for( i in 1: nrow(grid_search)) {
        ##     ## run GRACKLE NMF
        ##     print(round(i/nrow(grid_search),3))
        grackle <- GRACKLE(
            Y = dat$train_expression,
            net_similarity = as.matrix(g_adjacency),
            patient_similarity = pat_sim,
            diff_threshold = 1e-5,
            lambda_1 = grid_search$lambda_1[i],
            lambda_2 =  grid_search$lambda_2[i],
            k = num_groups,
            verbose = F,
            error_terms = F,
            beta = 0,
            iterations = 50)
        z = 50
        ## correspondence between selected W LV's and top loading gene modules
        while(any(is.na(grackle$H))) {
            z = z-5
        grackle <- GRACKLE(
            Y = dat$train_expression,
            net_similarity = as.matrix(g),
            patient_similarity = pat_sim,
            diff_threshold = 1e-4,
            lambda_1 = grid_search$lambda_1[i],
            lambda_2 =grid_search$lambda_2[i],
            k = k,
            verbose = F,
            beta = 0,
            error_terms = F,
            iterations = z)
        }

        ## correspondence between selected W LV's and top loading gene modules
        score <- evaluationWrapper(test_expression = dat$test_expression,
                                   test_metadata = dat$test_metadata,
                                   H_train = grackle$H,
                                   k = num_groups,
                                   clusters = clusters,
                                   aligned_clusters = large_clusters)


        ##                    scores[[i]] <- score
    } )

    ##                stopCluster(cl)

    grid_search$score <- unlist(scores)

    

    print(paste("avg GRACKLE score", round(mean(grid_search$score),3)))
    print(paste("top GRACKLE score", max(grid_search %>% filter(scores == max(grid_search$score)) %>% .$score)))
    

    ## nmf <- runNMF(dat$train_expression,num_groups, "lee", seed = 42)
    ## nmf_res <- evaluationWrapper(test_expression = dat$test_expression,
    ##                              test_metadata = dat$test_metadata,
    ##                              H_train = nmf$H,
    ##                              k = num_groups,
    ##                              clusters = clusters,
    ##                              aligned_clusters = large_clusters)
    nmf_res <- grid_search$score[1]
    print(paste("NMF score", nmf_res))

    
    grnmf <- runGRNMF(expression_matrix = dat$train_expression, k = num_groups, seed = 42, max_ = 300, alpha = 1)
    grnmf_res <- evaluationWrapper(test_expression = dat$test_expression,
                                   test_metadata = dat$test_metadata,
                                   H_train = grnmf$H,
                                   k = num_groups,
                                   clusters = clusters,
                                   aligned_clusters = large_clusters)

    print(paste("GRNMF score", grnmf_res))


    res[[d]] <- list(grid_search = grid_search,nmf_res = nmf_res, grnmf_res = grnmf_res, expression_noise = noise_sequence)
    ## res[[d]] <- list(grid_search = grid_search, transitivity = network_t, modularity = network_noise, expression_noise = noise_sequence[ns])



    cat("\n")
    
}





noise_simulation_results_file <- paste0("../GRACKLE_results/small_simulation_D/simulation_results_num_nodes_", num_nodes,"_num_groups_", num_groups,"_num_samples_", num_per_group,"_noise_", noise_sequence, ".RData")
print(paste("saving", noise_simulation_results_file))
save(results, file = noise_simulation_results_file)
cat("\n")



