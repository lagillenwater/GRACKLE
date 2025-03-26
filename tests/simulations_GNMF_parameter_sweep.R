
suppressMessages({
   suppressWarnings({
       library(tidyverse)
        library(igraph)
        library(devtools)
       library(parallel)
       library(optparse)
        load_all()
        library(reticulate)
        use_virtualenv("/mnt/grackle_env")
        library(tensorflow)
  })
})



breast_g <- get(load("../GRACKLE_data/data/breast_igraph_prob_1_cor_0_05.RData"))
### noise test
noise_sequence <- c(seq(.3,.7, .2))
t_sequence <- c(0,seq(.1,.3,.1))
num_nodes = 400
num_groups = 5
num_per_group = 20

rand_graph_file = paste0("../GRACKLE_results/small_simulation_A/random_graph_", num_nodes, "nodes_", Sys.Date(),".RData")


if(file.exists(rand_graph_file)) {
    print("loading random graph")
    load(rand_graph_file)
} else {
    print('generating random graph')
    g <- old_randomNetwork(prior_graph = breast_g, num_nodes = num_nodes)
    save(g,file = rand_graph_file)
}




for(x in t_sequence) {

    cat("\n")
    print(paste("t_sequence", x))

    if(x == 0) {
        new_g = g
    } else {
                new_g <- adjustTransitivity(g,target_transitivity = x)
    }
    print(new_g)
    network_t = transitivity(as_undirected(new_g), weights = NA)
    print(paste("network transivity: ", network_t))
        ## generate list of permuted networks with distinctly permuted modules
    set.seed(42)
    clusters <- cluster_louvain(as.undirected(new_g),weights = NA)
    membership <- clusters$membership
    membership_counts <- table(as.factor(membership))
    print(membership_counts)
    large_clusters <- names(membership_counts[order(membership_counts, decreasing = T)][1:5])
    cat("\n")
    network_noise <- round(modularity(clusters), 4)
    print(paste("prior network modularity: ", network_noise))

 

    noisy_expression_file <-  paste0("../GRACKLE_results/small_simulation_A/data/noisy_simulations_num_nodes_", num_nodes,"_num_groups_", num_groups,"_per_group_", num_per_group,"_t_",x,"_", Sys.Date(),".RData")

    
    if(file.exists(noisy_expression_file)) {
        print("loading noisy expression")
        load(noisy_expression_file)
    } else {
        set.seed(42)
        print("creating noisy expression")
        noisy_expression <- lapply(noise_sequence, function(z){
            ## simulate gene expression for the permuted networks
            sim_exp <- mclapply(large_clusters, function(x) {parallelSimulateExpression(g, iterations = 3, max_expression = 100,num_samples = num_per_group,select_nodes = which(membership == x),perturbation_constant = 2.8,noise_constant = z)}, mc.cores = 8)
                    ## combine expression data
            sim_exp <- lapply(sim_exp, as.data.frame)
            combined_exp <- do.call(rbind, sim_exp)
            rownames(combined_exp) <- paste0("sample", 1:nrow(combined_exp))
            
                    return(combined_exp)
        })

        names(noisy_expression) <- noise_sequence
        print(dim(noisy_expression[[1]]))
        print(nrow(noisy_expression[[1]]))
        metadata <- simulateMetadata(group_size = rep(num_per_group,num_groups), group_labels = paste0("subgroup", 1:num_groups), row_names = paste0("sample", 1:nrow(noisy_expression[[1]])))
        print("saving noisy expression")
        save(noisy_expression, metadata, file = noisy_expression_file)

    }
    
    mod =     round(modularity(cluster_louvain(as.undirected(new_g), weights = NA)), 1)

    nn_sequence <- seq(mod,mod - .3,-.1)

    for(y in nn_sequence) {
            if(y !=0) {
            print("adding network noise")
            new_g <- networkNoise(g,goal_modularity = y, membership = membership,large_clusters = large_clusters)
        } 

        cat("\n")
        network_noise <- round(modularity(cluster_louvain(as.undirected(new_g), weights = NA)), 4)
        print(paste("network modularity: ", network_noise))
        
        ##split the data into train and test
        g_adjacency <- as_adjacency_matrix(new_g)

        results <- list()
        
        for(ns in 1:length(noise_sequence)) {

            cat("\n")
            print(paste("noise_sequence", noise_sequence[ns]))
            res <- list()

            iterations <- 100
            
            for(d in 1:iterations) {

                cat("\n")
                print(paste0("running simulation for transitivity: ", network_t, " and modularity: ", network_noise, " at expression noise:", noise_sequence[ns]))
                print(paste("iteration", d, "out of", iterations))

                
                dat <- split_data(noisy_expression[[ns]], metadata , training_size = .7, seed = d)
##                print(head(rownames(dat$train_expression)))

                ## print(head(rownames(dat$train_expression)))

                
                min_vals <- apply(dat$train_expression,2, min)
                max_vals <- apply(dat$train_expression,2, max)

                ## min-max scale the input matrix
                dat$train_expression <-as.matrix( min_max_scale(dat$train_expression,min_vals,max_vals))
                dat$test_expression <- as.matrix(min_max_scale(dat$test_expression,min_vals,max_vals))


                grnmf_res <- data.frame(alpha = c(10^ 0, 10^1, 10^2, 10^3, 10^4), score = 0)
                ##grnmf_res <- data.frame(alpha = seq(.1), score = 0)
                for (a in grnmf_res$alpha) { 
                
                    grnmf <- runGRNMF(expression_matrix = dat$train_expression, k = num_groups, seed = 42, max_ = 300, alpha = a)
                    tmp_res <- evaluationWrapper(test_expression = dat$test_expression,
                                                   test_metadata = dat$test_metadata,
                                                   H_train = grnmf$H,
                                                   k = num_groups,
                                                   clusters = clusters,
                                                   aligned_clusters = large_clusters)
                    print(paste("GRNMF score for alpha", a, " = ", tmp_res))
                    grnmf_res[grnmf_res$alpha == a, "score"] <- tmp_res
                    } 

                


                res[[d]] <- list(grnmf_res = grnmf_res)
                ## res[[d]] <- list(grid_search = grid_search, transitivity = network_t, modularity = network_noise, expression_noise = noise_sequence[ns])



                cat("\n")
                
            }

            results[[ns]] <- res
            
        noise_simulation_results_file <- paste0("../GRACKLE_results/small_simulation_A/simulation_grnmf_sweep_results_num_nodes_", num_nodes,"_num_groups_", num_groups,"_num_samples_", num_per_group,"_t_", x, "_m_", y, "_2025_3_25.RData")
        print(paste("saving", noise_simulation_results_file))
        save(results, file = noise_simulation_results_file)

        cat("\n")

        }
            


        
    }
}





