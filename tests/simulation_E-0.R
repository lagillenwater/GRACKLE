 ## set seed
                                        #
## setwd("~/OneDrive - The University of Colorado Denver/Projects/GRACKLE/")
suppressMessages({
    suppressWarnings({
#        library(tidyverse)
        library(igraph)
        library(devtools)
  #      library(ggplot2)
 #       library(ComplexHeatmap)
        library(parallel)
        library(optparse)
        load_all()
        library(reticulate)
        use_virtualenv("/mnt/grackle_env")
        library(tensorflow)
        library(aricode)
        library(dynutils)
    })
})




                                        # Define the list of options
option_list <- list(
    make_option(c("-n", "--num_nodes"), type = "numeric",  help = "Number of nodes int he graph", metavar = "NUMBER"),
    make_option(c("-i", "--iterations"), type = "numeric", help = "Iterations to run", metavar = "NUMBER"),
    make_option(c("-g", "--num_per_group"), type = "numeric",  help = "Number of samples per group", metavar = "NUMBER"),
    make_option(c("-k", "--num_groups"), type = "numeric",  help = "Number of groups", metavar = "NUMBER"),
    make_option(c("-c", "--num_cores"), type = "numeric",  help = "Number of cores", metavar = "NUMBER")
    
    
)


## Parse the command line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

load( file = "../GRACKLE_data/data/Breast/TCGA/directed_breast_g_with_PAM50.RData")

num_nodes = opt$num_nodes
num_groups = opt$num_groups
num_per_group = opt$num_per_group

### noise test
noise_sequence <- c(seq(.9,.3, -.3))
t_sequence <- 0

## num_nodes = 1000
## num_groups = 5
## num_per_group = 20

rand_graph_file = paste0("../GRACKLE_results/small_simulation_A/random_graph_", num_nodes, "nodes.RData")


if(file.exists(rand_graph_file)) {
    print("loading random graph")
    load(rand_graph_file)
} else {
    print('generating random graph')
                                       # generate random breast_g# generate random nework 
    g <- randomNetwork(prior_graph = directed_breast_g_with_PAM50, num_nodes = num_nodes)
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
    
    network_t = transitivity(as_undirected(new_g), weights = NA)
    print(paste("network transivity: ", network_t))
    
    ## generate list of permuted networks with distinctly permuted modules
    set.seed(42)
    clusters <- cluster_louvain(as.undirected(new_g),weights = NA, resolution = 10)
    membership <- clusters$membership
    membership_counts <- table(as.factor(membership))
    large_clusters <- names(membership_counts[order(membership_counts, decreasing = T)][1:5])
    print(membership_counts[large_clusters])
    
    noisy_expression_file <-  paste0("../GRACKLE_results/small_simulation_A/data/noisy_simulations_num_nodes_", num_nodes,"_num_groups_", num_groups,"_per_group_", num_per_group,"_t_",x,".RData")

    
    if(file.exists(noisy_expression_file)) {
        print("loading noisy expression")
        load(noisy_expression_file)
    } else {

                ##split the data into train and test
        g_adjacency <- as_adjacency_matrix(new_g)
        set.seed(42)
        print("creating noisy expression")
        metadata <- simulateMetadata(group_size = rep(num_per_group,num_groups), group_labels = paste0("subgroup", 1:num_groups), row_names = paste0("sample", 1:(num_per_group*num_groups)))
        noisy_expression <- lapply(noise_sequence, function(z){
            ## simulate gene expression for the permuted networks
            n_row <- num_groups*num_per_group
            sim_exp <- matrix(nrow =n_row , ncol = num_nodes, runif(n_row*num_nodes, min = 0, max = 20))
            rownames(sim_exp) <- paste0("sample", 1:n_row)
            ## simulate up/down regulation of pathways
            for(i in 1:num_groups) {
##                print(i)
                end <- i * num_per_group
                start <- end - num_per_group + 1
                tochange <- sample(start:end, ceiling(length(start:end) * z))
                nodes <- clusters[large_clusters[i]][[1]]
                sim_exp[tochange,nodes] <- runif(tochange*length(nodes),min = 20, max = 25)
            }
            dist_matrix <- dist(sim_exp)
            hc <- hclust(dist_matrix)
            tmp <- cbind(metadata %>% rownames_to_column("sample") %>%pivot_longer(-sample) %>% filter(value ==1), sim = cutree(hc, k = num_groups))
            print(ARI(tmp$name, tmp$sim))
            return(sim_exp)
        })
        names(noisy_expression) <- noise_sequence
        print("saving noisy expression")
        save(noisy_expression, metadata, file = noisy_expression_file)

    }


    mod =     round(modularity(cluster_louvain(as.undirected(new_g), weights = NA)), 1)

    nn_sequence <- seq(.3,.1,-.1)

    for(y in nn_sequence) {
        
        print("adding network noise")
        new_g <- networkNoise(new_g,goal_modularity = y, membership = membership,large_clusters = large_clusters)
        
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

            iterations <- opt$iterations
            
            for(d in 1:iterations) {
               
                cat("\n")
                print(paste0("running simulation for transitivity: ", network_t, " and modularity: ", network_noise, " at expression noise:", noise_sequence[ns], " ----- iteration ", d, " out of ", iterations))
                dat <- split_data(noisy_expression[[ns]], metadata , training_size = .7, seed = d)
                pat_sim <- as.matrix(dat$train_metadata) %*% t(dat$train_metadata)
                min_vals <- apply(dat$train_expression,2, min)
                max_vals <- apply(dat$train_expression,2, max)
                ## min-max scale the input matrix
                dat$train_expression <-as.matrix( min_max_scale(dat$train_expression,min_vals,max_vals))
                dat$test_expression <- as.matrix(min_max_scale(dat$test_expression,min_vals,max_vals))
                i_seq <-seq(0,100,10)
                j_seq <-seq(0,20,2)
                print("GRACKLE grid search ...") 
                grid_search <- as.data.frame(expand.grid(i_seq,j_seq))
                names(grid_search) <- c("lambda_1", "lambda_2")
                grid_search$score <- 0
                scores <- lapply(1:nrow(grid_search), function(i) {
                ##     ## run GRACKLE NMF
##                    print(i)
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
                            z = 25
                    ## correspondence between selected W LV's and top loading gene modules
                    while(any(is.na(grackle$H)|any(is.infinite(grackle$H)))|(max(grackle$H) - min(grackle$H)) >1e21) {
                        z = z-5
                        grackle <- GRACKLE(
                            Y = dat$train_expression,
                             net_similarity = as.matrix(g_adjacency),
                            patient_similarity = pat_sim,
                            diff_threshold = 1e-4,
                            lambda_1 = grid_search$lambda_1[i],
                            lambda_2 =grid_search$lambda_2[i],
                            k = num_groups,
                            verbose = F,
                            beta = 0,
                            error_terms = F,
                            iterations = z)
                       print(paste("GRACKLE iterations:",z))
                    }
                    ## correspondence between selected W LV's and top loading gene modules
                    score <- evaluationWrapper(test_expression = dat$test_expression,
                                               test_metadata = dat$test_metadata,
                                               H_train = grackle$H,
                                               k = num_groups,
                                               clusters = clusters,
                                               aligned_clusters = large_clusters)
                })
                grid_search$score <- unlist(scores)
                print(paste("avg GRACKLE score", round(mean(grid_search$score),3)))
                print(paste("top GRACKLE score", max(grid_search %>% filter(scores == max(grid_search$score)) %>% .$score)))
                print(grid_search %>% filter(scores == max(grid_search$score)))
                nmf_res <- grid_search$score[1]
                print(paste("NMF score", nmf_res))
                grnmf <- runGRNMF(expression_matrix = dat$train_expression, k = num_groups, seed = 42, max_ = 100, alpha = 1)
                grnmf_res <- evaluationWrapper(test_expression = dat$test_expression,
                                               test_metadata = dat$test_metadata,
                                               H_train = grnmf$H,
                                               k = num_groups,
                                               clusters = clusters,
                                               aligned_clusters = large_clusters)
                print(paste("GRNMF score", grnmf_res))
                res[[d]] <- list(grid_search = grid_search,nmf_res = nmf_res, grnmf_res = grnmf_res, transitivity = network_t, modularity = network_noise, expression_noise = noise_sequence[ns])
                ## res[[d]] <- list(grid_search = grid_search, transitivity = network_t, modularity = network_noise, expression_noise = noise_sequence[ns])
                cat("\n")
            }

            results[[ns]] <- res
        }
            


        noise_simulation_results_file <- paste0("../GRACKLE_results/small_simulation_A/updated_simulation_results_num_nodes_", num_nodes,"_num_groups_", num_groups,"_num_samples_", num_per_group,"_t_", x, "_m_", y, ".RData")
        print(paste("saving", noise_simulation_results_file))
        save(results, file = noise_simulation_results_file)
        cat("\n")

        
    }
}





