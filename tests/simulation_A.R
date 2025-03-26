

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



breast_g <- get(load("../GRACKLE_data/data/breast_igraph_prob_1_cor_0_05.RData"))

num_nodes = opt$num_nodes
num_groups = opt$num_groups
num_per_group = opt$num_per_group

### noise test
noise_sequence <- c(seq(.3,.7, .2))
t_sequence <- c(0,seq(.1,.3,.1))

num_nodes = 400
num_groups = 5
num_per_group = 20


rand_graph_file = paste0("../GRACKLE_results/small_simulation_A/random_graph_", num_nodes, "nodes.RData")


if(file.exists(rand_graph_file)) {
    print("loading random graph")
    load(rand_graph_file)
} else {
    print('generating random graph')
                                       # generate random breast_g# generate random nework 
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


    noisy_expression_file <-  paste0("../GRACKLE_results/small_simulation_A/data/noisy_simulations_num_nodes_", num_nodes,"_num_groups_", num_groups,"_per_group_", num_per_group,"_t_",x,".RData")

    
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

            iterations <- opt$iterations
            
            for(d in 1:iterations) {

                cat("\n")
                print(paste0("running simulation for transitivity: ", network_t, " and modularity: ", network_noise, " at expression noise:", noise_sequence[ns]))
                print(paste("iteration", d, "out of", iterations))

                
                dat <- split_data(noisy_expression[[ns]], metadata , training_size = .7, seed = d)
##                print(head(rownames(dat$train_expression)))

                ## print(head(rownames(dat$train_expression)))

                
                pat_sim <- as.matrix(dat$train_metadata) %*% t(dat$train_metadata)

                
                min_vals <- apply(dat$train_expression,2, min)
                max_vals <- apply(dat$train_expression,2, max)

                ## min-max scale the input matrix
                dat$train_expression <-as.matrix( min_max_scale(dat$train_expression,min_vals,max_vals))
                dat$test_expression <- as.matrix(min_max_scale(dat$test_expression,min_vals,max_vals))

                i_seq <-seq(0,1,.1)
                j_seq <-seq(0,1,.1)

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
##                    print(round(i/nrow(grid_search),3))

                    print(i)
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
                        iterations = 300)

                            z = 50
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


                    ##                    scores[[i]] <- score
                    print(score)
                })

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





