## set seed
## set.seed(42)
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
        load_all()
    })
})

breast_g <- get(load("../GRACKLE_data/data/breast_igraph_prob_1_cor_0_05.RData"))

degrees <- degree(breast_g)

## pdf("./results/networks/degree_scatter.pdf")
## qplot(y  = seq_along(degrees), x= degrees, geom = "point", ylab = "" ) + theme_classic()
## dev.off()

num_nodes = 500
### noise test
noise_sequence <- c(seq(.3,.7,.1))

t_sequence <- c(0,seq(.1,.3,.1))
nn_sequence <- c(0,seq(.8,.6,-.1))

#noise_sequence <- .3
## t_sequence <- 0
## nn_sequence <- 0

num_groups = 5
num_per_group = 20

rand_graph_file = paste0("../GRACKLE_results/small_simulation_A/random_graph_", num_nodes, "nodes.RData")

if(file.exists(rand_graph_file)) {
    print("loading random graph")
    load(rand_graph_file)
} else {

                                        # generate random breast_g# generate random nework 
    g <- randomNetwork(prior_graph = breast_g, num_nodes = num_nodes)
    save(g,file = rand_graph_file)
}

## g_title <- paste0("T: ", round(transitivity(g,type = "global"),4))

## fixed_layout= layout_with_lgl(g)
## plot network
# pdf("./results/networks/network.pdf")
## plot(
##   g, 
##   vertex.size = 3,  # Adjust node size as needed
##   vertex.label = NA,  # Hide vertex labels for clarity
##   edge.width = 1,  # Adjust edge width as needed
##   edge.arrow.size = .15,
##   layout = fixed_layout,
##   main = g_title
  
## )
# dev.off()


## for( x in seq(.1,.3,.1)) {

##     new_g <- adjustTransitivity(g,target_transitivity = x, max_iter = 1e4)
  
##   g_title <- paste0("T: ", round(transitivity(new_g,type = "global"),4))
  
##   print(g_title)
    
##   plot(
##     new_g, 
##     vertex.size = 3,  # Adjust node size as needed
##     vertex.label = NA,  # Hide vertex labels for clarity
##     edge.width = 1,  # Adjust edge width as needed
##     edge.arrow.size = .15,
##     layout = fixed_layout,
##     main = g_title
##   )
## }
  
##   # pdf("./results/networks/network.pdf")
## plot(
##   new_g, 
##   vertex.size = 3,  # Adjust node size as needed
##   vertex.label = NA,  # Hide vertex labels for clarity
##   edge.width = 1,  # Adjust edge width as needed
##   edge.arrow.size = .15,
##   layout = layout_with_lgl
## )
# dev.off()


## permute modules

# Assign colors to each community
## V(g)$color <- membership(clusters)

## clusters_to_highlight <- lapply(large_clusters, function(x) clusters[[x]])
# Plot the graph with module overlay
#pdf("./results/networks/network_odules.pdf")

## g_title = paste("M:", round(modularity(clusters), 4))

## plot(
##   clusters, 
##   g, 
##   mark.groups = clusters_to_highlight,
##   vertex.size = 3,  # Adjust node size as needed
##   vertex.label = NA,  # Hide vertex labels for clarity
##   edge.width = 0.3,  # Adjust edge width as needed
##   edge.arrow.size = .15,
##   layout = fixed_layout, 
##   main = g_title
## )
#dev.off()

## for(x in seq(.7,.5,-.1,)) {
##   new_g <- networkNoise(g,goal_modularity = x, membership = membership,large_clusters = large_clusters)
  
##   g_title <- paste("M:", round(modularity(cluster_louvain(as.undirected(new_g), weights = NA)), 4))
  
##   print(g_title)
  
##   plot(
##     clusters, 
##     new_g, 
##     mark.groups = clusters_to_highlight,
##     vertex.size = 3,  # Adjust node size as needed
##     vertex.label = NA,  # Hide vertex labels for clarity
##     edge.width = 0.3,  # Adjust edge width as needed
##     edge.arrow.size = .15,
##     layout = fixed_layout, 
##     main = g_title
##   )
 
## }

## plot a single cluster
## cluster_id <- large_clusters[1]

## vertices_in_cluster <- which(membership(clusters) == cluster_id)

## # Create a subgraph with only the vertices in the selected cluster
## subgraph <- induced_subgraph(g, vids = vertices_in_cluster)


## # Plot the subgraph
## pdf("./results/networks/single_module.pdf")
## plot(
##   subgraph, 
##   vertex.size = 10,  # Adjust node size as needed
##   vertex.label = NA,  # Hide vertex labels for clarity
##   edge.width = 0.5,
##   edge.arrow.size = .5,# Adjust edge width as needed,
##   layout = layout_with_lgl
## )
## dev.off()




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

    clusters <- cluster_louvain(as.undirected(new_g),weights = NA)
    membership <- clusters$membership
    membership_counts <- table(as.factor(membership))
    print(membership_counts)
    large_clusters <- names(membership_counts[order(membership_counts, decreasing = T)][1:num_groups])

    noisy_expression_file <-  paste0("../GRACKLE_results/small_simulation_A/data/noisy_simulations_num_nodes_", num_nodes,"_num_groups_", num_groups,"_per_group_", num_per_group,"_t_",x,".RData")


    if(!file.exists(noisy_expression_file)) {

        print("creating noisy expression")
        noisy_expression <- mclapply(noise_sequence, function(z){

            ## simulate gene expression for the permuted networks
            sim_exp <- mclapply(large_clusters, function(x) parallelSimulateExpression(g, iterations = 3, max_expression = 100,num_samples = num_per_group,select_nodes = which(membership == x),perturbation_constant = 2.8,noise_constant = z), mc.cores = length(large_clusters))
            
            ## combine expression data
            sim_exp <- lapply(sim_exp, as.data.frame)
            combined_exp <- do.call(rbind, sim_exp)
            rownames(combined_exp) <- paste0("sample", 1:nrow(combined_exp))
            
            
            return(combined_exp)
        }, mc.cores = length(noise_sequence))

        names(noisy_expression) <- noise_sequence
        print("saving noisy expression")
        save(noisy_expression, file = noisy_expression_file)
        
    } else {
        print("loading noisy expression")
        load(noisy_expression_file)
    }

    metadata <- simulateMetadata(group_size = rep(num_per_group,num_groups), group_labels = paste0("subgroup", 1:num_groups), row_names = paste0("sample", 1:nrow(noisy_expression[[1]])))

        ## lapply(noise_sequence, function(x){
        ##   pdf(paste0("./results/simulations/plots/permuted_network_heatmap_",x,"_noise.pdf"))
        ##   print(plotSimulatedExpression(noisy_expression[[as.character(x)]],metadata))
        ##   dev.off()
        ## })

    
    for(y in nn_sequence) {


        if(y !=0) {
            new_g <- networkNoise(g,goal_modularity = y, membership = membership,large_clusters = large_clusters)
        } 

        cat("\n")
        print(paste0("runing simulation for transitivity: ", x, " and modularity: ", y))

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
                
                print(paste("iteration", d, "out of", iterations))
                dat <- split_data(noisy_expression[[ns]], metadata , training_size = .7)

                pat_sim <- as.matrix(dat$train_metadata) %*% t(dat$train_metadata)

                min_vals <- apply(dat$train_expression,2, min)
                max_vals <- apply(dat$train_expression,2, max)

                ## min-max scale the input matrix
                dat$train_expression <-as.matrix( min_max_scale(dat$train_expression,min_vals,max_vals))
                dat$test_expression <- as.matrix(min_max_scale(dat$test_expression,min_vals,max_vals))



                i_seq <-seq(0,1,.1)
                j_seq <-seq(0,1,.1)

                
                grid_search <- as.data.frame(expand.grid(i_seq,j_seq))
                names(grid_search) <- c("lambda_1", "lambda_2")
                grid_search$score <- 0

                cl <- makeCluster(detectCores() -1)
                clusterEvalQ(cl, .libPaths("/home/lucas/R/x86_64-pc-linux-gnu-library/4.1"))
                clusterEvalQ(cl,{

                    library(devtools)
                    load_all()
                })
                clusterExport(cl, varlist = c("pat_sim", "g_adjacency", "dat", "num_groups", "grid_search", "large_clusters", "clusters"))       
                
                scores <- parLapply(cl,1:nrow(grid_search), function(i) {

                    ## run GRACKLE NMF

                    grackle <- GRACKLE(
                        Y = dat$train_expression,
                        net_similarity = as.matrix(g_adjacency),
                        patient_similarity = pat_sim,
                        diff_threshold = 1e-6,
                        lambda_1 = grid_search$lambda_1[i],
                        lambda_2 =  grid_search$lambda_2[i],
                        k = num_groups,
                        verbose = F,
                        beta = .0)


                    ## correspondence between selected W LV's and top loading gene modules
                    score <- evaluationWrapper(test_expression = dat$test_expression,
                                               test_metadata = dat$test_metadata,
                                               H_train = grackle$H,
                                               k = num_groups,
                                               clusters = clusters,
                                               aligned_clusters = large_clusters)
                    return(score)
                })

                stopCluster(cl)

                grid_search$score <- unlist(scores)

                print(paste("avg GRACKLE score", round(mean(grid_search$score),3)))
                print(paste("top GRACKLE score", max(grid_search %>% filter(scores == max(grid_search$score)) %>% .$score)))
                

                nmf <- runNMF(dat$train_expression,num_groups, "lee", seed = 42)
                nmf_res <- evaluationWrapper(test_expression = dat$test_expression,
                                             test_metadata = dat$test_metadata,
                                             H_train = nmf$H,
                                             k = num_groups,
                                             clusters = clusters,
                                             aligned_clusters = large_clusters)

                print(paste("NMF score", nmf_res))

                
                grnmf <- runGRNMF(expression_matrix = dat$train_expression, k = num_groups, seed = 42, max_ = 300, alpha = .1)
                grnmf_res <- evaluationWrapper(test_expression = dat$test_expression,
                                               test_metadata = dat$test_metadata,
                                               H_train = grnmf$H,
                                               k = num_groups,
                                               clusters = clusters,
                                               aligned_clusters = large_clusters)

                print(paste("GRNMF score", grnmf_res))


                res[[d]] <- list(grid_search = grid_search,nmf_res = nmf_res, grnmf_res = grnmf_res, transitivity = network_t, modularity = network_noise, expression_noise = noise_sequence[ns])

                cat("\n")
                
            }

            results[[ns]] <- res
        }
            


        noise_simulation_results_file <- paste0("../GRACKLE_results/small_simulation_A/simulation_results_num_nodes_", num_nodes,"_num_groups_", num_groups,"_num_samples_", num_per_group,"_t_", x, "_m_", y, ".RData")
        print(paste("saving", noise_simulation_results_file))
        save(results, file = noise_simulation_results_file)
        cat("\n")

        
    }
}






        
