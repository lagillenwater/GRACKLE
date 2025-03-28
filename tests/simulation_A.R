# Load required libraries
library(tidyverse)
library(igraph)
library(devtools)
library(parallel)
library(optparse)
load_all()
library(reticulate)
use_virtualenv("/mnt/grackle_env")
library(tensorflow)

# Load breast graph data
breast_g <- get(load("../GRACKLE_data/data/breast_igraph_prob_1_cor_0_05.RData"))

# Define parameters
noise_sequence <- c(seq(0.3, 0.7, 0.2))
t_sequence <- c(0, seq(0.1, 0.3, 0.1))
num_nodes <- 400
num_groups <- 5
num_per_group <- 20

# File path for random graph
rand_graph_file <- paste0("../GRACKLE_results/small_simulation_A/random_graph_", num_nodes, "nodes.RData")

# Load or generate random graph
if (file.exists(rand_graph_file)) {
    print("Loading random graph")
    load(rand_graph_file)
} else {
    print("Generating random graph")
    g <- old_randomNetwork(prior_graph = breast_g, num_nodes = num_nodes)
    save(g, file = rand_graph_file)
}

# Iterate over transitivity sequence
for (x in t_sequence) {
    cat("\n")
    print(paste("t_sequence:", x))

    # Adjust transitivity
    if (x == 0) {
        new_g <- g
    } else {
        new_g <- adjustTransitivity(g, target_transitivity = x)
    }

    # Calculate network transitivity
    network_t <- transitivity(as.undirected(new_g), weights = NA)
    print(paste("Network transitivity:", network_t))

    # Cluster the network
    set.seed(42)
    clusters <- cluster_louvain(as.undirected(new_g), weights = NA)
    membership <- clusters$membership
    membership_counts <- table(as.factor(membership))
    print(membership_counts)

    # Identify large clusters
    large_clusters <- names(membership_counts[order(membership_counts, decreasing = TRUE)][1:5])

    # Calculate network modularity
    network_noise <- round(modularity(clusters), 4)
    print(paste("Prior network modularity:", network_noise))

    # File path for noisy expression data
    noisy_expression_file <- paste0(
        "../GRACKLE_results/small_simulation_A/data/noisy_simulations_num_nodes_",
        num_nodes, "_num_groups_", num_groups, "_per_group_", num_per_group, "_t_", x, ".RData"
    )

    # Load or generate noisy expression data
    if (file.exists(noisy_expression_file)) {
        print("Loading noisy expression")
        load(noisy_expression_file)
    } else {
        set.seed(42)
        print("Creating noisy expression")
        noisy_expression <- lapply(noise_sequence, function(z) {
            sim_exp <- mclapply(large_clusters, function(cluster) {
                parallelSimulateExpression(
                    g, iterations = 3, max_expression = 100,
                    num_samples = num_per_group, select_nodes = which(membership == cluster),
                    perturbation_constant = 2.8, noise_constant = z
                )
            }, mc.cores = 8)

            # Combine expression data
            sim_exp <- lapply(sim_exp, as.data.frame)
            combined_exp <- do.call(rbind, sim_exp)
            rownames(combined_exp) <- paste0("sample", 1:nrow(combined_exp))
            return(combined_exp)
        })

        names(noisy_expression) <- noise_sequence
        metadata <- simulateMetadata(
            group_size = rep(num_per_group, num_groups),
            group_labels = paste0("subgroup", 1:num_groups),
            row_names = paste0("sample", 1:nrow(noisy_expression[[1]]))
        )

        print("Saving noisy expression")
        save(noisy_expression, metadata, file = noisy_expression_file)
    }

    # Modularity sequence
    mod <- round(modularity(cluster_louvain(as.undirected(new_g), weights = NA)), 1)
    nn_sequence <- seq(mod, mod - 0.3, -0.1)

    # Iterate over modularity sequence
    for (y in nn_sequence) {
        if (y != 0) {
            print("Adding network noise")
            new_g <- networkNoise(g, goal_modularity = y, membership = membership, large_clusters = large_clusters)
        }

        # Calculate network modularity
        network_noise <- round(modularity(cluster_louvain(as.undirected(new_g), weights = NA)), 4)
        print(paste("Network modularity:", network_noise))

        # Split data into train and test
        g_adjacency <- as_adjacency_matrix(new_g)
        results <- list()

        # Iterate over noise sequence
        for (ns in seq_along(noise_sequence)) {
            cat("\n")
            print(paste("noise_sequence:", noise_sequence[ns]))
            res <- list()

            # Perform iterations
            for (d in seq_len(opt$iterations)) {
                cat("\n")
                print(paste0(
                    "Running simulation for transitivity: ", network_t,
                    " and modularity: ", network_noise,
                    " at expression noise: ", noise_sequence[ns]
                ))
                print(paste("Iteration", d, "out of", opt$iterations))

                # Split data
                dat <- split_data(noisy_expression[[ns]], metadata, training_size = 0.7, seed = d)

                # Min-max scale the input matrix
                min_vals <- apply(dat$train_expression, 2, min)
                max_vals <- apply(dat$train_expression, 2, max)
                dat$train_expression <- as.matrix(min_max_scale(dat$train_expression, min_vals, max_vals))
                dat$test_expression <- as.matrix(min_max_scale(dat$test_expression, min_vals, max_vals))

                # Perform GRACKLE grid search
                i_seq <- seq(0, 1, 0.1)
                j_seq <- seq(0, 1, 0.1)
                grid_search <- expand.grid(lambda_1 = i_seq, lambda_2 = j_seq)
                grid_search$score <- 0

                print("GRACKLE grid search ...")
                scores <- lapply(seq_len(nrow(grid_search)), function(i) {
                    grackle <- GRACKLE(
                        Y = dat$train_expression,
                        net_similarity = as.matrix(g_adjacency),
                        patient_similarity = as.matrix(dat$train_metadata) %*% t(dat$train_metadata),
                        diff_threshold = 1e-4,
                        lambda_1 = grid_search$lambda_1[i],
                        lambda_2 = grid_search$lambda_2[i],
                        k = num_groups,
                        verbose = FALSE,
                        error_terms = FALSE,
                        beta = 0,
                        iterations = 100
                    )

                    # Evaluate GRACKLE results
                    evaluationWrapper(
                        test_expression = dat$test_expression,
                        test_metadata = dat$test_metadata,
                        H_train = grackle$H,
                        k = num_groups,
                        clusters = clusters,
                        aligned_clusters = large_clusters
                    )
                })

                grid_search$score <- unlist(scores)
                print(paste("Avg GRACKLE score:", round(mean(grid_search$score), 3)))
                print(paste("Top GRACKLE score:", max(grid_search$score)))

                # Save results
                res[[d]] <- list(
                    grid_search = grid_search,
                    transitivity = network_t,
                    modularity = network_noise,
                    expression_noise = noise_sequence[ns]
                )
            }

            results[[ns]] <- res
        }

        # Save simulation results
        noise_simulation_results_file <- paste0(
            "../GRACKLE_results/small_simulation_A/updated_simulation_results_num_nodes_",
            num_nodes, "_num_groups_", num_groups, "_num_samples_", num_per_group,
            "_t_", x, "_m_", y, ".RData"
        )
        print(paste("Saving", noise_simulation_results_file))
        save(results, file = noise_simulation_results_file)
    }
}