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

# Load the preprocessed breast network
breast_graph <- get(load("../GRACKLE_data/data/breast_igraph_prob_1_cor_0_05.RData"))

# Define parameters
noise_sequence <- seq(0.3, 0.7, 0.2)
transitivity_sequence <- c(0, seq(0.1, 0.3, 0.1))
num_nodes <- 400
num_groups <- 5
num_per_group <- 20

# File path for random graph
random_graph_file <- paste0(
  "../GRACKLE_results/small_simulation_A/random_graph_", num_nodes, "nodes.RData"
)

# Load or generate random graph
if (file.exists(random_graph_file)) {
  message("Loading random graph")
  load(random_graph_file)
} else {
  message("Generating random graph")
  random_graph <- old_randomNetwork(prior_graph = breast_graph, num_nodes = num_nodes)
  save(random_graph, file = random_graph_file)
}

# Iterate over transitivity sequence
for (transitivity in transitivity_sequence) {
  cat("\n")
  message(paste("Transitivity sequence:", transitivity))

  # Adjust transitivity
  if (transitivity == 0) {
    adjusted_graph <- random_graph
  } else {
    adjusted_graph <- adjustTransitivity(random_graph, target_transitivity = transitivity)
  }

  # Calculate network transitivity
  network_transitivity <- transitivity(as.undirected(adjusted_graph), weights = NA)
  message(paste("Network transitivity:", network_transitivity))

  # Cluster the network
  set.seed(42)
  clusters <- cluster_louvain(as.undirected(adjusted_graph), weights = NA)
  membership <- clusters$membership
  membership_counts <- table(as.factor(membership))
  print(membership_counts)

  # Identify large clusters
  large_clusters <- names(membership_counts[order(membership_counts, decreasing = TRUE)][1:5])

  # Calculate network modularity
  network_modularity <- round(modularity(clusters), 4)
  message(paste("Prior network modularity:", network_modularity))

  # File path for noisy expression data
  noisy_expression_file <- paste0(
    "../GRACKLE_results/small_simulation_A/data/noisy_simulations_num_nodes_",
    num_nodes, "_num_groups_", num_groups, "_per_group_", num_per_group, "_t_", transitivity, ".RData"
  )

  # Load or generate noisy expression data
  if (file.exists(noisy_expression_file)) {
    message("Loading noisy expression")
    load(noisy_expression_file)
  } else {
    set.seed(42)
    message("Creating noisy expression")
    noisy_expression <- lapply(noise_sequence, function(noise) {
      simulated_expression <- mclapply(large_clusters, function(cluster) {
        parallelSimulateExpression(
          random_graph, iterations = 3, max_expression = 100,
          num_samples = num_per_group, select_nodes = which(membership == cluster),
          perturbation_constant = 2.8, noise_constant = noise
        )
      }, mc.cores = 8)

      # Combine expression data
      simulated_expression <- lapply(simulated_expression, as.data.frame)
      combined_expression <- do.call(rbind, simulated_expression)
      rownames(combined_expression) <- paste0("sample", seq_len(nrow(combined_expression)))
      combined_expression
    })

    names(noisy_expression) <- noise_sequence
    metadata <- simulateMetadata(
      group_size = rep(num_per_group, num_groups),
      group_labels = paste0("subgroup", seq_len(num_groups)),
      row_names = paste0("sample", seq_len(nrow(noisy_expression[[1]])))
    )

    message("Saving noisy expression")
    save(noisy_expression, metadata, file = noisy_expression_file)
  }

  # Modularity sequence
  modularity <- round(modularity(cluster_louvain(as.undirected(adjusted_graph), weights = NA)), 1)
  modularity_sequence <- seq(modularity, modularity - 0.3, -0.1)

  # Iterate over modularity sequence
  for (modularity_value in modularity_sequence) {
    if (modularity_value != 0) {
      message("Adding network noise")
      adjusted_graph <- networkNoise(
        random_graph, goal_modularity = modularity_value,
        membership = membership, large_clusters = large_clusters
      )
    }

    # Calculate network modularity
    network_modularity <- round(modularity(cluster_louvain(as.undirected(adjusted_graph), weights = NA)), 4)
    message(paste("Network modularity:", network_modularity))

    # Split data into train and test
    graph_adjacency <- as_adjacency_matrix(adjusted_graph)
    results <- list()

    # Iterate over noise sequence
    for (noise_index in seq_along(noise_sequence)) {
      cat("\n")
      message(paste("Noise sequence:", noise_sequence[noise_index]))
      res <- list()

      # Perform iterations
      for (iteration in seq_len(opt$iterations)) {
        cat("\n")
        message(paste0(
          "Running simulation for transitivity: ", network_transitivity,
          " and modularity: ", network_modularity,
          " at expression noise: ", noise_sequence[noise_index]
        ))
        message(paste("Iteration", iteration, "out of", opt$iterations))

        # Split data
        data <- split_data(
          noisy_expression[[noise_index]], metadata,
          training_size = 0.7, seed = iteration
        )

        # Min-max scale the input matrix
        min_vals <- apply(data$train_expression, 2, min)
        max_vals <- apply(data$train_expression, 2, max)
        data$train_expression <- as.matrix(min_max_scale(data$train_expression, min_vals, max_vals))
        data$test_expression <- as.matrix(min_max_scale(data$test_expression, min_vals, max_vals))

        # Perform GRACKLE grid search
        lambda_1_sequence <- seq(0, 1, 0.1)
        lambda_2_sequence <- seq(0, 1, 0.1)
        grid_search <- expand.grid(lambda_1 = lambda_1_sequence, lambda_2 = lambda_2_sequence)
        grid_search$score <- 0

        message("GRACKLE grid search ...")
        scores <- lapply(seq_len(nrow(grid_search)), function(i) {
          grackle <- GRACKLE(
            Y = data$train_expression,
            net_similarity = as.matrix(graph_adjacency),
            patient_similarity = as.matrix(data$train_metadata) %*% t(data$train_metadata),
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
            test_expression = data$test_expression,
            test_metadata = data$test_metadata,
            H_train = grackle$H,
            k = num_groups,
            clusters = clusters,
            aligned_clusters = large_clusters
          )
        })

        grid_search$score <- unlist(scores)
        message(paste("Avg GRACKLE score:", round(mean(grid_search$score), 3)))
        message(paste("Top GRACKLE score:", max(grid_search$score)))

        # Save results
        res[[iteration]] <- list(
          grid_search = grid_search,
          transitivity = network_transitivity,
          modularity = network_modularity,
          expression_noise = noise_sequence[noise_index]
        )
      }

      results[[noise_index]] <- res
    }

    # Save simulation results
    simulation_results_file <- paste0(
      "../GRACKLE_results/small_simulation_A/updated_simulation_results_num_nodes_",
      num_nodes, "_num_groups_", num_groups, "_num_samples_", num_per_group,
      "_t_", transitivity, "_m_", modularity_value, ".RData"
    )
    message(paste("Saving", simulation_results_file))
    save(results, file = simulation_results_file)
  }
}