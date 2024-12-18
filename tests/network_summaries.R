library(igraph)


breast_g <- get(load("./data/breast_igraph_prob_1_cor_0_05.RData"))

summary(breast_g)
edge_density(breast_g)
table(as.factor(greedy_vertex_coloring(breast_g)))
hist(degree(breast_g, mode = "in"))
hist(degree(breast_g, mode = "out"))
is_connected(breast_g, mode = "weak")
summary(betweenness(breast_g, directed = TRUE, weights = NA))
hist(eigen_centrality(breast_g, directed = TRUE, weights = NA)$vector)

breast_clusters <- cluster_louvain(as.undirected(breast_g), weights = NA)
table(membership(breast_clusters))


blood_g <- get(load("./data/blood_igraph_prob_1_cor_0_05.RData"))

summary(blood_g)
edge_density(blood_g)
table(as.factor(greedy_vertex_coloring(blood_g)))
hist(degree(blood_g, mode = "in"))
hist(degree(blood_g, mode = "out"))
is_connected(blood_g, mode = "weak")
summary(betweenness(blood_g, directed = TRUE, weights = NA))
hist(eigen_centrality(blood_g, directed = TRUE, weights = NA)$vector)

blood_clusters <- cluster_infomap(as.undirected( blood_g), e.weights = NA, v.weights = NULL, nb.trials = 10)
table(clusters$membership)
