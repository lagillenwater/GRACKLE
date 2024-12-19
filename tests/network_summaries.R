library(igraph)


breast_g <- get(load("./data/breast_igraph_prob_0_cor_1.RData"))

summary(breast_g)
edge_density(breast_g)
table(as.factor(greedy_vertex_coloring(breast_g)))
hist(degree(breast_g, mode = "in"))
hist(degree(breast_g, mode = "out"))
is_connected(breast_g, mode = "weak")
mean(betweenness(breast_g, directed = TRUE, weights = NA))
hist(eigen_centrality(breast_g, directed = TRUE, weights = NA)$vector)

transitivity(breast_g, type = "global")
mean(transitivity(breast_g, type = "local"), na.rm = T)

breast_clusters <-cluster_louvain(as.undirected(breast_g),weights = NA,resolution = 1.5)
table(membership(breast_clusters))

## plot a single cluster

cluster_id <- 3

vertices_in_cluster <- which(membership(breast_clusters) == cluster_id)
vertices_in_cluster <- which(V(breast_g) %in% vertices_in_cluster)
# Create a subgraph with only the vertices in the selected cluster
subgraph <- induced_subgraph(breast_g, vids = vertices_in_cluster)

mean(betweenness(subgraph, weights = NA))
mean(closeness(subgraph, weights = NA,mode = "all"))
transitivity(subgraph, type = "global")
mean(transitivity(subgraph, type = "local"))
# Plot the subgraph
plot(
  subgraph, 
  vertex.size = 10,  # Adjust node size as needed
  vertex.label = NA,  # Hide vertex labels for clarity
  edge.width = 0.5,
  edge.arrow.size = .1,# Adjust edge width as needed,
  layout = layout_with_lgl
)



blood_g <- get(load("./data/blood_igraph_prob_1_cor_0_05.RData"))

summary(blood_g)
edge_density(blood_g)
table(as.factor(greedy_vertex_coloring(blood_g)))
hist(degree(blood_g, mode = "in"))
hist(degree(blood_g, mode = "out"))
is_connected(blood_g, mode = "weak")
summary(betweenness(blood_g, directed = TRUE, weights = NA))
hist(eigen_centrality(blood_g, directed = TRUE, weights = NA)$vector)

transitivity(blood_g,type = "global")
test <- transitivity(blood_g, type = "local")
mean(test, na.rm = T)

blood_clusters <- cluster_infomap(as.undirected( blood_g), e.weights = NA, v.weights = NULL, nb.trials = 10)
table(clusters$membership)
