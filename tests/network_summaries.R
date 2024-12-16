library(igraph)


load("./data/Breast/directed_breast_igraph.RData")

summary(breast_g)
edge_density(breast_g)
table(as.factor(greedy_vertex_coloring(breast_g)))

