### Simulation study A

## set seed
set.seed(42)

# Assign random edge weights based on the distribution of the breast tissue network
## read in the breast adjacency matrix and expression data
breast_network <- read.csv("./data/Breast_Network.csv",row.names = 1)
breast_edge_quantiles <- quantile(as.matrix(breast_network))


# generate random nework 
g <- randomNetwork(n = 100, breast_edge_quantiles)

deg <- degree(g)
deg <- deg[order(deg)]
hist(deg)



# generate list of permuted networks
permuted_networks <- lapply(1:10, function(x) permuteNetwork(g))

# check that degree distribution in permuted networks matches
lapply(permuted_networks,  function(x){
  tmp_deg <- degree(x)
  tmp_deg <- tmp_deg[order(tmp_deg)]
  identical(deg,tmp_deg)
  
})



breast_expression <- read.csv("./data/Breast_Expression.csv", row.names = 1)


# simulate gene expression for the permuted networks
sim_exp <- lapply(permuted_networks, function(x) simulateExpression(x,breast_edge_quantiles, iterations = 50, max_expression = 3000, num_samples = 5))

## check the quantiles of the simulated expression data against observed expression data
breast_exp_quantiles <- quantile(as.matrix(breast_expression))
lapply(sim_exp, quantile)


## calculate the network edge similarity for the permuted networks
edges <- lapply(permuted_networks, E)
edges <- lapply(edges, as_ids)
weights <- lapply(permuted_networks, function(x) E(x)$op)

permuted_edge_lists <- lapply(1:10, function(x) as.data.frame(cbind(edge = edges[[x]], weight = weights[[x]])))

edge_list <- permuted_edge_lists %>% reduce(full_join, by = "edge")

edge_list[is.na(edge_list)] <- 0

edge_list <- edge_list %>%
  column_to_rownames("edge")

names(edge_list) <- paste0("network_",1:10)

network_sim <- similarityCalc(as.matrix(t(edge_list)))
heatmap(network_sim)
