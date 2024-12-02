### Simulation study A

## set seed
set.seed(42)



# test <- read.csv('./data/Breast_Mammary_Tissue_TF_BindTFDF.csv')
# tfnet <- read.delim("./data/TF-RISet.tsv")
# 
# trrust <- read.delim("./data/trrust_rawdata.human.tsv", header = F)
# trrust <- trrust %>%
#   mutate(weight = ifelse(V3 == "Repression",-1,
#                          ifelse(V3 == "Activation",1,0))) %>%
#   rename(from = V1, to = V2) %>%
#   select(c(from,to,weight)) %>%
#   filter(weight != 0)
# 
# trrust_g <- graph_from_data_frame(trrust)
# ## in and out degree of breast network
# in_degrees <- degree(trrust_g,mode = "in")
# quantile(in_degrees)
# out_degrees <- degree(trrust_g,mode = "out")
# quantile(out_degrees)

# Assign random edge weights based on the distribution of in and out degree over regulatory relationships in GRAND Breast network
breast_expression <- read.csv("./data/Breast_expression.csv")
breast_network <- read.csv("./data/Breast_network.csv")

## pivot longer
breast_network_long <- breast_network %>%
  pivot_longer(cols = -X)

## Set interactions below threshold to 0
breast_network_long <- breast_network_long %>%
  mutate(value = ifelse(value < 1,0,value))

# generate random nework 
g <- degree.sequence.game(out_degrees, in_degrees, method = "simple")
sum(trrust$weight == 1)/nrow(trrust)
sm <- sample(c(1,-1), ecount(g), rep=TRUE, p=c(.62,.38))
E(g)$op <- sm

## components
comp <- components(g, mode = "strong")
community <- cluster_infomap(g)

# generate list of permuted networks
permuted_networks <- lapply(1:10, function(x) permuteNetwork(g, breast_edge_quantiles))

# check that degree distribution in permuted networks matches
lapply(permuted_networks,  function(x){
  tmp_deg <- degree(x)
  tmp_deg <- tmp_deg[order(tmp_deg)]
  identical(deg,tmp_deg)
  
})



breast_expression <- read.csv("./data/Breast_Expression.csv", row.names = 1)


# simulate gene expression for the permuted networks
sim_exp <- lapply(permuted_networks, function(x) simulateExpression(x,breast_edge_quantiles, iterations = 50, max_expression = 3000, num_samples = 5))

## combine expression data
sim_exp <- lapply(sim_exp, as.data.frame)
sim_exp <- lapply(sim_exp, function(x) {
  x %>% rownames_to_column("gene")
})
combined_exp <- sim_exp %>% reduce(full_join, by = "gene")

combined_exp <- combined_exp %>%
  column_to_rownames("gene")

names(combined_exp) <- paste0("sample", 1:ncol(combined_exp))

## install networkzoo
BiocManager::install("netZooR")


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
