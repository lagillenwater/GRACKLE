## set seed
set.seed(42)
setwd("~/OneDrive - The University of Colorado Denver/Projects/GRACKLE/")
library(tidyverse)
library(igraph)
library(devtools)
library(ggplot2)
library(ComplexHeatmap)
load_all()

load( "./results/simulation_A/data/NMF_model_results.RData")

noise_sequence <- seq(0,.8,.1)

merged_grid_search_scores <- lapply(noise_simulation_results, function(x){
  grid_searches <- lapply(x, function(y) x[[1]]$grid_search)
  merged_grid_searches <- grid_searches %>%  reduce(full_join, by = c("lambda_1", "lambda_2"))
  merged_scores <- merged_grid_searches %>%
    rowwise() %>%
    mutate(avg_score = mean(c_across(-c("lambda_1", "lambda_2")), na.rm = TRUE)) %>%
    select(lambda_1, lambda_2, avg_score)
                                 
})

names(merged_grid_search_scores) <- noise_sequence

merged_nmf_scores <- lapply(noise_simulation_results, function(x){
  nmf_scores <- lapply(x, function(y) x[[1]]$nmf_res)
  mean(unlist(nmf_scores), na.rm = T)
  
})


merged_grnmf_scores <- lapply(noise_simulation_results, function(x){
  grnmf_scores <- lapply(x, function(y) x[[1]]$grnmf_res)
  mean(unlist(grnmf_scores), na.rm = T)
  
})


## plot grid search results
lapply(1:9, function(x) {
  p1 <- ggplot(merged_grid_search_scores[[x]], aes(x = lambda_1, y= lambda_2, fill = avg_score)) +
    geom_tile() +
    scale_fill_continuous(limits = c(0,1)) + 
    theme_classic()  
  ggsave(filename = paste0("./results/simulations/scores/heatmap_avg_grid_search_", noise_sequence[x], "_noise.pdf"), plot = p1)
})


top_grackle_scores <- lapply(merged_grid_search_scores, function(x) {
  max(x$avg_score)
})


benchmark_data <- cbind(noise_percentage = noise_sequence, GRACKLE = unlist(top_grackle_scores), NMF = unlist(merged_nmf_scores), GRNMF = unlist(merged_grnmf_scores))

benchmark_data <- benchmark_data %>%
  as.data.frame() %>%
  pivot_longer(-noise_percentage, values_to = "LV_alignment", names_to = "Algorithm")


p1 <- ggplot(benchmark_data, aes(x = noise_percentage, LV_alignment, color =  Algorithm)) +
  geom_line(linewidth = 3) + 
  theme_classic()
ggsave("./results/simulation_A/benchmarks/simulation_1_benchmarks.pdf", p1)

