## set seed
set.seed(42)
setwd("~/OneDrive - The University of Colorado Denver/Projects/GRACKLE/")
library(tidyverse)
library(igraph)
library(devtools)
library(ggplot2)
library(ComplexHeatmap)
load_all()


result_files <- list.files(path= "../GRACKLE_results/small_simulation_A/", pattern= "simulation_results_num_nodes_400_num_groups_5", full.names = T)


results_list <- lapply(result_files, function(x) get(load(x)))

full_results <- lapply(results_list, function(x) {
  
  tmp <- lapply(x, function(y){
    grid_searches <- lapply(y, function(z) z$grid_search)
    avg_grackle <- lapply(grid_searches, function(z) mean(z$score))
    top_grackle <- lapply(grid_searches, function(z) {
      z %>%
        filter(score == max(score)) %>%
        distinct(score)
    })
    nmf_results <- lapply(y, function(z) z$nmf_res)
    grnmf_results <- lapply(y, function(z) z$grnmf_res)
    netnmf_results <-  lapply(grid_searches, function(z) {
      z %>%
        filter(lambda_1 ==0 & lambda_2 ==1) %>%
        .$score
    })
    modularity <- lapply(y, function(z) z$modularity)
    transitivity <- lapply(y, function(z) z$transitivity)
    expression_noise <- lapply(y, function(z) z$expression_noise)
    
    return(data.frame(avg_grackle = unlist(avg_grackle), 
             top_grackle = unlist(top_grackle),
             nmf_results = unlist(nmf_results),
             grnmf_results= unlist(grnmf_results),
             netnmf_results = unlist(netnmf_results),
             modularity = unlist(modularity),
             transitivity = unlist(transitivity),
             expression_noise=unlist(expression_noise)))
  })
  bind_rows(tmp)
  
})

full_results <- bind_rows(full_results)

full_results <- full_results %>%
  filter(expression_noise %in% c(.3,.5,.7))


ggplot()




merged_grid_search_scores <- lapply(noise_simulation_results, function(x){
  grid_searches <- lapply(x, function(y) x[[1]]$grid_search)
  merged_grid_searches <- grid_searches %>%  reduce(full_join, by = c("lambda_1", "lambda_2"))
  merged_scores <- merged_grid_searches %>%
    rowwise()  %>%
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


avg_grackle_scores <- lapply(merged_grid_search_scores, function(x) {
 mean(x$avg_score)
})


benchmark_data <- cbind(noise_percentage = noise_sequence, TOP_GRACKLE = unlist(top_grackle_scores), AVG_GRACKLE= unlist(avg_grackle_scores), NMF = unlist(merged_nmf_scores), GRNMF = unlist(merged_grnmf_scores))

benchmark_data <- benchmark_data %>%
  as.data.frame() %>%
  pivot_longer(-noise_percentage, values_to = "LV_alignment", names_to = "Algorithm")


p1 <- ggplot(benchmark_data, aes(x = noise_percentage, LV_alignment, color =  Algorithm)) +
  geom_line(linewidth = 3) + 
  theme_classic()
ggsave("./results/simulation_A/benchmarks/simulation_1_benchmarks.pdf", p1)


p2 <- ggplot(benchmark_data, aes(x = noise_percentage, LV_alignment, color =  Algorithm)) +
  geom_line(linewidth = 3) + 
  theme_classic()
ggsave("./results/simulation_A/benchmarks/simulation_1_benchmarks.pdf", p1)



