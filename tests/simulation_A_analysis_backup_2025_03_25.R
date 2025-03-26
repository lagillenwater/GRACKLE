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



transitivity_analysis <- lapply(results_list, function(x) {
  tmp <- lapply(x, function(y) {
    grid_searches <- bind_rows(lapply(y, function(z) z$grid_search))
        transitivity <- unlist(lapply(y, function(z) round(z$transitivity,1)))
  grid_searches$transitivity <- transitivity        
  return(grid_searches)
  })
})

transitivity_analysis <- bind_rows(transitivity_analysis)

transitivity_analysis <- transitivity_analysis %>%
  filter(lambda_2 ==1 & lambda_1 ==0)

p1 <- ggplot(transitivity_analysis, aes(x=as.factor(transitivity), y = score)) +
  geom_boxplot()
  

gnmf_result_files <- list.files(path= "./results_2025_03_25/", pattern= "simulation_grnmf", full.names = T)

gnmf_results_list <- lapply(gnmf_result_files, function(x) get(load(x)))


full_results <- lapply(results_list, function(x) {
  
  tmp <- lapply(x, function(y){
    
    grid_searches <- lapply(y, function(z) z$grid_search)
    avg_grackle <- lapply(grid_searches, function(z) {
      z %>%
          summarize(mean(score))
    })
   top_grackle <- lapply(grid_searches, function(z) {
     z %>%
       filter(score == max(score)) %>%
       distinct(score) %>%
      .$score
   })
    nmf_results <- lapply(y, function(z) {z$nmf_res    })
    gnmf_results <- lapply(y, function(z) z$grnmf_res)
    pr_nmf_results <-  lapply(grid_searches, function(z) {
      z %>%
        filter(lambda_1 ==0) %>%
        filter(score == max(score)) %>%
        distinct(score) %>%
        .$score
    })
    
    avg_pr_nmf_results <- lapply(grid_searches, function(z) {
      z %>%
        filter(lambda_1 == 0) %>%
        summarise(mean(score))
    })
    
    modularity <- lapply(y, function(z) round(z$modularity, 1))
    transitivity <- lapply(y, function(z) round(z$transitivity,1))
    expression_noise <- lapply(y, function(z) z$expression_noise)
    
    dat <- data.frame(avg_grackle = unlist(avg_grackle), 
               top_grackle = unlist(top_grackle),
               nmf_results = unlist(nmf_results),
               gnmf_results= unlist(gnmf_results),
               avg_pr_nmf = unlist(avg_pr_nmf_results),
               pr_nmf_results = unlist(pr_nmf_results),
               modularity = unlist(modularity),
               transitivity = unlist(transitivity),
               expression_noise=unlist(expression_noise))
    
    
    return(colMeans(dat))
  })
  bind_rows(tmp)
  
})

full_results <- bind_rows(full_results)

full_results <- full_results %>%
  filter(expression_noise %in% c(.3,.5,.7))


gnmf_results <- lapply(gnmf_results_list, function(x) {
  tmp <- lapply(x, function(y) {
    top_gnmf <- lapply(y, function(z) {
      z$grnmf_res %>%
        filter(score == max(score)) %>%
        distinct(score)
    })
    avg_gnmf <- lapply(y, function(z) {
      z$grnmf_res %>%
        summarise(mean(score))
        
    })
    dat <- data.frame(
      avg_gnmf <- unlist(avg_gnmf),
      top_gnmf <- unlist(top_gnmf)
    )
    
    return(colMeans(dat))
  })
  bind_rows(tmp)
})

gnmf_results <- bind_rows(gnmf_results)
names(gnmf_results) <- c("avg_gnmf", "top_gnmf")

full_results <- full_results %>%
  select(-gnmf_results)

full_results <- cbind(full_results, gnmf_results)

top_results <- full_results %>%
  as.data.frame() %>%
  mutate(row_num = row_number()) %>%
  pivot_longer(cols = -c(row_num,transitivity,modularity,expression_noise))

top_results$name <-  factor(top_results$name, levels = c( "avg_grackle","top_grackle","avg_gnmf", "top_gnmf", "avg_pr_nmf", "pr_nmf_results","nmf_results" ))

tmp_top_results <- top_results %>%
  filter(name %in%  c( "avg_grackle","avg_gnmf", "avg_pr_nmf","nmf_results" ))
p1 <- ggplot(tmp_top_results, aes(x = as.factor(expression_noise), y = value, fill = name)) +
  geom_violin(alpha = .5, position = "dodge",trim = T, scale = "width") +
  geom_boxplot(position = position_dodge(.9), alpha = .5, width= .5,outliers = FALSE) +
  theme_classic() +
  labs( x= "Expression Noise", y = "Average Accuracy") +
  theme(legend.position = "bottom")
 
ggsave(p1, filename = "../GRACKLE_results/average_violin_overall_2025_3_25.pdf", height = 5, width = 5, units = "in")

tmp_top_results <- top_results %>%
  filter(name %in%  c( "top_grackle", "top_gnmf", "pr_nmf_results","nmf_results" ))
p1 <- ggplot(tmp_top_results, aes(x = as.factor(expression_noise), y = value, fill = name)) +
  geom_violin(alpha = .5, position = "dodge",trim = T, scale = "width") +
  geom_boxplot(position = position_dodge(.9), alpha = .5, width= .5,outliers = FALSE) +
  theme_classic() +
  labs( x= "Expression Noise", y = "Average Accuracy") +
  theme(legend.position = "bottom")

ggsave(p1, filename = "../GRACKLE_results/top_violin_overall_2025_3_25.pdf",  height = 5, width = 5, units = "in")


for(i in seq(0.3,.7,.2)){
  tmp <- tmp_top_results %>%
    filter(expression_noise == i)
  
  tmp$modularity <- factor(tmp$modularity, levels = c(.9,.8,.7,.6))
  
   p1<- ggplot(tmp, aes(x = as.factor(modularity), y = value, fill = name)) +
    geom_violin(alpha = .5, position = position_dodge(.7),trim = T, scale = "width") +
     geom_boxplot(position = position_dodge(.7), alpha = .5, width= .5,outliers = FALSE) +
    theme_classic() +
    labs( x= "Modularity", y = "Average Accuracy", title = paste("Expression noise ", i)) +
    theme(legend.position = "bottom")
  
  ggsave(p1, filename = paste0("../GRACKLE_results/mod_violin_noise_",i,"_2025_3_25.pdf"),  height = 5, width = 5, units = "in")
}


for(i in seq(.3,.7,.2)){
   tmp <- tmp_top_results %>%
    filter(expression_noise == i)
  p1<- ggplot(tmp, aes(x = as.factor(transitivity), y = value, fill = name)) +
    geom_violin(alpha = .5, position = position_dodge(.7),trim = T, scale = "width") +
    geom_boxplot(position = position_dodge(.7), alpha = .5, width= .5,outliers = FALSE) +
    theme_classic() +
    labs( x= "Transitivity", y = "Average Accuracy", title = paste("Expression noise ", i)) +
    theme(legend.position = "bottom")
  
  ggsave(p1, filename = paste0("../GRACKLE_results/trans_violin_noise_",i,"_2025_3_25.pdf"),  height = 5, width = 5, units = "in")
  
}


full_results <- lapply(results_list, function(x) {
  
  tmp <- lapply(x, function(y){
    
    grid_searches <- lapply(y, function(z) z$grid_search)
    merged_grid_searches <- grid_searches %>% reduce(full_join, by = c("lambda_1", "lambda_2"))
    merged_grid_searches <- merged_grid_searches %>%
      rowwise()  %>%
      mutate(avg_score = mean(c_across(-c("lambda_1", "lambda_2")), na.rm = TRUE)) %>%
      select(lambda_1,lambda_2, avg_score)
   
    expression_noise <- lapply(y, function(z) z$expression_noise)
    
    merged_grid_searches$expression_noise <- expression_noise[[1]]
    
    
    return(merged_grid_searches)
  })
  bind_rows(tmp)
  
})

full_results <- bind_rows(full_results)

full_results <- full_results %>%
  filter(expression_noise %in% c(.3,.5,.7))

overall <- full_results %>%
  group_by(lambda_1,lambda_2) %>%
  summarise(avg = mean(avg_score),  .groups = "keep" )


p1 <- ggplot(overall, aes(x = lambda_1, y = lambda_2, fill = avg)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid ="white", high = "red", midpoint = mean(overall$avg)) +
  labs(title = "Overall", x = "lambda_1", y = "lambda_2", fill = "Average Accuracy") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(p1, filename =  "../GRACKLE_results/overall_heatmap.pdf")


for(i in seq(.3,.7,.2)){
    tmp <- full_results %>%
      filter(expression_noise == i) %>%
    group_by(lambda_1,lambda_2) %>%
      summarise(avg = mean(avg_score),  .groups = "keep" )
    
 p1 <-  ggplot(tmp, aes(x = lambda_1, y = lambda_2, fill = avg)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid ="white", high = "red", midpoint = mean(tmp$avg)) +
      labs(title = paste("Expression noise", i) , x = "lambda_1", y = "lambda_2", fill = "Average Accuracy") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    ggsave(p1, filename = paste0("../GRACKLE_results/heatmap_noise_",i,".pdf"))
}

# merged_grid_search_scores <- lapply(noise_simulation_results, function(x){
#   grid_searches <- lapply(x, function(y) x[[1]]$grid_search)
#   merged_grid_searches <- grid_searches %>%  reduce(full_join, by = c("lambda_1", "lambda_2"))
#   merged_scores <- merged_grid_searches %>%
#     rowwise()  %>%
#     mutate(avg_score = mean(c_across(-c("lambda_1", "lambda_2")), na.rm = TRUE)) %>%
#     select(lambda_1, lambda_2, avg_score)
#                                  
# })


# names(merged_grid_search_scores) <- noise_sequence
# 
# merged_nmf_scores <- lapply(noise_simulation_results, function(x){
#   nmf_scores <- lapply(x, function(y) x[[1]]$nmf_res)
#   mean(unlist(nmf_scores), na.rm = T)
#   
# })
# 
# 
# merged_grnmf_scores <- lapply(noise_simulation_results, function(x){
#   grnmf_scores <- lapply(x, function(y) x[[1]]$grnmf_res)
#   mean(unlist(grnmf_scores), na.rm = T)
#   
# })
# 
# 
# ## plot grid search results
# lapply(1:9, function(x) {
#   p1 <- ggplot(merged_grid_search_scores[[x]], aes(x = lambda_1, y= lambda_2, fill = avg_score)) +
#     geom_tile() +
#     scale_fill_continuous(limits = c(0,1)) + 
#     theme_classic()  
#   ggsave(filename = paste0("./results/simulations/scores/heatmap_avg_grid_search_", noise_sequence[x], "_noise.pdf"), plot = p1)
# })
# 
# 
# top_grackle_scores <- lapply(merged_grid_search_scores, function(x) {
#   max(x$avg_score)
# })
# 
# 
# avg_grackle_scores <- lapply(merged_grid_search_scores, function(x) {
#  mean(x$avg_score)
# })
# 
# 
# benchmark_data <- cbind(noise_percentage = noise_sequence, TOP_GRACKLE = unlist(top_grackle_scores), AVG_GRACKLE= unlist(avg_grackle_scores), NMF = unlist(merged_nmf_scores), GRNMF = unlist(merged_grnmf_scores))
# 
# benchmark_data <- benchmark_data %>%
#   as.data.frame() %>%
#   pivot_longer(-noise_percentage, values_to = "LV_alignment", names_to = "Algorithm")
# 
# 
# p1 <- ggplot(benchmark_data, aes(x = noise_percentage, LV_alignment, color =  Algorithm)) +
#   geom_line(linewidth = 3) + 
#   theme_classic()
# ggsave("./results/simulation_A/benchmarks/simulation_1_benchmarks.pdf", p1)
# 
# 
# p2 <- ggplot(benchmark_data, aes(x = noise_percentage, LV_alignment, color =  Algorithm)) +
#   geom_line(linewidth = 3) + 
#   theme_classic()
# ggsave("./results/simulation_A/benchmarks/simulation_1_benchmarks.pdf", p1)
# 


