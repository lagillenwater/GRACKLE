set.seed(42)
setwd("~/OneDrive - The University of Colorado Denver/Projects/GRACKLE/")
library(tidyverse)
library(devtools)
library(ggplot2)
load_all()

## load results



gnmf_results <- get(load("./results_2025_03_25/TCGA_breast_with_PAM50_results_GNMF_2025-03-25.RData"))

results <- get(load("../GRACKLE_results/TCGA_breast/TCGA_breast_with_PAM50_results.RData"))

unlist(lapply(results, function(x) x[[3]]))
unlist(lapply(gnmf_results, function(x) max(x[[1]]$score)))

grid_searches <- lapply(results, function(x) {x[[1]]})
merged_grid_searches <- grid_searches %>% reduce(full_join, by = c("lambda_1", "lambda_2"))
merged_grid_searches <- merged_grid_searches %>%
  rowwise()  %>%
  mutate(avg_score = mean(c_across(-c("lambda_1", "lambda_2")), na.rm = TRUE)) %>%
  select(lambda_1,lambda_2, avg_score)


p1 <- ggplot(merged_grid_searches, aes(x = lambda_1, y = lambda_2, fill = avg_score)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid ="white", high = "red", midpoint = mean(merged_grid_searches$avg_score)) +
  labs(title = "Alignment to BRCA PAM50 Subtypes", x = "lambda_1", y = "lambda_2", fill = "Average \n ARI") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(p1, file = "../GRACKLE_results/TCGA_breast/PAM50_parameter_heatmap.pdf", width = 5, height= 5,units = "in")


top_results <- lapply(1:length(results), function(x) {
    grid_search <- results[[x]][[1]]
  GRACKLE_score <- grid_search %>%
    filter(!(lambda_1 == 0 & lambda_2 == 0))%>%
    filter(score == max(score)) %>%
    distinct(score)
  nmf_score <- results[[x]][[2]]
  
  grnmf_score <- gnmf_results[[x]][[1]] %>%
    filter(score == max(score)) %>%
    distinct(score)
  
  netnmf_score <- grid_search %>%
    filter(lambda_1 == 0) %>%
    filter(score == max(score)) %>%
    distinct(score)
    
  tmp <- c("GRACKLE" = GRACKLE_score, "NMF" = nmf_score, "GNMF" = grnmf_score, "pr-NMF" = netnmf_score)
})


top_results <- bind_rows(top_results) 
names(top_results) <- c("GRACKLE" , "NMF" , "GNMF" , "pr-NMF" )

## load the random results
load("../GRACKLE_results/TCGA_breast/TCGA_breast_random_patients_to_PAM50_results.RData")

top_results <- cbind(top_results,randomGRACKLE =unlist(res))

top_results <- top_results %>%
  as.data.frame() %>%
  mutate(iteration = row_number()) %>%
  pivot_longer(cols = -iteration)

top_results$name <- factor(top_results$name, levels = c("GRACKLE" ,"randomGRACKLE", "NMF" , "GNMF" , "pr-NMF" ) )

p2 <- ggplot(top_results, aes(x = name, y = value)) +
  geom_violin(alpha = .5, aes(fill = name)) +
  geom_boxplot(alpha = .5, width= .1,outliers = FALSE)+
  labs(title = "Alignment to BRCA PAM50 Subtypes", x = "Algorithm", y = "ARI") +
    theme_classic() + 
theme(axis.text.x  = element_blank(),
        legend.position = "bottom")  +
  ylim(c(0,.6))


ggsave(p2, file = paste0("../GRACKLE_results/TCGA_breast/PAM50_benchmarking_", Sys.Date(), ".pdf"), width = 5, height= 5,units = "in")


### Methylation Data
gnmf_results <- get(load("./results_2025_03_25/TCGA_breast_with_PAM50_results_GNMF_2025-03-25.RData"))
results <- get(load("../GRACKLE_results/TCGA_breast/TCGA_breast_Methylation_to_PAM50_results.RData"))

unlist(lapply(results, function(x) x[[3]]))
unlist(lapply(gnmf_results, function(x) max(x[[1]]$score)))


grid_searches <- lapply(results, function(x) {x[[1]]})
merged_grid_searches <- grid_searches %>% reduce(full_join, by = c("lambda_1", "lambda_2"))
merged_grid_searches <- merged_grid_searches %>%
  rowwise()  %>%
  mutate(avg_score = mean(c_across(-c("lambda_1", "lambda_2")), na.rm = TRUE)) %>%
  select(lambda_1,lambda_2, avg_score)


p1 <- ggplot(merged_grid_searches, aes(x = lambda_1, y = lambda_2, fill = avg_score)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid ="white", high = "red", midpoint = mean(merged_grid_searches$avg_score)) +
  labs(title = "Alignment to BRCA PAM50 Subtypes", x = "lambda_1", y = "lambda_2", fill = "Average \n ARI") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(p1, file = "../GRACKLE_results/TCGA_breast/Methylation_parameter_heatmap.pdf", width = 5, height= 5,units = "in")


top_results <- lapply(1:length(results), function(x) {
  grid_search <- results[[x]][[1]]
  GRACKLE_score <- grid_search %>%
    filter(!(lambda_1 == 0 & lambda_2 == 0))%>%
    filter(score == max(score)) %>%
    distinct(score)
  nmf_score <- results[[x]][[2]]
  
  grnmf_score <- gnmf_results[[x]][[1]] %>%
    filter(score == max(score)) %>%
    distinct(score)
  
  netnmf_score <- grid_search %>%
    filter(lambda_1 == 0) %>%
    filter(score == max(score)) %>%
    distinct(score)
  
  tmp <- c("GRACKLE" = GRACKLE_score, "NMF" = nmf_score, "GNMF" = grnmf_score, "pr-NMF" = netnmf_score)
})


top_results <- bind_rows(top_results) 
names(top_results) <- c("GRACKLE" , "NMF" , "GNMF" , "pr-NMF" )

## load the random results
load("../GRACKLE_results/TCGA_breast/TCGA_breast_random_patients_Methylation_to_PAM50_results.RData")

top_results <- cbind(top_results,randomGRACKLE =unlist(res))

top_results <- top_results %>%
  as.data.frame() %>%
  mutate(iteration = row_number()) %>%
  pivot_longer(cols = -iteration)

top_results$name <- factor(top_results$name, levels = c("GRACKLE" ,"randomGRACKLE", "NMF" , "GNMF" , "pr-NMF" ) )

p2 <- ggplot(top_results, aes(x = name, y = value)) +
  geom_violin(alpha = .5, aes(fill = name)) +
  geom_boxplot(alpha = .5, width= .1,outliers = FALSE)+
  labs(title = "Alignment to BRCA PAM50 Subtypes", x = "Algorithm", y = "ARI") +
  theme_classic() + 
  theme(axis.text.x  = element_blank(),
        legend.position = "bottom")  +
  ylim(c(0,.6))


ggsave(p2, file = paste0("../GRACKLE_results/TCGA_breast/Methylation_benchmarking_", Sys.Date(), ".pdf"), width = 5, height= 5,units = "in")

## Multiple independent data
gnmf_results <- get(load("./results_2025_03_25/TCGA_breast_with_PAM50_results_GNMF_2025-03-25.RData"))
results <- get(load("../GRACKLE_results/TCGA_breast/TCGA_breast_Independent_to_PAM50_results.RData"))

unlist(lapply(results, function(x) x[[3]]))
unlist(lapply(gnmf_results, function(x) max(x[[1]]$score)))


grid_searches <- lapply(results, function(x) {x[[1]]})
merged_grid_searches <- grid_searches %>% reduce(full_join, by = c("lambda_1", "lambda_2"))
merged_grid_searches <- merged_grid_searches %>%
  rowwise()  %>%
  mutate(avg_score = mean(c_across(-c("lambda_1", "lambda_2")), na.rm = TRUE)) %>%
  select(lambda_1,lambda_2, avg_score)


p1 <- ggplot(merged_grid_searches, aes(x = lambda_1, y = lambda_2, fill = avg_score)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid ="white", high = "red", midpoint = mean(merged_grid_searches$avg_score)) +
  labs(title = "Alignment to BRCA PAM50 Subtypes", x = "lambda_1", y = "lambda_2", fill = "Average \n ARI") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(p1, file = "../GRACKLE_results/TCGA_breast/Methylation_parameter_heatmap.pdf", width = 5, height= 5,units = "in")


top_results <- lapply(1:length(results), function(x) {
  grid_search <- results[[x]][[1]]
  GRACKLE_score <- grid_search %>%
    filter(!(lambda_1 == 0 & lambda_2 == 0))%>%
    filter(score == max(score)) %>%
    distinct(score)
  nmf_score <- results[[x]][[2]]
  
  grnmf_score <- gnmf_results[[x]][[1]] %>%
    filter(score == max(score)) %>%
    distinct(score)
  
  netnmf_score <- grid_search %>%
    filter(lambda_1 == 0) %>%
    filter(score == max(score)) %>%
    distinct(score)
  
  tmp <- c("GRACKLE" = GRACKLE_score, "NMF" = nmf_score, "GNMF" = grnmf_score, "pr-NMF" = netnmf_score)
})


top_results <- bind_rows(top_results) 
names(top_results) <- c("GRACKLE" , "NMF" , "GNMF" , "pr-NMF" )

## load the random results
load("../GRACKLE_results/TCGA_breast/TCGA_breast_random_patients_Independent_to_PAM50_results.RData")

top_results <- cbind(top_results,randomGRACKLE =unlist(res))

top_results <- top_results %>%
  as.data.frame() %>%
  mutate(iteration = row_number()) %>%
  pivot_longer(cols = -iteration)

top_results$name <- factor(top_results$name, levels = c("GRACKLE" ,"randomGRACKLE", "NMF" , "GNMF" , "pr-NMF" ) )

p2 <- ggplot(top_results, aes(x = name, y = value)) +
  geom_violin(alpha = .5, aes(fill = name)) +
  geom_boxplot(alpha = .5, width= .1,outliers = FALSE)+
  labs(title = "Alignment to BRCA PAM50 Subtypes", x = "Algorithm", y = "ARI") +
  theme_classic() + 
  theme(axis.text.x  = element_blank(),
        legend.position = "bottom")  +
  ylim(c(0,.6))


ggsave(p2, file = paste0("../GRACKLE_results/TCGA_breast/Independent_benchmarking_", Sys.Date(), ".pdf"), width = 5, height= 5,units = "in")

### Random Graph

## load results

results <- get(load("../GRACKLE_results/TCGA_breast/TCGA_breast_random_patients_to_PAM50_results.RData"))

noise <- seq(.1,1,.1)
noise_results <- lapply(1:length(results),function(x) {
  tmp <- lapply(results[[x]], function(y){
  y[[1]]
  })
  data.frame(noise = noise[x], scores = unlist(tmp))
})


noise_results <- bind_rows(noise_results)

ggplot(noise_results, aes(x = as.factor(noise), y = scores))+
  geom_boxplot() +
  geom_hline(yintercept = .223) +
  theme_classic() +
  labs(x= "noise")





