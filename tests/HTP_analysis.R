set.seed(42)
setwd("~/OneDrive - The University of Colorado Denver/Projects/GRACKLE/")
library(tidyverse)
library(devtools)
library(ggplot2)
load_all()



file.list <- list.files(".", pattern ="HTP_results_updated_1_23", full.names = T )

results <- lapply(file.list[c(1,2,3)],function(x) get(load(x)))

results <- lapply(results, function(x) {
  x[18] <- NULL
  return(x)
})

all_results <- lapply(results, function(x) {
    grid_searches <- lapply(x, function(y) {y[[1]]})
    grackle_grids <- lapply(grid_searches, function(z) {
      z %>%
        filter(!((lambda_1 == 0 & lambda_2 ==0) | (lambda_1 == 0 & lambda_2 == 20)))
      return(z)
    })
    avg_scores <- lapply(grackle_grids, function(z) mean(z$score))
    
    top_scores <- lapply(grackle_grids, function(z) {
      z %>%
        filter(score == max(score)) %>%
        distinct(score)
    })
        
        
    
    nmf <- lapply(grid_searches, function(z) {
      z %>%
        filter(lambda_1 == 0 & lambda_2 ==0) %>%
        select(score)
    })
        
    grnmf <- lapply(x, function(y) {floor(y[[2]])})
    
    netnmf <-lapply(grid_searches, function(z) {
      z %>% 
         filter(lambda_1 == 0 & lambda_2 == 20) %>%
        select(score)
    })
      
    k <- x[[1]]$k
    
    tmp <- data.frame(
      avg_GRACKLE = unlist(avg_scores),
      top_GRACKLE = unlist(top_scores),
      nmf = unlist(nmf),
      netnmf = unlist(netnmf),
     grnmf = unlist(grnmf) 
    )
   tmp$k = k
   return(tmp)
})


top_results <- bind_rows(all_results) 

top_results <- top_results %>%
  as.data.frame() %>%
  pivot_longer(cols = -k)

top_results$name <- factor(top_results$name, levels = c("top_GRACKLE", "avg_GRACKLE", "grnmf", "netnmf", "nmf"))



p1 <- ggplot(top_results, aes(x = as.factor(k), y = value, fill = name)) +
  geom_violin(alpha = .5, position = "dodge",trim = F, scale = "width") +
  geom_boxplot(position = position_dodge(.9), alpha = .5, width= .3,outliers = FALSE) +
  theme_classic() +
  labs( x= "Rank (k)", y = "Number of enriched conditions") +
  theme(legend.position = "bottom")
ggsave(p1, file = "../GRACKLE_results/HTP_benchmarking_by_k.pdf", width = 8, height= 4,units = "in")


# p2 <- ggplot(top_results, aes(x = k, y = value, color = name, group = name)) +
#   geom_line(linewidth = 2) +
#   geom_point(size = 3) +
#   labs(title = "Number of Enriched Conditions",
#        x = "Rank",
#        y = "Enriched Conditions",
#        color = "Algorithm") +
#   theme_classic() +
#   theme(legend.position = "bottom")

#ggsave(p2, file = "../GRACKLE_results/HTP_benchmarking.pdf", width = 6, height= 4,units = "in")


## heatmap 

grid_searches <- lapply(results[[2]], function(x) {x[[1]]})
merged_grid_searches <- grid_searches %>% reduce(full_join, by = c("lambda_1", "lambda_2"))
merged_grid_searches <- merged_grid_searches %>%
  rowwise()  %>%
  mutate(avg_score = mean(c_across(-c("lambda_1", "lambda_2")), na.rm = TRUE)) %>%
  select(lambda_1,lambda_2, avg_score)
merged_grid_searches %>% filter(avg_score > 4)

p1 <- ggplot(merged_grid_searches%>% filter(lambda_2 == 0), aes(x = lambda_1, y = avg_score)) +
  geom_line(linewidth = 2) +
  geom_point(size = 3) +
  labs(title = "Number of Enriched Conditions for K = 6",
       x = "Lambda 1",
       y = "Enriched Conditions") +
  theme_classic() 
ggsave(p1, file = "../GRACKLE_results/HTP_6.pdf", width = 4, height= 5,units = "in")


 
p2 <- ggplot(merged_grid_searches, aes(x = lambda_1, y = lambda_2, fill = avg_score)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid ="white", high = "red", midpoint = mean(merged_grid_searches$avg_score)) +
  labs(title = "Rank = 4", x = "lambda_1", y = "lambda_2", fill = "enriched conditions") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(p2, file = "../GRACKLE_results/HTP_4.pdf", width = 5, height= 5,units = "in")



#### LV by condition
load("HTP_k_4_l_1_70_l_2_6.RData")
load("HTP_k_3_l_1_100_l_2_6.RData")
props <- lapply(inner_scores, function(x) {
  tmp <- x[,2] / rowSums(x)
  names(tmp) <- paste0("LV",1:length(tmp))
  return(tmp)
})

props_tab <- as.data.frame(bind_rows(props))
rownames(props_tab) <- names(props)

props_tab <- props_tab %>%
  rownames_to_column("Condition") %>%
  pivot_longer(cols = -Condition)

p3 = ggplot(props_tab, aes(x = name, y = Condition, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid ="white", high = "red", midpoint = mean(props_tab$value)) +
  labs(title = "Rank =4", x = "", y = "Condition", fill = "Proportion Cases") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(p3, file = "../GRACKLE_results/HTP_4_proportion_cases.pdf", width = 5, height= 5,units = "in")


## download gene sets
## load hallmark gene sets                                                                                                   
library(msigdbr)                                                                                
pathways <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "HPO")              





H <- read.csv("grackle_H_HTP_k_4_1_23.csv")


hypothyroidism <- pathways[which(grepl("HYPOTHYR", pathways$gs_name)),]
H_tmp <- H[order(H$LV1, decreasing = T),]
H_tmp <- H_tmp %>%
  mutate(rank = row_number())
H_tmp <- H_tmp %>%
  left_join(hypothyroidism, by = c("X"="human_gene_symbol")) %>%
  mutate(highlight = ifelse(is.na(gs_name),  NA,"Hypothyroidism"))

p1 <- ggplot(H_tmp, aes(x = rank, y = LV1)) +
  geom_point(aes(alpha = ifelse(is.na(highlight), 0.1, 1)), color = "grey", size = 3) +  # Plot all points with transparency
  geom_point(data = subset(H_tmp, !is.na(highlight)), aes(color = highlight), size = 4) +  # Plot highlighted points on top
  scale_color_manual(values = c("Hypothyroidism" = "red")) +
  scale_alpha_identity() +
  theme_classic()
ggsave(p1,file = "HTP_H_1.pdf",width = 8, height= 5,units = "in")


Sleep_Apnea <- pathways[which(grepl("SLEEP_APNEA", pathways$gs_name)),]
H_tmp <- H[order(H$LV3, decreasing = T),]
H_tmp <- H_tmp %>%
  mutate(rank = row_number())
H_tmp <- H_tmp %>%
  left_join(Sleep_Apnea, by = c("X"="human_gene_symbol")) %>%
  mutate(highlight = ifelse(is.na(gs_name),  NA,"Sleep_Apnea")) %>%
  distinct(X, .keep_all = T)

p2 <- ggplot(H_tmp, aes(x = rank, y = LV3)) +
  geom_point(aes(alpha = ifelse(is.na(highlight), 0.5, 1)), color = "grey", size = 3) +  # Plot all points with transparency
  geom_point(data = subset(H_tmp, !is.na(highlight)), aes(color = highlight), size = 4) +  # Plot highlighted points on top
  scale_color_manual(values = c("Sleep_Apnea" = "red")) +
  scale_alpha_identity() +
  theme_classic() +
  theme(legend.position = "bottom")
ggsave(p2,file = "HTP_H_3.pdf",width = 6, height= 5,units = "in")

Autoimmune <- pathways[which(grepl("AUTOIMMUNE", pathways$gs_name)),]
H_tmp <- H[order(H$LV3, decreasing = T),]
H_tmp <- H_tmp %>%
  mutate(rank = row_number())
H_tmp <- H_tmp %>%
  left_join(Autoimmune, by = c("X"="human_gene_symbol")) %>%
  mutate(highlight = ifelse(is.na(gs_name),  NA,"Autoimmune"))

p3 <- ggplot(H_tmp, aes(x = rank, y = LV3)) +
  geom_point(aes(alpha = ifelse(is.na(highlight), 0.5, 1)), color = "grey", size = 3) +  # Plot all points with transparency
  geom_point(data = subset(H_tmp, !is.na(highlight)), aes(color = highlight), size = 4) +  # Plot highlighted points on top
  scale_color_manual(values = c("Autoimmune" = "red")) +
  scale_alpha_identity() +
  theme_classic()
ggsave(p3,file = "HTP_H_3.pdf",width = 8, height= 5,units = "in")



plot(H[,2],H[,3])
heatmap(as.matrix(apply(H[,2:5],2,scale)))

TCGA <- read.csv("~/Downloads/PAM50_gene_LVs.csv")
heatmap(as.matrix(TCGA[,2:6]))
plot(TCGA[,2],TCGA[,3])
meth <- read.csv("~/Downloads/PAM50_gene_LVs_Methylation_sim.csv")
heatmap(as.matrix(meth[,2:6]))
plot(meth[,2],meth[,3])
ind <- read.csv("~/Downloads/PAM50_gene_LVs_Independent_sim.csv")
heatmap(as.matrix(ind[,2:6]))
plot(ind[,2],ind[,3])

