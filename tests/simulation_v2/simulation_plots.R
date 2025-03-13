library(ggplot2)
library(tidyverse)
load("simulation_v2_results_summary_2.RData")

results_toplot <- results_summary %>%
    pivot_longer(cols = -c(lambda_1,lambda_2))

ggplot(results_toplot, aes( x= lambda_1, y = value, color = name)) +
    geom_violin()





tmp_toplot <- results_summary %>%
    ungroup() %>%
    select(lambda_1,  mean_train_w_ari, mean_test_w_ari) %>%
    pivot_longer(cols = -lambda_1)

ggplot(tmp_toplot, aes( x= as.factor(lambda_1), y = value, color = name)) +
    geom_violin(trim = T)


load("simulation_v2_results_k_15.RData")

tmp_toplot <- results %>%
    as.data.frame() %>%
    select(lambda_1,train_w_ari,test_w_ari) %>%
    pivot_longer(cols = -lambda_1)

ggplot(tmp_toplot, aes( x= as.factor(lambda_1), y = value, fill = name)) +
    geom_violin(trim = T)


load("simulation_v2_results_summary_k_15.RData")

results_toplot <- results_summary %>%
    select(lambda_1,lambda_2,mean_train_w_ari) %>%
    pivot_longer(cols = -c(lambda_1,lambda_2))
p1 <- ggplot(results_toplot, aes(x = as.factor(lambda_1), y = as.factor(lambda_2), fill  = value)) +
    geom_tile()+
    geom_text(aes(label = round(value,2)), color = "white") +
    scale_fill_gradient(low = "blue", high = "red") +
    theme_minimal() +
    labs(title = "Train W ARI", x = "Lambda 1", y = "Lambda 2")
results_toplot <- results_summary %>%
    select(lambda_1,lambda_2,mean_test_w_ari) %>%
    pivot_longer(cols = -c(lambda_1,lambda_2))
p2 <- ggplot(results_toplot, aes(x = as.factor(lambda_1), y = as.factor(lambda_2), fill  = value)) +
    geom_tile()+
    geom_text(aes(label = round(value,2)), color = "white") +
    scale_fill_gradient(low = "blue", high = "red") +
    theme_minimal()+
    labs(title = "Test W ARI", x = "Lambda 1", y = "Lambda 2")
results_toplot <- results_summary %>%
    select(lambda_1,lambda_2,mean_train_h_ari) %>%
    pivot_longer(cols = -c(lambda_1,lambda_2))
p3 <- ggplot(results_toplot, aes(x = as.factor(lambda_1), y = as.factor(lambda_2), fill  = value)) +
    geom_tile()+
    geom_text(aes(label = round(value,2)), color = "white") +
    scale_fill_gradient(low = "blue", high = "red") +
    theme_minimal()+
    labs(title = "Train H ARI", x = "Lambda 1", y = "Lambda 2")
library(gridExtra)
grid.arrange(p1,p2,p3)
