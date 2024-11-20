## Script for running the simulation study based on DREAM 4 data.
library(devtools)
load_all()

library(NMF)

set.seed(42)
## aggregate the gene expression data
train <- aggregateExpression("./DREAM4/Size 10/DREAM4 training data/")
## impute the missing data with random values in the same distribution as observed values
train$expression <- noisyImputation(train$expression, .3)
original_net <- aggregateNetworks("./DREAM4/Size 10/DREAM4 gold standards/")

## split the data into train and test
dat <- split_data(train$expression, train$metadata, training_size = .7)



pat_sim <- similarity_calc(as.matrix(dat$train_metadata))

## min-max scale the input matrix
Y <- apply(dat$train_expression, 2, min_max_scale)

res <- sapply(seq(.2,1,.2), function(x){
  print(x)
 
  

  #  heatmap(as.matrix(train$expression))
  
  # 
  # train$expression <- imputeExpression(train$expression)
  # heatmap(as.matrix(train$expression))
  # 
  # train$expression[is.na(train$expression)] <- 0
  # heatmap(as.matrix(train$expression))
  ## aggregate the distinct gold standard networks
 
  net <- permuteNetwork(original_net,x)
  
    
 
 
  
  ## run GRACKLE NMF
  g_res <- GRACKLE(
    Y = Y,
    net_similarity = as.matrix(net),
    patient_similarity = pat_sim,
    diff_threshold = 1e-6,
    lambda_1 = 0,
    lambda_2 = .1,
    k = 5, 
    verbose = F,
    beta = 0)
  
  ## project the test data into W_test
  W_test <- project_W(apply(dat$test_expression,2,min_max_scale), g_res$H)
  
  ## evaluate sample loadings
  top_sample_LVs <- sampleLoadingsEvaluation(W_test,dat$test_metadata)
  
  ## evaluate gene loadings
  top_loadings <- geneLoadingsEvaluation(g_res$H, 10)
  
  ## correspondence between selected W LV's and top loading gene modules
  g_corr <- sum(unlist(lapply(1:5, function(x) {
    identical(top_sample_LVs$top[x], top_loadings$top[x])
  })))/5
  
  
  # ## plotting the matrices
  ##   heatmap(g_res$H,main = "GRACKLE")
  ##   heatmap(g_res$W, main = "GRACKLE")
  
  
  
  l_res <- nmf(Y,5,'lee', seed = 42)
  # heatmap(basis(l_res), main = "NMF")
  # heatmap(coef(l_res), main= "NMF")
  
  W_test <-  project_W(apply(dat$test_expression,2,min_max_scale), coef(l_res))
  ## evaluate sample loadings
  top_sample_LVs <- sampleLoadingsEvaluation(W_test,dat$test_metadata)
  
  ## evaluate gene loadings
  top_loadings <- geneLoadingsEvaluation(coef(l_res), 10)
  
  ## correspondence between selected W LV's and top loading gene modules
  l_corr <- sum(unlist(lapply(1:5, function(x) {
    identical(top_sample_LVs$top[x], top_loadings$top[x])
  })))/5
  
  
  
  return(list(g = g_corr, l = l_corr))
})  


toplot <- data.frame(permuted_percentage = seq(.2,1,.2), grackle = as.numeric(res[1,]), nmf = as.numeric(res[2,]))
library(tidyverse)
toplot <- toplot %>%
  pivot_longer(cols = -permuted_percentage)

ggplot(toplot, aes(x = permuted_percentage, y = value, color = name)) +
  ylim(c(0,1))+
  geom_line(linewidth = 2) +
  theme_classic()

