## Script for running the simulation study based on DREAM 4 data.
library(devtools)
load_all()


net <- aggregateNetworks("./DREAM4/Size 10/DREAM4 gold standards/")
net_sim <-  similarity_calc(as.matrix(net))


## split the data into train and test
dat <- split_data(train$expression, train$metadata, training_size = .7)

set.seed(42)

res <- sapply(seq(.1,.9,.2), function(x){
  print(x)
  ## aggregate the gene expression data
  train <- aggregateExpression("./DREAM4/Size 10/DREAM4 training data/")
  
  ## impute the missing data with random values in the same distribution as observed values
  dat$train_expression <- noisyImputation(dat$train_expression, x)
  dat$test_expression <- noisyImputation(dat$test_expression, x)
  
  # 
  # train$expression <- imputeExpression(train$expression)
  # heatmap(as.matrix(train$expression))
  # 
  # train$expression[is.na(train$expression)] <- 0
  # heatmap(as.matrix(train$expression))
  ## aggregate the distinct gold standard networks
 
 
  

  pat_sim <- similarity_calc(as.matrix(dat$train_metadata))
  # pat_sim[pat_sim < 1] <- 0
   
  ## min-max scale the input matrix
  Y <- apply(dat$train_expression, 2, min_max_scale)

  i_seq <-seq(0,1,.1) 
  j_seq <-seq(0,1,.1)
  
  grid_search <- as.data.frame(expand.grid(i_seq,j_seq))
  names(grid_search) <- c("lambda_1", "lambda_2")
  grid_search$score <- 0
  
  for(i in 1:nrow(grid_search)){
    ## run GRACKLE NMF
    g_res <- GRACKLE(
        Y = Y,
        net_similarity = net_sim,
        patient_similarity = pat_sim,
        diff_threshold = 1e-6,
        lambda_1 = grid_search$lambda_1[i],
        lambda_2 = grid_search$lambda_2[i],
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
    grid_search$score[i] <- sum(unlist(lapply(1:5, function(x) {
      identical(top_sample_LVs$top[x], top_loadings$top[x])
    })))/5
  }

 print( ggplot(grid_search, aes(x = lambda_1, y= lambda_2, fill = score)) +
    geom_tile() +
    ggtitle(x))
  
  
  
  
  l_res <- nmf(Y,5,'lee', seed = 42)
  
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


toplot <- data.frame(noise = seq(.1,.9,.2),  nmf = as.numeric(res[2,]))
library(tidyverse)
toplot <- toplot %>%
  pivot_longer(cols = -noise)

ggplot(toplot, aes(x = noise, y = value, color = name)) +
  ylim(c(0,1))+
  geom_line(linewidth = 2) +
  theme_classic()

