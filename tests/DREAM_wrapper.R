## Script for running the simulation study based on DREAM 4 data.
library(devtools)
load_all()

## aggregate the gene expression data
train <- aggregateExpression("/DREAM4/Size 10/DREAM4 training data/")

## impute the missing data with random values in the same distribution as observed values
train$expression <- imputeExpression(train$expression)

## aggregate the distinct gold standard networks
net <- aggregateNetworks("/DREAM4/Size 10/DREAM4 gold standards/")

## split the data into train and test
dat <- split_data(train$expression, train$metadata, training_size = .7)

## calculate the similarity matrices
#net_sim <- similarity_calc(as.matrix(net))
pat_sim <- similarity_calc(as.matrix(dat$train_metadata))

## min-max scale the input matrix
Y <- apply(dat$train_expression, 2, min_max_scale)

## run GRACKLE NMF
g_res <- GRACKLE(
    Y = Y,
    net_similarity = as.matrix(net),
    patient_similarity = pat_sim,
    diff_threshold = 1e-5,
    lambda_1 = 1,
    lambda_2 = 1,
    k = 10, 
    verbose = F)
  
## project the test data into W_test
W_test <- project_W(dat$test_expression, g_res$H, 10)

## test for predictive accuracy
coefficients <- categoricalPrediction(W_test, dat$test_metadata)
top_coefficients <- coefficients$top
## evaluate gene loadings
top_loadings <- geneLoadingsEvaluation(g_res$H, 10)


