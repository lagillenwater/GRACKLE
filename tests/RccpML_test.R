


library(RcppML)
A <- Matrix::rsparsematrix(400, 2000, 0.1, rand.x = runif) # sparse Matrix::dgCMatrix

system.time({
    model <- RcppML::nmf(A, k = 10)
})

h0 <- predict(model, A)
evaluate(model, A) # calculate mean squared error



load( "../GRACKLE_data/data/Breast/TCGA/Breast_filtered_gene_expression_with_PAM50.RData")
library(Matrix)
matrix_exp <-Matrix(as.matrix(expression_data), sparse = T)
nnzero(matrix_exp)
sum(matrix_exp==0)


system.time({
    model <- RcppML::nmf(matrix_exp, k = 10)
})

w <- model@w
h <- model@h

system.time({
    model <- nmf(as.matrix(matrix_exp), method = "lee", rank = 5, seed = 42)
})


library(Rcpp)

A <- Matrix::rsparsematrix(2000, 2000, 0.8, rand.x = runif) # sparse Matrix::dgCMatrix

sourceCpp("./src/nmf.cpp")

model <- nmf(as.matrix(A), 10,400)
h <- model$H
w <- model$W
print(summary(w))
print(summary(t(h)))


A[1:10,1:10]
tmp <- w%*%h
tmp[1:10,1:10]


library(NMF)
model <- NMF::nmf(as.matrix(A), method = "lee", rank = 5, seed = 42)
h <- coef(model)
w <- basis(model)
print(summary(w))
print(summary(t(h)))





grn <-  crossprod(A)
pat <- tcrossprod(A)

d_g <- diag(rowSums(grn))
d_p <- diag(rowSums(pat))

l_g <- d_g - grn
l_p <- d_p - pat


model <- nmf_with_graph_regularization(as.matrix(A),as.matrix(l_p), as.matrix(l_g),10, 1,1,100)


