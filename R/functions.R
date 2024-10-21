#library(matrixcalc)

GRNMF <- function(
    Y,
    input_grn,
    patient_similarity,
    diff_threshold = 1e-6,
    lambda_1 = .5,
    lambda_2 = .5,
    k = 10, 
    svd_res,
    beta = .5, 
    learning_rate = .01) {
  

  
  D_p <- diag(rowSums(patient_similarity))
  D_g <- diag(rowSums(input_grn)) 
  
  L_p <- D_p - patient_similarity
  L_g <- D_g - input_grn
  
  init_list <- initializeMatrices(Y,svd_res,k)
  
   ## Error vectors
  error_vec <- numeric()
  reconstruction_error_vec <- numeric()
  pat_sim_error_vec <- numeric()
  grn_error_vec <- numeric()
  
  H_diff_vec <- numeric()
  W_diff_vec <- numeric()
    
   W = init_list$W
   H = init_list$H
  
 # W <- matrix(runif(nrow(Y) * k), nrow = nrow(Y), ncol = k)
  #H <- matrix(runif(k * ncol(Y)), nrow = k, ncol = ncol(Y))
  
   # W <- apply(W,2,min_max_scale)
   # H <- t(apply(H,1,min_max_scale))


  
  reconstruction_error <- sum((Y-W%*%H)^2)
  pat_sim_error <- lambda_1 * sum(diag(t(W) %*% L_p %*%W))
  grn_error <-lambda_2 * sum(diag(H  %*% L_g %*% t(H)))
  
   
  # message("reconstruction error=", round(reconstruction_error,2), " | patient similarity error=", round(pat_sim_error,2), " | grn error=", round(grn_error,2))
  
  reconstruction_error_vec <- c(reconstruction_error_vec, reconstruction_error)
  pat_sim_error_vec <- c(pat_sim_error_vec, pat_sim_error)
  grn_error_vec <- c(grn_error_vec, grn_error)
  
  # norm_recon_error <- reconstruction_error_vec/max(reconstruction_error_vec)
  # norm_pat_sim_error <- pat_sim_error_vec/max(pat_sim_error_vec)
  # norm_grn_error <- grn_error_vec/max(grn_error_vec)
  
  H_diff <- 1
  W_diff <- 1
  
  
 
  while((H_diff > diff_threshold) | (W_diff >diff_threshold)) {
## for( i in 1:100){ 
    oldW <- W
    oldH <- H
    
            
    # Update W
    W <- W  * (Y %*% t(H) + lambda_1 *  L_p %*% W) / (W %*% H %*% t(H) + lambda_1 * D_p %*% W + .Machine$double.eps)

    # Update H
    H <- H  *  (t(W) %*% Y + lambda_2 * H %*% L_g ) / (t(W) %*% W %*% H + lambda_2 * H %*% D_g + .Machine$double.eps)
    
    
    reconstruction_error <- sum((Y-W%*%H)^2)
    pat_sim_error <- lambda_1 * sum(diag(t(W) %*% L_p %*%W))
    grn_error <-lambda_2 * sum(diag(H  %*% L_g %*% t(H)))
    
 #   message("reconstruction error=", round(reconstruction_error,2), " | patient similarity error=", round(pat_sim_error,2), " | grn error=", round(grn_error,2))
    
    
    reconstruction_error_vec <- c(reconstruction_error_vec, reconstruction_error)
    pat_sim_error_vec <- c(pat_sim_error_vec, pat_sim_error)
    grn_error_vec <- c(grn_error_vec, grn_error)
        
    # norm_recon_error <- reconstruction_error_vec/max(reconstruction_error_vec)
    # norm_pat_sim_error <- pat_sim_error_vec/max(pat_sim_error_vec)
    # norm_grn_error <- grn_error_vec/max(grn_error_vec)
    
    
  #   if(tail(norm_pat_sim_error,n = 1) < tail(norm_grn_error, n = 1)) {
  #     lambda_1 <- lambda_1 * (1-learning_rate)
  #   } else {
  #     
  #     lambda_2 <- lambda_2 * (1-learning_rate)
  # }
  #   print(paste("lambda_1", lambda_1))
  #   print(paste("lambda_2", lambda_2))
   
    W_diff <- round(sum((W - oldW)^2)/sum(W^2),7)
    W_diff_vec <-c(W_diff_vec, W_diff)

    H_diff <- round(sum((H - oldH)^2)/sum(H^2),7)
    H_diff_vec <-c(H_diff_vec, H_diff)
    
  }
  
  
 
    
  out <- list(residual=(Y-W%*%H), H=H, W=W,error = reconstruction_error_vec,H_diff = H_diff_vec, W_diff= W_diff_vec, pat_sim_error_vec = pat_sim_error_vec, grn_error_vec=grn_error_vec )
  
  
  return(out)
  
}

library(nnls)
project_W <- function(test_data,H_train,k) {
  # Initialize W_test
  W_test <- matrix(0, nrow=nrow(test_data), ncol=k)
  
  # Solve for W_test using NNLS
  for (i in 1:nrow(test_data)) {
    nnls_result <- nnls(t(H_train), test_data[i, ])
    W_test[i, ] <- coef(nnls_result)
  }
  
  return(W_test)
}

project_H <- function(test_data,W_train,k){
  # Initialize H_test
  H_test <- matrix(0, nrow=k, ncol=ncol(test_data))
  
  # Solve for H_test using NNLS
  for (i in 1:ncol(test_data)) {
    nnls_result <- nnls(W_train, test_data[ ,i])
    H_test[, i] <- coef(nnls_result)
  }
  
  return(H_test)
  
}

updateW <- function(W,H,Y, patient_similarity, l1, beta ) {
  degree_matrix <- diag(rowSums(patient_similarity))
  numerator <- W * (Y %*% t(H) + l1*patient_similarity %*% W)
  denominator <- W + degree_matrix%*%W + beta*W %*%H%*%t(H)
  updatedW <- (numerator/(denominator + .Machine$double.eps))
  updatedW[updatedW <0] <- 0
  return(updatedW)
}



updateH <- function(W,H,Y, input_grn, l2, beta ) {
  degree_matrix <- diag(rowSums(input_grn))
  numerator <- H * (t(t(Y) %*% W) + H %*%input_grn*l2)
  denominator <- H + H %*% degree_matrix + beta*H + t(W)%*%W%*%H
  updatedH <- (numerator/(denominator + .Machine$double.eps))
  updatedH[updatedW <0] <- 0
  return(updatedH)
}


initializeMatrices <- function(Y,svd_res,k) {
  
  # initialize matrices
  message("Initializing W & H")
  
  U <- svd_res$u
  D <- diag(svd_res$d)
  V <- svd_res$v
  
  W <- abs(U[, 1:k] %*% sqrt(D[1:k, 1:k]))
  H <- abs(sqrt(D[1:k, 1:k]) %*% t(V[, 1:k]))
  
  
  # W <- matrix(0, nrow=nrow(Y), ncol=k)
  # H <- matrix(0, nrow=k, ncol=ncol(Y))
  # 
  # 
  # W[, 1] <- sqrt(S[1]) * abs(U[, 1])
  # H[1, ] <- sqrt(S[1]) * abs(V[, 1])
  # 
  # for (i in 2:k) {
  #   u <- U[, i]
  #   v <- V[, i]
  #   u_pos <- pmax(u, 0)
  #   u_neg <- pmax(-u, 0)
  #   v_pos <- pmax(v, 0)
  #   v_neg <- pmax(-v, 0)
  #   
  #   norm_u_pos <- sqrt(sum(u_pos^2))
  #   norm_v_pos <- sqrt(sum(v_pos^2))
  #   norm_u_neg <- sqrt(sum(u_neg^2))
  #   norm_v_neg <- sqrt(sum(v_neg^2))
  #   
  #   if (norm_u_pos * norm_v_pos >= norm_u_neg * norm_v_neg) {
  #     W[, i] <- sqrt(S[i]) * u_pos / norm_u_pos
  #     H[i, ] <- sqrt(S[i]) * v_pos / norm_v_pos
  #   } else {
  #     W[, i] <- sqrt(S[i]) * u_neg / norm_u_neg
  #     H[i, ] <- sqrt(S[i]) * v_neg / norm_v_neg
  #   }
  # }
  # 
  

  return(list(W=W,H=H))
}

graphLaplacian <- function(causal_graph){
  degree_matrix <- diag(rowSums(causal_graph))
  laplacian <- degree_matrix - causal_graph
  return(laplacian)
}

findK <- function(svd_res){
 
   # Find k via variance explained
  singular_values = svd_res$d
  variance_explained = singular_values^2/sum(singular_values^2)
  cumulative_variance_explained = cumsum(variance_explained)
  k = min(which(cumulative_variance_explained > .9))
  
  message("k = ",k, "   (# = cummulative variance greater than .9")
  return(k)
}


partialCor <- function(mat) {
  cor_matrix <- cor(mat)
  
  precision_matrix <- solve(cor_matrix)
  
  partial_cor_matrix <- -precision_matrix / sqrt(diag(precision_matrix) %*% t(diag(precision_matrix)))
  
  # Set the diagonal to 1
  diag(partial_cor_matrix) <- 1
  
  
  return(partial_cor_matrix)
  
  
}

causalGraph <- function(C,W,k,percent_threshold) {
 
  cond_W <- cbind(W,  C[,condition])
  names(cond_W)[ncol(cond_W)] <- condition
  cond_partial <- partialCor(cond_W)
  
  exp_cond_W <- cbind(W,C)
  exp_cond_partial <- partialCor(exp_cond_W)
  
  cond_diff <- (abs(cond_partial[1:k,k+1]) - abs(exp_cond_partial[1:k,k+2]))/abs(cond_partial[1:k,k+1])
  
  causal_lvs <- which(cond_diff>percent_threshold)
  message("causal LVs = ", paste(causal_lvs,collapse = "  ") )
  
  cor_diff <- (abs(cond_partial[1:k,1:k]) - abs(exp_cond_partial[1:k,1:k]))/abs(cond_partial[1:k,1:k])
  
  cor_diff[cor_diff < percent_threshold] <- 0
      
  return(list(graph = cor_diff, causal_lvs = causal_lvs))
 
}

calcX <- function(C,W,k) {
        
  # CW <- Cp%*%W
  # CWr <- apply(CW,2, sum)
  # 
  # ## find the causal LVs
  # clv <- which(CWr > median(CWr))
   
  tmp_W <- cbind(C,W)
  
  cor_matrix <- cor(tmp_W)

  precision_matrix <- solve(cor_matrix)

  partial_cor_matrix <- -precision_matrix / sqrt(diag(precision_matrix) %*% t(diag(precision_matrix)))

  # Set the diagonal to 1
  diag(partial_cor_matrix) <- 1
  
  cor_diff <- (abs(cor_matrix) - abs(partial_cor_matrix))/abs(cor_matrix)
  
  causal_lvs <- cor_diff

  X <- cor_diff[,1:2]
  X[X< .5] <- 0
  
  # cor_matrix[1:2,1:2]
  # partial_cor_matrix[1:2,1:2]
  # 
#   X = matrix(0, nrow = 2, ncol = ncol(W))
  # 
  # for(i in 1:k) {
  #   single_tmp_W <-tmp_W[,c(1:2,i)] 
  # 
  #   res <- glmnet(y =  single_tmp_W[,3],x = single_tmp_W[,1:2], alpha = .9, intercept = T, standardize = F)
  # 
  #   cv_res <-  cv.glmnet(y =  single_tmp_W[,3],x = single_tmp_W[,1:2], alpha = .9, intercept = T, standardize = F)
  # 
  #   lambda_min <- cv_res$lambda.min
  # 
  #   res_coef <- coef(res, s = lambda_min)[2:3]
  #   
  #  # message(paste(round(res_coef,2), sep = " "))
  #   
  # #  if(all(res_coef > 0)){
  #     X[,i-2] <- coef(res, s = lambda_min)[2:3]
  #  #  } else {
  #  #    X[,i-2] < 0
  #  # }
  #  
  # }
  #  
  
  return(X)
  
}

# vector min_max scaling
min_max_scale <- function(x) {
  return ((x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T)))
}

# graph min_max scaling
matrix_min_max_scale <- function(matrix) {
  min_val <- min(matrix)
  max_val <- max(matrix)
  scaled_matrix <- (matrix - min_val) / (max_val - min_val)
  return(scaled_matrix)
}

# Normalize W and H such that their product remains unchanged
normalize_WH <- function(W, H) {
  norm_factors <- sqrt(colSums(W^2))
  W <- sweep(W, 2, norm_factors, "/")
  H <- sweep(H, 1, norm_factors, "*")
  return(list(W = W, H = H))
}



entrezToHGNC <- function(gene_ids,entrez_ids) {
      hgnc_ids <- gene_ids %>%
        filter(entrez_id %in% entrez_ids) %>%
        arrange(match(entrez_id,entrez_ids)) %>%
        select(symbol,entrez_id)
      return(hgnc_ids)
  }

entrezToHGNC <- function(gene_ids,entrez_ids) {
  hgnc_ids <- gene_ids %>%
    filter(entrez_id %in% entrez_ids) %>%
    arrange(match(entrez_id,entrez_ids)) %>%
    select(symbol,entrez_id)
  return(hgnc_ids)
}

download_gtex <- function(project_info) {
  data <- create_rse(project_info)
  
  metadata <- colData(data)
  
  return(data)
  
}
library(Rcpp)
cppFunction('
NumericMatrix euclideanDist(NumericMatrix X) {
  int n = X.nrow(), m = X.ncol();
  NumericMatrix out(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      double sum = 0;
      int count = 0;
      for (int k = 0; k < m; k++) {
        if (!NumericVector::is_na(X(i, k)) && !NumericVector::is_na(X(j, k))) {
          sum += pow(X(i, k) - X(j, k), 2);
          count++;
        }
      }
      if (count > 0) {
        out(i, j) = sum / count; // Normalize by the number of valid comparisons
      } else {
        out(i, j) = NA_REAL; // Assign NA if no valid comparisons
      }
      out(j, i) = out(i, j);
    }
  }
  return out;
}')


mean_absolute_deviation <- function(x) {
  mean(abs(x - mean(x)))
}


W <- results$W
project_into_W_diffproject_into_LV <- function(W,test_data, l2){
  W_hat <- (t(W)%*% test_data) /(t(W)%*%W + diag(ncol(W)))
}


library(fgsea)

GS_enrich <- function(pathways, rnk){
  simple_res <- fgsea::fgseaSimple(pathways = pathways,
                                   stats    = rnk,
                                   minSize  = 15,
                                   maxSize  = 500,
                                   nperm = 1000
  )
  
  simple_res <- simple_res %>%
    arrange(pval)
  
  return(simple_res)
}

GS_wrapper <- function(pathways ,H_matrix, sig_LVs) {
  enrich_results <- lapply(sig_LVs, function(x) {
    rnk <- scale(H_matrix[x,],center = T, scale = T)[,1]
    rnk <- rnk[order(rnk, decreasing = T)]
    GS_enrich(pathways,rnk)
  })
}

### Old functions

# causalNMF <- function(
#     Y,
#     exposure,
#     condition,
#     metadata,
#     diff_threshold,
#     percent_threshold,
#     beta,
#     lambda) {
#   
#   # Initialize W and H 
#   set.seed(123)
#   
#   # filter those with missing values for the condition
#   # set C matrix based on exposure and condition
#   C <- metadata %>%
#     select(c("LabID", exposure, condition)) %>%
#     na.omit()
#   
#   Y = Y[C$LabID,]
#   
#   C <- as.matrix(C[,c(exposure,condition)])
#   
#   # scale and center matrix
#   Y = apply(Y,2,min_max_scale)
#   Y = scale(Y, center = TRUE, scale= F )
#   
#   svd_res <- svd(Y)
#   
#   k <- findK(svd_res)
#   
#   init_list <- initializeMatrices(Y,svd_res,k)
#   
#   W = init_list$W
#   H = init_list$H
#   
#   #  message(paste0("best error based Frobenius norm = ", mean((Y-W%*%H)^2)))
#   
#   H_diff <-1
#   error_vec <- numeric()
#   H_diff_vec <- numeric()
#   
#   while(H_diff > diff_threshold) {
#     
#     causal_graph  <- causalGraph(C,W,k,percent_threshold)    
#     W <- updateW(W,H,causal_graph$graph,Y, beta = beta, lambda = lambda)
#     
#     oldH <- H
#     H <- updateH(W,H,Y,beta = beta)
#     
#     err0 <- round(sum((Y-W%*%H)^2) + sum(diag(graphLaplacian(causal_graph$graph))),4) 
#     error_vec = c(error_vec,err0)
#     H_diff <- round(sum((H - oldH)^2)/sum(H^2),7)
#     H_diff_vec <-c(H_diff_vec, H_diff)
#     
#     
#     message("error = ", err0, "   H_diff = ", H_diff)  
#     
#   }
#   
#   colnames(causal_graph$graph) <- paste0("LV",1:k)
#   rownames(causal_graph$graph) <- paste0("LV",1:k)
#   
#   out <- list(residual=(Y-W%*%H), H=H, W=W, causal_graph = causal_graph$graph,error = error_vec,H_diff = H_diff_vec, causal_lvs = causal_graph$causal_lvs )
#   
#   
#   return(out)
#   
# }
# 
# updateW <- function(W,H,causal_graph,Y,beta,lambda) {
#   degree_matrix <- diag(rowSums(causal_graph))
#   numerator <- Y %*% t(H) + lambda*W %*% causal_graph
#   denominator <- beta*W + lambda*W%*%degree_matrix + W %*%H%*%t(H)
#   updatedW <- W * (numerator/(denominator + .Machine$double.eps))
#   updatedW[updatedW <0] <- 0
#   return(updatedW)
# }
# 
# updateH <- function(W,H,Y,beta) {
#   numerator <- t(t(Y) %*% W) + H
#   denominator <- beta*H + t(W)%*%W%*%H
#   updatedH <- H * (numerator/(denominator+.Machine$double.eps))
#   updatedH[updatedH < 0] <-0
#   return(updatedH)
# }


