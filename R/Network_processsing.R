## Functions for processing network and expression data

#' alignNetwork
#' 
#' 
#' @description
#' alignNetwork is a function for reading in, filtering, and otherwise processing GRNs downloaded from GRAND for use in GRACKLE.
#' 
#' @param expression_file GTEX expression file, downloaded from GRAND. (For example, see https://grand.networkmedicine.org/tissues/Breast_tissue/)
#' @param adjacency_file Network adjacency downloaded from GRAND. (For exampple, see https://grand.networkmedicine.org/tissues/Breast_tissue/)
#' @param probability_threshold Number to threshold the filtering of the graph edge weights. (Default =1)
#' @param tissue String of the tissue type for the GRN and expression data
#' @return 
#' @export

alignNetwork <- function(expression_file, adjacency_file, probability_threshold = 1, tissue){
  
  filtered_expression_file <- paste0("./data/", tissue, "_filtered_expression_prob_",probability_threshold,".RData")
  filtered_network_file <- paste0("./data/", tissue, "_filtered_network_prob_", probability_threshold,".RData")

  if(file.exists(filtered_expression_file)){
    load(filtered_expression_file)
    load(filtered_network_file)
  } else {
                                           
    expression <- read.csv(expression_file)
    network <- read.csv(adjacency_file)
    
    ## pivot longer
    network_long <- network %>%
      pivot_longer(cols = -X)
    
    # ## Remove interactions below a threshold of 1 This threshold is based on the evidence of interaction, regardless of directionality.
    network_long <- network_long %>%
      filter(value >= probability_threshold)
    
    # change ensembl ID's to HGNC ID's
    network_symbol_map <- ensemblToHGNC(network_long$name)
    expression_symbol_map <- ensemblToHGNC(expression$X)
    
    network_long <- network_long %>%
      left_join(network_symbol_map, by = c("name" = "ensembl_gene_id"), relationship = "many-to-many")
    
    network_long <- network_long %>%
      filter(!is.na(symbol))
    
    network_long <- network_long %>%
      dplyr::select(X,symbol,value) %>%
      rename(from = X,to = symbol,probability = value)
    
    expression <- expression %>%
      left_join(expression_symbol_map, by = c("X" = "ensembl_gene_id"), relationship = "many-to-many")
    
    expression <- expression %>%
      filter(!is.na(symbol)) %>%
      select(-X)
    
    ## filter expression and network
    unique_genes <- unique(union(network_long$from,network_long$to))
    
    filtered_expression <- expression %>%
      filter(symbol %in% unique_genes) %>%
      column_to_rownames('symbol') %>%
      t() %>%
      as.data.frame()
    
    filtered_network_long <- network_long %>%
      filter(from %in% names(filtered_expression) & to %in% names(filtered_expression))
    
    save(filtered_breast_network_long, file = filtered_network_file)  
    save(filtered_expression, file = filtered_expression_file)
  } 
  
  directed_network_file = paste0("./data/", tissue, "_directed_network_", probability_threshold,".RData")
  
  if(file.exists(directed_network_file)){
    load(directed_network_file)
  } else {
   #Calculate the correlation coefficients between pairs of genes to determine
   directed_network <- makeDirected(filtered_expression,filtered_network_long)
   save(directed_network, file = directed_network_file)
  }
   
   # 
   # 
   # directed_network <- directed_network %>% filter(abs(correlation) > .2)
   # 
   # # some specific modifications for sgnesR. Expression simulator only takes 1 and -1 weights. 
   # directed_network$weight <- ifelse(directed_network$weight > 0,1,-1)
   # directed_network <- directed_network %>%
   #   select(from,to,weight)
   # 
   # 
   # correlations_filtered_expression <- filtered_expression %>%
   #   dplyr::select( union(directed_network$from, directed_network$to))
   # save(correlations_filtered_expression, file = "./data/Breast/correlation_filtered_expression.RData")
   # 
   # # ## Turn into igraph
   # breast_g <- graph_from_data_frame(directed_network)
   # save(breast_g, file = "./data/Breast/directed_breast_igraph.RData")
       
}


#' makeDirected
#' 
#' @description makedDirected is a function for  calculating the pairwise correlations coefficients between genes and assigning direction to the network based on correlation. 
#' @importFrom parallel detectCores  mclapply
#' @importFrom ppcor pcor
#' @importFrom dplyr filter %>% select arrange mutate
#' @param expression Data Frame of gene expression data
#' @param network Data Frame of the network edgelist with column names 'from', 'to', 'weight'.
#' @return Dataframe of network edgelist with directed weights
#' @export
makeDirected <- function(expression,network) {
  
  ## arrange network by targets
  network <- network %>%
    arrange(to)
  ## identify the target genes
  targets <- unique(network$to)
  
  ## find number of cores for parallelization
  num_cores <- detectCores() - 1
  correlations <- mclapply( targets, function(i) {
    ## filter network by targets
    tmp_network <- network %>%
      filter(to == i)
    tmp_expression <- expression %>%
      dplyr::select(c(tmp_network$from,i))
    
    ## calculate the partial correlations
    partial_correlations <- suppressWarnings(pcor(tmp_expression))
    coefs <- partial_correlations$estimate
    p_values <- partial_correlations$p.value
    
    colnames(coefs) <- colnames(p_values) <-  names(tmp_expression)
    rownames(coefs) <- rownames(p_values) <- names(tmp_expression)
        
    coef_res <- coefs  %>%
      as.data.frame  %>%
      dplyr::select(i)
    
    p_value_res <- p_values %>%
      as.data.frame() %>%
      dplyr::select(i)
    
    res <- cbind(coef_res,p_value_res)
    names(res) <- c("rho", "p")
    
    if(!(i %in% tmp_network$from)) {
      res <- res %>%
        filter(rownames(.) != i)
    }
    
    return(res)
    
  },mc.cores = num_cores)
  
  network$correlation <-unlist(sapply(correlations, function(x) x$rho))
  
  network$pvalue <- unlist(sapply(correlations, function(x) x$p))
  network <- network %>%
    mutate(weight = ifelse(correlation > 0, probability,-probability))
  
  return(network)
}
