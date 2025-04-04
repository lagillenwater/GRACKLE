## Functions for processing network and expression data

## Functions for processing network and expression data

#' processNetwork
#' @description processNetwork is a function for reading in the external network and expression files and filtering for the intersection of genes. 
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @param expression_file GTEX expression file, downloaded from GRAND. (For example, see https://grand.networkmedicine.org/tissues/Breast_tissue/)
#' @param adjacency_file Network adjacency downloaded from GRAND. (For exampple, see https://grand.networkmedicine.org/tissues/Breast_tissue/)
#' @param tissue String of the tissue type for the GRN and expression data
#' #' @export

processNetwork <- function(expression_file, adjacency_file, tissue) {
  
  filtered_expression_file <- paste0("../GRACKLE_data/data/", tissue, "_filtered_expression",".RData")
  filtered_network_file <- paste0("../GRACKLE_data/data/", tissue, "_filtered_network" ,".RData")

    print(filtered_expression_file)
    
  if(file.exists(filtered_expression_file)) {
    print("found filtered files")
    load(filtered_expression_file)
    load(filtered_network_file)
  } else {
    print("reading in expression and network files")            
    expression <- read.csv(expression_file)
    network <- read.csv(adjacency_file)

    print("filtering")
    ## pivot longer
    network_long <- network %>%
      pivot_longer(cols = -X)
    
    
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

    print("saving files")
    save(filtered_network_long, file = filtered_network_file)  
    save(filtered_expression, file = filtered_expression_file)
    print("done")
  }
    print("done")
} 
  



#' alignNetwork
#' 
#' 
#' @description
#' alignNetwork is a function for reading in, filtering, and otherwise processing GRNs downloaded from GRAND for use in GRACKLE.
#' 
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @param expression_file GTEX expression file, downloaded from GRAND. (For example, see https://grand.networkmedicine.org/tissues/Breast_tissue/)
#' @param adjacency_file Network adjacency downloaded from GRAND. (For exampple, see https://grand.networkmedicine.org/tissues/Breast_tissue/)
#' @param probability_threshold Number to threshold the filtering of the graph edge weights. (Default =1)
#' @param tissue String of the tissue type for the GRN and expression data
#' @param correlation_threshold Number representing the threshold for significant correlations for signing the network. (Default = .05)
#' @export

alignNetwork <- function(expression_file, adjacency_file, probability_threshold = 1, tissue, correlation_threshold = .05){
  
  filtered_expression_file <- paste0("./data/", tissue, "_filtered_expression_prob_",probability_threshold,".RData")
  filtered_network_file <- paste0("./data/", tissue, "_filtered_network_prob_", probability_threshold,".RData")

  if(file.exists(filtered_expression_file)){
    print("found filtered files")
    load(filtered_expression_file)
    load(filtered_network_file)

  } else {
    print("reading in expression and network files")            
    expression <- read.csv(expression_file)
    network <- read.csv(adjacency_file)
    
    print("filtering")
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
    print("saving files")
    save(filtered_network_long, file = filtered_network_file)  
    save(filtered_expression, file = filtered_expression_file)
    print("done")
  } 
  
  directed_network_file = paste0("./data/", tissue, "_directed_network_", probability_threshold,".RData")
  
  if(file.exists(directed_network_file)){
    print("found directed network")
    load(directed_network_file)
  } else {
    print("creating directed network")
   #Calculate the correlation coefficients between pairs of genes to determine
   directed_network <- makeDirected(filtered_expression,filtered_network_long)
   save(directed_network, file = directed_network_file)
  }

    igraph_file = paste0( "./data/", tissue,"_igraph_prob_", probability_threshold, "_cor_",gsub("[.]", "_", as.character(correlation_threshold)),".RData")

    if(file.exists(igraph_file)) {
        print("found igraph file")
    } else {
    
        directed_network <- directed_network %>% filter(pvalue < correlation_threshold)
        print("filtering by correlation")
        
        ##   some specific modifications for sgnesR. Expression simulator only takes 1 and -1 weights. 
        directed_network$weight <- ifelse(directed_network$weight > 0,1,-1)
        directed_network <- directed_network %>%
            select(from,to,weight)
        
        
        correlations_filtered_expression <- filtered_expression %>%
            dplyr::select( union(directed_network$from, directed_network$to))
        save(correlations_filtered_expression, file = paste0( "./data/", tissue,"_correlation_filtered_expression_prob_", probability_threshold, "_cor",correlation_threshold,".RData"))


                                        # ## Turn into igraph
        g <- graph_from_data_frame(directed_network)
        save(g, file = igraph_file)

             print("saving igraph object")
    }
    print("done")
    
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

   n_targets <- length(targets)
  ## find number of cores for parallelization
#  num_cores <- 4
    correlations <- lapply(targets, function(i) {

        n <- which(targets == i)
        if(n %% 10 == 0) {    print(n/n_targets)}
        
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
    
  }) #,mc.cores = num_cores)
  
  network$correlation <-unlist(sapply(correlations, function(x) x$rho))
  
  network$pvalue <- unlist(sapply(correlations, function(x) x$p))
  network <- network %>%
    mutate(weight = ifelse(correlation > 0, probability,-probability))
  
  return(network)
}
