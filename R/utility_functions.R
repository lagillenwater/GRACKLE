## Miscelaneous utility functions


#' min_max_scale
#' 
#'  min_max_scale is a function for scaling a vector of values between 0 and 1
#'
#' @param data A data frame of numeric values
#' @param min_vals A numeric vector of minimum values for each column
#' @param max_vals A numeric vector of maximum values for each column
#' @param axis A numeric value representing the axis to scale (1 = rows, 2 = columns) (Default = 2)
#' @return scaled min max data frame
#' @export
min_max_scale <- function(data, min_vals, max_vals, axis = 2) {
  scaled_data <- sweep(data, axis, min_vals, "-")
  scaled_data <- sweep(scaled_data, axis, max_vals - min_vals, "/")
return(scaled_data)
}


#' normalize
#' 
#'  normalize is a function for scaling matrices
#'
#' @param matrix input matrix to normalize
#' @return normalized matrix
#' @export
normalize <- function(matrix) {
  return(matrix / sqrt(rowSums(matrix^2)))
}


#'  ensemblToHGNC
#' 
#'  ensemblToHGNC is a function for converting entrez ID's to HGNC ID's
#'
#' @importFrom utils download.file
#' @importFrom dplyr filter arrange select %>%
#' @param ensembl_ids  Character vector of entrez ID's to convert to HGNC gene ID's
#' @return hgnc_ids   Character vector of HGNC ID's
#' @export

ensemblToHGNC <- function(ensembl_ids) {
  hgnc_url <-  "https://storage.googleapis.com/public-download-files/hgnc/archive/archive/monthly/tsv/hgnc_complete_set_2024-10-01.txt"
  destfile <- tempfile()  # Create a temporary file
  
  # Download the file
  download.file(hgnc_url, destfile)
  
  # Read the file into a dataframe
  hgnc_df <- read.csv(destfile, header = TRUE, sep = "\t")  
  
  # filter and return ID's
  hgnc_ids <- hgnc_df %>%
    filter(ensembl_gene_id %in% ensembl_ids) %>%
    arrange(match(ensembl_gene_id,ensembl_ids)) %>%
    select(symbol,ensembl_gene_id)
    
  return(hgnc_ids)
}


#'  entrezToHGNC
#' 
#'  entrezToHGNC is a function for converting entrez ID's to HGNC ID's
#'
#' @importFrom utils download.file
#' @importFrom dplyr filter arrange select %>%
#' @param entrez_ids  Character vector of entrez ID's to convert to HGNC gene ID's
#' @return hgnc_ids   Character vector of HGNC ID's
#' @export
entrezToHGNC <- function(entrez_ids) {
  hgnc_url <-  "https://storage.googleapis.com/public-download-files/hgnc/archive/archive/monthly/tsv/hgnc_complete_set_2024-10-01.txt"
  destfile <- tempfile()  # Create a temporary file
  
  # Download the file
  download.file(hgnc_url, destfile)
  
  # Read the file into a dataframe
  hgnc_df <- read.csv(destfile, header = TRUE, sep = "\t")  
  
  # filter and return ID's
  hgnc_ids <- hgnc_df %>%
    filter(entrez_id %in% entrez_ids ) %>%
    arrange(match(entrez_id,entrez_ids)) %>%
    select(symbol,entrez_id)
    
  return(hgnc_ids)
}


