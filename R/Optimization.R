## Functions for optimization

#'  gridSearch
#'
#'  gridSearch is a function for performing a grid search to optimize the two main parameters in the GRACKLE algorithm, lambda_1 and lambda_2
#'  
#' @param lambda_1_sequence A numeric vector of variables to test for lambda_1 in the grid search. Default is seq(0,1,.1).
#' @param lambda_2_sequence A numeric vector of variables to test for lambda_2 in the grid search. Default is seq(0,1,.1).
#' @param ... Other options used to control the behavior of the GRACKLE algorithms. Passed on to [GRACKLE:GRACKLE()]
#' @param method A string representing the nmf algorithm to run. See nmf for options. Default is 'lee'.
#' @param seed Specifies the starting point or seeding method. See the NMF package for details.  Default is 42.
#' @return A list containing the following components:
#' \item{H}{The matrix \code{H} obtained from the NMF calculation.}
#' \item{W}{The matrix \code{W} obtained from the NMF calculation.}
#' @export
