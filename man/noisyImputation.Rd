% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DREAM_data.R
\name{noisyImputation}
\alias{noisyImputation}
\title{noisyImputation}
\usage{
noisyImputation(expression_matrix, scaling_factor = 0.1)
}
\arguments{
\item{expression_matrix}{Data.frame of simulated gene expression data.}

\item{scaling_factor}{A numeric value to scale the imputed matrix by. The default is .1.}
}
\value{
expression_matrix Imputed matrix
}
\description{
noisyImputation is a function for imputing data with noise while maintaining the original correlation substructure
}
