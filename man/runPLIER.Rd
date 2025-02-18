% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Benchmarks.R
\name{runPLIER}
\alias{runPLIER}
\title{runPLIER}
\usage{
runPLIER(expression_matrix, k)
}
\arguments{
\item{expression_matrix}{Matrix of gene expression data (samples as rows, genes as columns)}

\item{k}{A numeric variable representing the number of latent variables (aka, rank). Default is 5.}
}
\description{
runPLIER implements a version of PLIER for comparison.
}
