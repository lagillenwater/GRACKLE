% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_functions.R
\name{networkNoise}
\alias{networkNoise}
\title{networkNoise}
\usage{
networkNoise(graph, goal_modularity, membership, large_clusters)
}
\arguments{
\item{graph}{igraph graph object}

\item{goal_modularity}{The goal modularity to acheive in the permuted graph}

\item{membership}{Vertex membership of igraph object}

\item{large_clusters}{list of large_cluster ids}
}
\value{
permuted graph
}
\description{
networkNoise is a function for increasing the noise of the graph
}
