% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models.R
\name{model.sw}
\alias{model.sw}
\title{Simulate network via small world model}
\usage{
model.sw(p = 300, dim = 1, nei = 3, rw.prob = 0.05, ...)
}
\arguments{
\item{p}{The number of vertices in the graph}

\item{dim}{The dimension of the starting lattice}

\item{nei}{The neighborhood within which the vertices of the lattice will be connected}

\item{rw.prob}{The probability for rewiring}

\item{...}{Additional arguments passed to /code{igraph::sample_smallworld}}
}
\value{
An igraph object
}
\description{
Simulate network via small world model
}
