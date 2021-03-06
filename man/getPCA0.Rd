% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getPCA0.R
\name{getPCA0}
\alias{getPCA0}
\title{PCA using Jacobi algorithm or block algorithm.}
\usage{
getPCA0(x, center = TRUE, scale = TRUE)
}
\arguments{
\item{x}{a real nxp matrix}

\item{center}{a logical value indicating whether the variables should be shifted to be zero centered}

\item{scale}{a logical value indicating whether the variables should be scaled to have unit variance}
}
\value{
a list with the coordinates of the variables (var.coord) and the idividuals (Y), the variance of each component and its percentage and the matrix of the principal components.
}
\description{
Computation of the PCA of a matrix using function svd and the whole matrix.
}
\details{
Principal components and variance explained for each component.
}
\examples{
getPCA0(iris[,-5])
}
