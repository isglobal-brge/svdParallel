##' PCA using Jacobi algorithm 
##'
##' Principal components and variance explained for each component 
##' using two-sided Jacobi algorithm
##' @title PCA using Jacobi algorithm
##' @param x a real nxp matrix
##' @param center a logical value indicating whether the variables should be shifted to be zero centered
##' @param scale a logical value indicating whether the variables should be scaled to have unit variance  
##' @param tol a small positive error tolerance. Default is machine tolerance
##' @param ... other parameters to be passed from JacobiR
##' @export pcaJacobiR
##' @examples
##' pcaJacobiR(iris[,-5])
##' @return a list of two components, variance of each component and the matrix of the principal components.
##' It also plots the map of the variables and the individuals in the two first principal components.

getPCA <- function(x, center = TRUE, scale = TRUE, 
                   method = "Jacobi",
                   tol=.Machine$double.eps, ...){
  
  x.scale <- scale(x, center, scale)
  x.norm <- x.scale/sqrt(nrow(x.scale)-1)
  
  
  svdX <- svdJacobiR(x.norm, ...)
 
  #------------ Variables ----------------#
  variance <- (svdX$d)^2/sum((svdX$d)^2)
  var.coord <- sweep(svdX$v,2,svdX$d,FUN="*")
  rownames(var.coord) <- colnames(x)
  colnames(var.coord) <- paste("PC",1:ncol(var.coord),sep="")
  var.cos2 <- var.coord**2
  var.qual <- t(apply(var.cos2, 1 ,cumsum))
  
  
  #------------ Individuals ----------------#
  Y <- x%*%svdX$v*sqrt(nrow(x))
  colnames(Y) <- paste("PC",1:ncol(var.contr),sep="")
  ind.contr <- apply(Y**2, 2, function(v) v*100/sum(v))
  ind.iner <- apply(Y**2, 1, sum)
  ind.dist <- sqrt(ind.iner)
  ind.cos2 <- sweep(Y**2, 1, ind.iner, FUN="/")
  ind.qual <- t(apply(ind.cos2, 1, cumsum))
  
  ans <- list(var.corrd, Y=Y, var = svdX$d^2, components = svdX$v)
  class(ans) <- "getPCA"
  return(ans)
}
