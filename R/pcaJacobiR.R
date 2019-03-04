##' PCA using Jacobi algorithm 
##'
##' Principal components and variance explained for each component 
##' using two-sided Jacobi algorithm
##' @title PCA using Jacobi algorithm
##' @param x a real nxp matrix
##' @param tol a small positive error tolerance. Default is machine tolerance
##' @param ... other parameters to be passed from JacobiR
##' @param center a logical value indicating whether the variables should be shifted to be zero centered
##' @param scale a logical value indicating whether the variables should be scaled to have unit variance  
##' @export pcaJacobiR
##' @examples
##' pcaJacobiR(iris[,-5])
##' all.equal(pcaJacobiR(iris[,-5])$components, prcomp(iris[,-5])$rotation)
##' @return a list of two components, sdev with the standard deviations of the principal components and components, the matix of the principal components.

pcaJacobiR <- function(x, tol=.Machine$double.eps, center = TRUE, scale = FALSE,...){
  if(center & !scale){
    x = scale(x, center=TRUE, scale=FALSE)
  }
  if(!center & scale){
    x = scale(x, center=FALSE, scale=TRUE)
  }
  if(center & scale){
    x = scale(x, center=TRUE, scale=TRUE)
  }
  svdX <- svdJacobiR(x,...)
  
  ans <- list(sdev = svdX$d/sqrt(nrow(x)-1), components = svdX$v)
  return(ans)
}
