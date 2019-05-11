##' Computation of the PCA of a matrix using function svd and the whole matrix.
##'
##' Principal components and variance explained for each component.
##' @title PCA using Jacobi algorithm or block algorithm.
##' @param x a real nxp matrix
##' @param center a logical value indicating whether the variables should be shifted to be zero centered
##' @param scale a logical value indicating whether the variables should be scaled to have unit variance  
##' @export getPCA0
##' @examples
##' getPCA0(iris[,-5])
##' @return a list with the coordinates of the variables (var.coord) and the idividuals (Y), the variance of each component and its percentage and the matrix of the principal components.



getPCA0 <- function(x, center = TRUE, scale = TRUE){
  
  
  x <- as.matrix(x)
  x.scale <- scale(x, center, scale)
  x.norm <- x.scale/sqrt(nrow(x.scale)-1)
  
  svdX <- svd(x.norm)
  
  #------------ Variables ----------------#
  variance <- (svdX$d)^2/sum((svdX$d)^2)
  var.coord <- sweep(svdX$v,2,svdX$d,FUN="*")
  rownames(var.coord) <- colnames(x)
  colnames(var.coord) <- paste("PC",1:ncol(var.coord),sep="")
  var.cos2 <- var.coord**2
  var.qual <- t(apply(var.cos2, 1 ,cumsum))
  
  
  #------------ Individuals ----------------#
  if(is.matrix(x)){
    Y <- x.norm%*%svdX$v*sqrt(nrow(x.norm))
  }
  else{
    Y <- as.matrix(Reduce(rbind,x.norm))%*%svdX$v*sqrt(nrow(Reduce(rbind,x.norm)))
  }
  colnames(Y) <- paste("PC",1:ncol(var.coord),sep="")
  ind.contr <- apply(Y**2, 2, function(v) v*100/sum(v))
  ind.iner <- apply(Y**2, 1, sum)
  ind.dist <- sqrt(ind.iner)
  ind.cos2 <- sweep(Y**2, 1, ind.iner, FUN="/")
  ind.qual <- t(apply(ind.cos2, 1, cumsum))
  
  ans <- list(varcoord=var.coord, Y=Y, var = svdX$d^2, percvar=variance, components = svdX$v)
  class(ans) <- "getPCA"
  return(ans)
}