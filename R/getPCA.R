##' Computation of the PCA of a matrix using different methods: Jacobi algorithm or block algorithm.
##'
##' Principal components and variance explained for each component 
##' using two-sided Jacobi algorithm or block algorithm.
##' @title PCA using Jacobi algorithm or block algorithm.
##' @param x a real nxp matrix
##' @param method selects the method with which the function will compute the SVD. Can be: blockSVD, generalBlockSVD, Jacobi and JacobiR.
##' @param center a logical value indicating whether the variables should be shifted to be zero centered
##' @param scale a logical value indicating whether the variables should be scaled to have unit variance  
##' @param tol a small positive error tolerance. Default is machine tolerance.
##' @param ... other parameters to be passed from the function of the method used.
##' @export getPCA
##' @examples
##' getPCA(iris[,-5])
##' @return a list with the coordinates of the variables (var.coord) and the idividuals (Y), the variance of each component and its percentage and the matrix of the principal components.

getPCA <- function(x, center = TRUE, scale = TRUE, 
                   method = "blockSVD",
                   tol=.Machine$double.eps, ...){
  if(!is.matrix(x)){
    x = as.matrix(x)
  }
  
  x.scale <- scale(x, center, scale)
  x.norm <- x.scale/sqrt(nrow(x.scale)-1)
  
  method.pca <- charmatch(method, c("blockSVD", "generalBlockSVD", "Jacobi", "JacobiR"))
  if (is.na(method.pca))
    stop("type must be 'blockSVD', 'generalBlockSVD', 'Jacobi' or 'JacobiR'")
  
  if (method.pca==1) {
    svdX <- blockSVD(x.norm, ...)
  }
  if (method.pca==2) {
    svdX <- generalBlockSVD(x.norm, ...)
  }
  if (method.pca==3) {
    svdX <- svdJacobi(x.norm, ...)
  }
  if (method.pca==4) {
    svdX <- svdJacobiR(x.norm, ...)
  }
 
  #------------ Variables ----------------#
  variance <- (svdX$d)^2/sum((svdX$d)^2)
  var.coord <- sweep(svdX$v,2,svdX$d,FUN="*")
  rownames(var.coord) <- colnames(x)
  colnames(var.coord) <- paste("PC",1:ncol(var.coord),sep="")
  var.cos2 <- var.coord**2
  var.qual <- t(apply(var.cos2, 1 ,cumsum))
  
  
  #------------ Individuals ----------------#
  Y <- x%*%svdX$v*sqrt(nrow(x))
  colnames(Y) <- paste("PC",1:ncol(var.coord),sep="")
  ind.contr <- apply(Y**2, 2, function(v) v*100/sum(v))
  ind.iner <- apply(Y**2, 1, sum)
  ind.dist <- sqrt(ind.iner)
  ind.cos2 <- sweep(Y**2, 1, ind.iner, FUN="/")
  ind.qual <- t(apply(ind.cos2, 1, cumsum))
  
  ans <- list(varcoord=var.coord, Y=Y, var = svdX$d^2, percvar=variance, components = svdX$v, method = method)
  class(ans) <- "getPCA"
  return(ans)
}
