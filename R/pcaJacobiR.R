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
##' @return a list of two components, variance of each component and the matrix of the principal components.
##' It also plots the map of the variables and the individuals in the two first principal components.

pcaJacobiR <- function(x, tol=.Machine$double.eps, center = TRUE, scale = TRUE,...){
  if(center & !scale){
    x = scale(x, center=TRUE, scale=FALSE)
  }
  if(!center & scale){
    x = scale(x, center=FALSE, scale=TRUE)
  }
  if(center & scale){
    x = scale(x, center=TRUE, scale=TRUE)
  }
  x <- x/sqrt(nrow(x)-1)
  svdX <- svdJacobiR(x,...)
  
  #Variables:
  variance <- (svdX$d)^2/sum((svdX$d)^2)
  var.coord <- sweep(svdX$v,2,svdX$d,FUN="*")
  rownames(var.coord) <- colnames(x)
  colnames(var.coord) <- paste("PC",1:ncol(var.coord),sep="")
  var.cos2 <- var.coord**2
  var.qual <- t(apply(var.cos2, 1 ,cumsum))
  par(cex=.7)
  plot(var.coord[,1:2],xlim=c(-1,1),ylim=c(-1,1),main="Variables",
       type="n",xlab=paste("PC1 (",round(variance[1]*100,2),"%)"),ylab=paste("PC1 (",round(variance[2]*100,2),"%)"))
  abline(h=0); abline(v=0)
  noms <- abbreviate(rownames(var.coord))
  text(var.coord[,1:2],adj=c(1,-1),labels=noms)
  arrows(rep(0,10),rep(0,10),var.coord[,1],var.coord[,2],length=0.1,angle=20)
  curve(sqrt(1-x^2),-1,1,add=T)
  curve(-sqrt(1-x^2),-1,1,add=T)
  
  #individuals
  Y <- x%*%svdX$v*sqrt(nrow(x))
  colnames(Y) <- paste("PC",1:ncol(var.contr),sep="")
  ind.contr <- apply(Y**2, 2, function(v) v*100/sum(v))
  ind.iner <- apply(Y**2, 1, sum)
  ind.dist <- sqrt(ind.iner)
  ind.cos2 <- sweep(Y**2, 1, ind.iner, FUN="/")
  ind.qual <- t(apply(ind.cos2, 1, cumsum))
  par(cex=.7)
  plot(var.coord[,1:2],xlim=c(min(Y[,1])-0.5,max(Y[,1]+0.5)),ylim=c(min(Y[,2])-0.5,max(Y[,2]+0.5)), main="Individuals",
       type="n",xlab=paste("PC1 (",round(variance[1]*100,2),"%)"),ylab=paste("PC1 (",round(variance[2]*100,2),"%)"))
  abline(h=0); abline(v=0)
  points(Y[,1],Y[,2], pch=16)
  
  ans <- list(var = svdX$d^2, components = svdX$v)
  return(ans)
}
