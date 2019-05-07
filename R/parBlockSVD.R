##' SVD using an incremental SVD algorithm in parallel.
##'
##' Singular values and left singular vectors of a real nxp matrix.
##' @title SVD using an incremental SVD algorithm.
##' @param x a real nxp matrix
##' @param mc.cores number of partitions of the matrix, number of svd computed in parallel.
##' @param ncomponents number of components to be computed. Default 2
##' @param method "svd" for complete svd, or "irlba" to compute only the first ncomponents.
##' @export parBlockSVD
##' @examples
##' (V <- (matrix(1:30, nrow=5, ncol=6)))
##' blockSVD(V,k=2)
##' all.equal(blockSVD(V,k=2, method="svd")$u, base::svd(V)$u)
##' all.equal(blockSVD(V,k=2, method="irlba")$u, base::svd(V)$u[,1:2])
##' @return a list of two components with the singular values and left singular vectors of the matrix



parBlockSVD <- function(x, ncomponents=2, mc.cores=2,
                        method="svd"){
  
  method.block <- charmatch(method, c("svd", "irlba"))
  
  svdPartial <- function(x){
    ss <- svd(x)
    ans <- sweep(ss$u, 2, FUN="*", ss$d)
    return(ans)
  }
  
  if(is.list(x) && !is.data.frame(x)){
    xt <- lapply(x,function(v) t(v)) 
    xx <- parLapply(cl,xt,fun=svdPartial)
    X <- Reduce(cbind, xx)
    
    if(method.block == 1){
      s <- svd(X)
      ans <- list(d = s$d, v = s$u)
    }
    
    if(method.block == 2){
      s <- irlba(X, nv=0, nu=ncomponents)
      ans <- list(d = s$d, v = s$u)
    }
  }
  else{
    x <- as.matrix(x)
    if(nrow(x)>=ncol(x)){
      xt <- t(x)
      p <- ncol(xt)
      cols <- seq(1,p,1)
      index <- split(cols, cut(cols,breaks=mc.cores))
      xx <- lapply(index, function(v,index) v[,index], v=xt) 
      ll <- parLapply(cl,xx, fun=svdPartial)
      X <- Reduce(cbind,ll)
      
      if(method.block == 1){
        s <- svd(X)
        ans <- list(d = s$d, v = s$u)
      }
      
      if(method.block == 2){
        s <- irlba(X, nv=0, nu=ncomponents)
        ans <- list(d = s$d, v = s$u)
      }
    }
    else{
      p <- ncol(x)
      cols <- seq(1,p,1)
      index <- split(cols, cut(cols,breaks=mc.cores))
      xx <- lapply(index, function(v,index) v[,index], v=x) 
      ll <- parLapply(cl,xx,fun=svdPartial)
      X <- Reduce(cbind,ll)
      
      if(method.block == 1){
        s <- svd(X)
        v <- crossprod(x,s$u)
        v <- sweep(v,2,s$d,FUN="/")
        ans <- list(d = s$d, v = v)
      }
      
      if(method.block == 2){
        s <- irlba(X, nv=0, nu=ncomponents)
        v <- crossprod(x,s$u)
        v <- sweep(v,2,s$d,FUN="/")
        ans <- list(d = s$d, v = v)
      }
    }
    
  }
  
  return(ans)
}