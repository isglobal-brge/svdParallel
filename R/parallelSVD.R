##' SVD using an incremental SVD algorithm.
##'
##' Singular values and left singular vectors of a real nxp matrix 
##' @title SVD using an incremental SVD algorithm.
##' @param x a real nxp matrix
##' @param k number of partitions of the matrix, number of svd computed in parallel.
##' @param ncomponents number of components to be computed. Default 2
##' @param method "svd" for complete svd, or "irlba" to compute only the first ncomponents.
##' @export parallelSVD
##' @examples
##' (V <- (matrix(1:30, nrow=5, ncol=6)))
##' parallelSVD(V,k=2)
##' all.equal(parallelSVD(V,k=2, method="svd")$u, base::svd(V)$u)
##' all.equal(parallelSVD(V,k=2, method="irlba")$u, base::svd(V)$u[,1:2])
##' @return a list of two components with the singular values and left singular vectors of the matrix


parallelSVD <- function(x, ncomponents=2, mc.cores,
                        method=c("svd", "irlba")){
  
  if(charmatch(method, c("svd", "irlba")) == 1){
    p <- ncol(x)
    cols <- seq(1,p,1)
    index <- split(cols, cut(cols,breaks=mc.cores))
    xx <- lapply(index, function(v,index) v[,index], v=x) 
    
    svdPartial <- function(x){
      ss <- svd(x)
      ans <- sweep(ss$u, 2, FUN="*", ss$d)
      return(ans)
    }  
    
    ll <- lapply(xx, svdPartial)
    X <- Reduce(cbind,ll)
    s <- svd(X)
    ans <- list(d = s$d, u = s$u)
    ans
  }
  if(charmatch(method, c("svd", "irlba")) == 2){
    partitions <- 5
    p <- ncol(x)
    cols <- seq(1,p,1)
    index <- split(cols, cut(cols,breaks=mc.cores))
    xx <- lapply(index, function(v,index) v[,index], v=x) 
    
    svdPartial <- function(x){
      ss <- svd(x)
      ans <- sweep(ss$u, 2, FUN="*", ss$d)
      return(ans)
    }  
    
    ll <- lapply(xx, svdPartial)
    X <- Reduce(cbind,ll)
    s <- irlba(X, nv=0, nu=ncomponents)
    ans <- list(d = s$d, u = s$u)
    ans
  }
  return(ans)
}

