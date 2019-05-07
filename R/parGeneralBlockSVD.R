##' SVD using an incremental SVD algorithm in parallel.
##'
##' Singular values and left singular vectors of a real nxp matrix 
##' @title SVD using an incremental SVD algorithm.
##' @param x a real nxp matrix
##' @param k number of local SVDs to concatenate at each level 
##' @param q number of levels
##' @export parGeneralBlockSVD
##' @examples
##' (V <- (matrix(1:30, nrow=5, ncol=6)))
##' generalBlockSVD(V,2,3)
##' all.equal(generalBlockSVD(V,2,3)$u, base::svd(V)$u)
##' @return a list of three components with the singular values and left and right singular vectors of the matrix

parGeneralBlockSVD <- function(A, k=2, q=1){
  
  t <- 0
  if(nrow(A)>ncol(A)){
    A <- t(A)
    t <- 1
  }
  
  n <- nrow(A)
  p <- ncol(A)
  M <- k^q
  if(!is.matrix(A))
    A <- as.matrix(A)
  if(M>p)
    stop("k^q must not be greater than the number of columns of the matrix")
  
  Ai <- list()
  cols <- seq(1,p,1)
  index <- split(cols, cut(cols,breaks=M))
  Ai <- lapply(index, function(v,index) v[,index], v=A) 
  
  for(j in 1:q){
    svdj <- parLapply(cl,Ai,fun=svd)
    Ai <- list()
    for(i in 1:(length(svdj)/k)){
      Ai[[i]] <- sweep(svdj[[(i-1)*k+1]]$u, 2, FUN="*",svdj[[(i-1)*k+1]]$d)
      for(l in 1:(k-1)){
        Ai[[i]] <- cbind(Ai[[i]],sweep(svdj[[(i-1)*k+1+l]]$u,2, FUN="*",svdj[[(i-1)*k+1+l]]$d))
      }
    }	
  }
  
  if(length(Ai)>1){
    Ai <- Reduce(cbind,Ai)
  }
  
  svdA <- svd(Ai[[1]])
  v <- crossprod(A,svdA$u)
  v <- sweep(v,2,svdA$d,FUN="/")
  
  if(t==1){
    ans <- list(d=svdA$d, u=v, v=svdA$u)
  }
  else{
    ans <- list(d=svdA$d, u=svdA$u, v=v)  
  }
  
  return(ans)
}