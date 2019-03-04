##' SVD using Jacobi algorithm 
##'
##' Singular values, right singular vectors and left singular vectors of a real nxp matrix 
##' using two-sided Jacobi algorithm
##' @title SVD using Jacobi algorithm
##' @param x a real nxp matrix
##' @param ... other parameters to be passed from JacobiR
##' @export svdJacobiR
##' @examples
##' (V <- (matrix(1:30, nrow=5, ncol=6)))
##' svdJacobiR(V)
##' all.equal(svdJacobiR(V)$v, base::svd(V)$v)
##' @return a list of three components as for \code{base::svd}


svdJacobiR <- function(x, ...){
  if(!is.matrix(x)){
    x = as.matrix(x)
  }
  n <- nrow(x)
  p <- ncol(x)
  if(p<n){
    A <- crossprod(x)
    jacobi <- JacobiR(A, ...)
    d <- sqrt(abs(jacobi$values))
    v <- jacobi$vectors
    u <- crossprod(t(x),v)
    u <- sweep(u, 2, d, FUN="/")
    ans <- list(d=d, v=v, u=u)
  }
  else{
    A <- tcrossprod(x)
    jacobi <- JacobiR(A, ...)
    d <- sqrt(abs(jacobi$values))
    u <- jacobi$vectors
    v <- crossprod(x,u)
    v <- sweep(v, 2, d, FUN="/")
    ans <- list(d=d, v=v, u=u)
  }

  return(ans)
}
