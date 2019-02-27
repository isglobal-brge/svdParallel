##' SVD using Jacobi algorithm 
##'
##' Singular values, right singular vectors and left singular vectors of a real nxp matrix 
##' using two-sided Jacobi algorithm
##' @title SVD using Jacobi algorithm
##' @param x a real nxp matrix
##' @param tol a small positive error tolerance. Default is machine tolerance
##' @export svdJacobi
##' @examples
##' (V <- (matrix(1:30, nrow=5, ncol=6)))
##' svdJacobi(V)
##' all.equal(Jacobi(V)$v, base::svd(V)$v)
##' @return a list of three components as for \code{base::svd}


svdJacobi <- function(x, tol=.Machine$double.eps){
  n <- nrow(x)
  p <- ncol(x)
  if(p<n){
    A <- crossprod(x)
    jacobi <- Jacobi(A,tol)
    d <- sqrt(abs(jacobi$values))
    v <- jacobi$vectors
    u <- x%*%v
    u <- sweep(u,2,d,FUN="/")
    ans <- list(d=d, v=v, u=u)
  }
  else{
    A <- tcrossprod(x)
    jacobi <- Jacobi(A,tol)
    d <- sqrt(abs(jacobi$values))
    u <- jacobi$vectors
    v <- t(x)%*%u
    v <- sweep(v,2,d,FUN="/")
    ans <- list(d=d, v=v, u=u)
  }

  return(ans)
}