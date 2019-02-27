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
  n <- nrow(A)
  p <- ncol(A)
  if(p<n){
    A <- crossprod(x)
  }
  else{
    A <- tcrossprod(x)
  }
  jacobi <- Jacobi(A,tol)
  d <- sqrt(jacobi$values)
  v <- jacobi$vectors
  u <- A%*%v
  u <- sweep(u,2,d,FUN="/")
  ans <- list(d=d, v=v, u=u)
  return(ans)
}