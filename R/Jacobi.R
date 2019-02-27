##' SVD using Jacoxi algorhtm 
##'
##' Eigenvalues and eigenvectores of a real symmetric matrix using 
##' two-sided Jacoxi algorithm
##' @title SVD using Jacoxi algorhtm
##' @param x a real symmetric matrix
##' @param tol a small positive error tolerance. Default is machine tolerance
##' @export Jacobi
##' @examples
##' (V <- crossprod(matrix(1:25, 5)))
##' Jacobi(V)
##' all.equal(Jacobi(V)$values, base::eigen(V)$values)
##' @return a list of two components as for \code{base::eigen}


Jacobi <- function(x, tol=.Machine$double.eps){
  n <- nrow(x)
  
  while(offA(x)>tol){
    for(k in 1:(n-1)){
      f <- floor((n-k)/2)
      if(f!=(n-k)/2){
        f <- f+1
      }
      
      is <- seq(1,f,1)
      js <- n-k+2-is
      
      for(i in 1:length(is)){
        if(x[is[i],js[i]]!=0){
          x <- changes(x,is[i],js[i])
        }
      }
      
      if(k>2){
        is2 <- seq(n-k+2,n-floor(k/2),1)
        js2 <- 2*n-k+2-is2
        for(i in 1:length(is2)){
          if(x[is2[i],js2[i]]!=0){
            x <- changes(x,is2[i],js2[i])
          }            
        }
      }
    }
    
    f <- floor(n/2)
    
    if(f!=n/2){
      f <- f+1
    }
    
    is3 <- seq(2,f,1)
    js3 <- n+2-is3
    
    for(i in 1:length(is3)){
      if(x[is3[i],js3[i]]!=0){
        x <- changes(x,is3[i],js3[i])
      }            
    }
  }
  ans <- list(values=diag(x))
  ans
}



