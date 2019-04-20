##' SVD using parallel Jacobi algorithm 
##'
##' Eigenvalues and eigenvectores of a real symmetric matrix using 
##' two-sided Jacobi algorithm in parallel.
##' It needs doParallel library for parallelization.
##' @title SVD using Jacobi algorithm.
##' @param x a real symmetric matrix
##' @param tol a small positive error tolerance. Default is machine tolerance
##' @export parJacobi
##' @examples
##' (V <- crossprod(matrix(1:25, 5)))
##' ncores <- detectCores() - 1
##' registerDoParallel(cores=ncores)
##' cl <- makeCluster(ncores)
##' parJacobi(V)
##' stopCluster(cl)
##' @return a list of two components as for \code{base::eigen}



parJacobi <- function(x, tol=.Machine$double.eps){
  n <- nrow(x)
  matV <- diag(n)
  ls <- list()
  while(offA(x)>tol){
    for(k in 1:2){
      f <- floor((n-k)/2)
      if(f!=(n-k)/2){
        f <- f+1
      }
      
      is <- seq(1,f,1)
      js <- n-k+2-is
      for(i in 1:length(is)){
        ls[[i]] <- c(is[i],js[i])
      }
      
      res <- parLapply(cl,ls,fun=Vij, A=x, V=matV)
      #res <- mclapply(ls,FUN=ijBsct, A=x, V=V, mc.cores=detectCores())
      
      for(i in 1:length(res)){
        x <- parChanges(x,res[[i]]$ij,res[[i]]$sct)
        matV[,res[[i]]$ij[1]] <- res[[i]]$Vs$i
        matV[,res[[i]]$ij[2]] <- res[[i]]$Vs$j
      }
    }
    
    for(k in 3:(n-1)){
      f <- floor((n-k)/2)
      if(f!=(n-k)/2){
        f <- f+1
      }
      
      is <- seq(1,f,1)
      js <- n-k+2-is
      is2 <- seq(n-k+2,n-floor(k/2),1)
      js2 <- 2*n-k+2-is2
      is <- c(is, is2)
      js <- c(js, js2)
      for(i in 1:length(is)){
        ls[[i]] <- c(is[i],js[i])
      }
      res <- parLapply(cl,ls,fun=Vij, A=x, V=matV)
      #res <- mclapply(ls,FUN=ijBsct, A=x, V=V, mc.cores=detectCores())
      for(i in 1:length(res)){
        x <- parChanges(x,res[[i]]$ij,res[[i]]$sct)
        matV[,res[[i]]$ij[1]] <- res[[i]]$Vs$i
        matV[,res[[i]]$ij[2]] <- res[[i]]$Vs$j
      }
      
    }
    
    f <- floor(n/2)
    
    if(f!=n/2){
      f <- f+1
    }
    
    is3 <- seq(2,f,1)
    js3 <- n+2-is3
    for(i in 1:length(is3)){
      ls[[i]] <- c(is3[i],js3[i])
    }
    res <- parLapply(cl,ls,fun=Vij, A=x, V=matV)
    #res <- mclapply(ls,FUN=ijBsct, A=x, V=V, mc.cores=detectCores())
    for(i in 1:length(res)){
      x <- parChanges(x,res[[i]]$ij,res[[i]]$sct)
      matV[,res[[i]]$ij[1]] <- res[[i]]$Vs$i
      matV[,res[[i]]$ij[2]] <- res[[i]]$Vs$j
    }
  }
  ans <- list(values=diag(x), vectors = matV)
  ans
}