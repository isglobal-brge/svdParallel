
fsct <- function(A,i,j){
  theta <- 0.5*atan2(2*A[i,j],A[i,i]-A[j,j])
  c <- cos(theta)
  s <- sin(theta)
  t <- tan(theta)
  ans <- c(s, c, t)
  return(ans)
}

rotsct <- function(A,i,j){
  tau <- (A[i,i]-A[j,j])/(2*A[i,j])
  sol <- solve2(tau)
  t <- sol[which.min(abs(sol))]
  c <- 1/sqrt(1+t^2)
  s <- c*t  
  ans <- c(s, c, t)
  return(ans)
}

X <- matrix(rnorm(100), ncol=10)
XX <- crossprod(X)

microbenchmark::microbenchmark(fsct=fsct(XX,1,4),
                               rotsct=rotsct(XX,1,4))
