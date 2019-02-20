##' @title Auxiliar functions (not to be called by the user)
##' @param x tau value
##' @param A a real symmetric matrix
##' @param i ith row
##' @param A jth column


offA <- function(A){
  sel <- lower.tri(A)
  ans <- sum(A[sel]^2)
  return(ans)
}

solve2 <- function(x){
  ss <- sqrt(x^2 + 1)
  x1 <- -x + ss
  x2 <- -x -ss
  ans <- c(x1,x2)
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

changes <- function(A,i,j){
  n <- nrow(A)
  sct <- rotsct(A,i,j)
  B <- A  
 
  B[i,] <- B[,i] <- sct[2]*A[i,]+sct[1]*A[j,]
  B[j,] <- B[,j] <- -sct[1]*A[i,]+sct[2]*A[j,]
  
  B[i,i]<- A[i,i]+sct[3]*A[i,j]
  B[j,j] <- A[j,j]-sct[3]*A[i,j]
  
  B[i,j]<- 0 
  B[j,i] <- 0
  
  return(B)
}




