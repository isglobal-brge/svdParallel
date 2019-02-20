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


fsct <- function(A,i,j){
  theta <- 0.5*atan2(2*A[i,j],A[i,i]-A[j,j])
  c <- cos(theta)
  s <- sin(theta)
  t <- tan(theta)
  ans <- c(s, c, t)
  return(ans)
}

changes <- function(A,i,j,sct){
  n <- nrow(A)
  B <- A  
  
  B[i,] <- B[,i] <- sct[2]*A[i,]+sct[1]*A[j,]
  B[j,] <- B[,j] <- -sct[1]*A[i,]+sct[2]*A[j,]
  
  B[i,i]<- A[i,i]+sct[3]*A[i,j]
  B[j,j] <- A[j,j]-sct[3]*A[i,j]
  
  B[i,j]<- 0 
  B[j,i] <- 0
  
  return(B)
}

buildV <- function(V,i,j,sct){
  Vi <- V[, i]
  V[, i] <- sct[2]*Vi + sct[1]*V[, j]
  V[, j] <- -sct[1]*Vi + sct[2]*V[, j]
  
  return(V)
}





