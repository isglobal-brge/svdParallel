##' @title Auxiliar functions (not to be called by the user)
##' @param x tau value
##' @param A a real symmetric matrix
##' @param V a real symmetrix matrix
##' @param i ith row
##' @param j jth column
##' @param ind vector with i and j indexes: c(i,j)
##' @param sct vector with (sin(theta),cos(theta),tan(theta)) of the Jacobi rotation



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


fsct <- function(A,i,j){
  theta <- 0.5*atan2(2*A[i,j],A[i,i]-A[j,j])
  c <- cos(theta)
  s <- sin(theta)
  t <- tan(theta)
  ans <- c(s, c, t)
  return(ans)
}


buildV <- function(A,i,j,sct){
  Ai <- A[, i]
  A[, i] <- sct[2]*Ai + sct[1]*A[, j]
  A[, j] <- -sct[1]*Ai + sct[2]*A[, j]
  
  return(V)
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

Vij <- function(ind, A, V){
  
  theta <- 0.5*atan2(2*A[ind[1],ind[2]],A[ind[1],ind[1]]-A[ind[2],ind[2]])
  c <- cos(theta)
  s <- sin(theta)
  t <- tan(theta)
  sct <- c(s, c, t)
  
  Vi <- V[, ind[1]]
  V[, ind[1]] <- sct[2]*Vi + sct[1]*V[, ind[2]]
  V[, ind[2]] <- -sct[1]*Vi + sct[2]*V[, ind[2]]
  
  Vs <- list(i=V[,ind[1]],j=V[,ind[2]])
  
  ans <- list(ij=ind, Vs=Vs, sct=sct)
  return(ans)
  
}

parChanges <- function(A,ind,sct){
  n <- nrow(A)
  B <- A  
  
  B[ind[1],] <- B[,ind[1]] <- sct[2]*A[ind[1],]+sct[1]*A[ind[2],]
  B[ind[2],] <- B[,ind[2]] <- -sct[1]*A[ind[1],]+sct[2]*A[ind[2],]
  
  B[ind[1],ind[1]]<- A[ind[1],ind[1]]+sct[3]*A[ind[1],ind[2]]
  B[ind[2],ind[2]] <- A[ind[2],ind[2]]-sct[3]*A[ind[1],ind[2]]
  
  B[ind[1],ind[2]]<- 0 
  B[ind[2],ind[1]] <- 0
  
  return(B)
}










