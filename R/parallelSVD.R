parallelSVD <- function(x, ncomponents, k, 
                        method=c("svd", "irlba"), ...){
  
  xx <- 
    
  svdPartial <- function(x){
    ss <- svd(x)
    ans <- sweep(ss$u, 2, FUN="*", ss$d)
    return(ans)
  }  
  
  ll <- lapply(x, svd)
  X <- Reduce(ll, cbind)
  ans <- svd(X)
  ans
  
}




x1  <- matrix(rnorm(10000), nrow=1000)
x2  <- matrix(rnorm(10000), nrow=1000)

s1 <- svd(t(x1))
s2 <- svd(t(x2))


m1 <- sweep(s1$u, 2, FUN="*", s1$d)
m2 <- sweep(s2$u, 2, FUN="*", s2$d)
m <- cbind(m1,m2)

ss1 <- svd(m)

ss2 <- svd(cbind(t(x1), t(x2)))

X <- cbind(t(x1), t(x2))

ss3 <- blockSVD(X, q=1, k=nrow(x1))
