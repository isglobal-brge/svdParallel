A = matrix(c(1,2,3,4,5,6,7,8,2,9,10,11,12,13,14,15,3,10,16,17,18,19,20,21,4,11,
             17,22,23,24,25,26,5,12,18,23,27,28,29,30,6,13,19,24,28,31,32,33,7,14,20,
             25,29,32,34,35,8,15,21,26,30,33,35,36),nrow=8)

X = matrix(rnorm(90),nrow=10)

blockSVD(X, 3, 2)$d
svd(X)$d
library(irlba)
p <- parallelSVD(X,mc.cores = 3,method = "irlba")

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

