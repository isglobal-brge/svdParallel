X <- matrix(rnorm(100), ncol=10)
XX <- crossprod(X)

identical(svd(XX)$d, Jacobi(XX))

microbenchmark::microbenchmark(jacobi=Jacobi(XX),
                               svd=svd(XX),
                               jacobiR=JacobiR(XX, only.values=FALSE))
