X <- matrix(rnorm(100), ncol=10)
XX <- crossprod(X)

identical(svd(XX)$d, Jacobi(XX)$values)

microbenchmark::microbenchmark(jacobi=Jacobi(XX),
                               svd=svd(XX),
                               jacobiR=JacobiR(XX, only.values=TRUE))

JacobiR(X)
svdJacobiR(X)
svdJacobiR(iris)

pcaJ <- pcaJacobi(iris[,-5])
pcaJR <- pcaJacobiR(iris[,-5])

all.equal(pcaJ, pcaJR, check.attributes = FALSE)

pcaJ
pcaJR
svdJR <- svdJacobiR(scale(iris[,-5],center=TRUE,scale=FALSE))
svdJR$v
