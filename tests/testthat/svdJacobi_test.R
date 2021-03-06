library(doParallel)
A <- matrix(1:30, nrow=6)

mat = matrix(c(1,2,3,4,5,6,7,8,2,9,10,11,12,13,14,15,3,10,16,17,18,19,20,21,4,11,
               17,22,23,24,25,26,5,12,18,23,27,28,29,30,6,13,19,24,28,31,32,33,7,14,20,
               25,29,32,34,35,8,15,21,26,30,33,35,36),nrow=8)
B <- matrix(1:50,nrow=5)


svdJ <- svdJacobi(B)
svdC <- svd(B)
parsvdJ <- parsvdJacobi(B)

svdJ$d
svdC$d
parsvdJ$d

svdJ$v
svdC$v
parsvdJ$v

svdJ$u
svdC$u

#Testing times: only Jacobi functions (for symmetric matrices):
X <- matrix(rnorm(100), ncol=10)
XX <- crossprod(X)
ncores <- detectCores() - 1
registerDoParallel(cores=ncores)
cl <- makeCluster(ncores)
microbenchmark::microbenchmark(jacobi=Jacobi(XX),
                               parjacobi=parJacobi(XX),
                               svd=svd(XX),
                               jacobiR=JacobiR(XX, only.values=TRUE))
stopCluster(cl)
#Testing times: SVD for a general matrix:

microbenchmark::microbenchmark(svdJacobi=svdJacobi(X),
                               parsvdJacobi=parsvdJacobi(X),
                               svd=svd(X))

