#Examples in parallel

ncores <- detectCores() - 1
registerDoParallel(cores=ncores)
cl <- makeCluster(ncores)
microbenchmark::microbenchmark("getPCA2"=getPCA(t(genes), mc.cores = 2, center=F, scale=F),
                               "getPCA42"=getPCA(t(genes), method="generalBlockSVD", k=4, q=2, center=F, scale=F),
                               "parblock"=agetPCA(t(genes), method = "parBlock"),
                               "parGeneral"=agetPCA(t(genes), method="parGeneral", k=4, q=2),
                               "PCA"=PCA(t(genes), graph=FALSE), times=20)
stopCluster(cl)