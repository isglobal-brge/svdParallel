# Application of the blockSVD. breast data.
#

breast <- read.csv("~/GitHub/svdParallel/breast.csv", quote="'")
head(breast)
ncol(breast)

num_breast <- breast[,3:32]
head(num_breast)

dim(num_breast)

data1 <- num_breast[1:200,]
data2 <- num_breast[201:400,]
data3 <- num_breast[401:569,]

data_l <- list(data1,data2,data3)
length(data_l)

gpcaL <- getPCA(data_l, center=F, scale=F)
plot(gpcaL, type="individuals")

#The diagnosis of breast tissues (M = malignant, B = benign)
m <- breast[,2] =="M"
b <- breast[,2] =="B"
plot(gpcaL$Y[m,c(1,2)], main="Individuals",
     type="p", col="blue",xlab="PC1", ylab="PC2", xlim=c(-5000,0))
points(gpcaL$Y[b,c(1,2)], col="red")
