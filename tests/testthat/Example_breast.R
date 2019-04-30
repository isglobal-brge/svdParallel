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

gpcaL <- getPCA(data_l, center=T, scale=T)
plot(gpcaL)
plot(gpcaL, type="individuals")

#The diagnosis of breast tissues (M = malignant, B = benign)
m <- breast[,2] =="M"
b <- breast[,2] =="B"
plot(gpcaL$Y[m,c(1,2)], main="Individuals",
     type="p", col="blue",xlab="PC1 (44.27%)", ylab="PC2(18.97%)", xlim=c(-18,12), ylim=c(-14,9))
points(gpcaL$Y[b,c(1,2)], col="red")
abline(a=0,b=0, lty=3)
abline(v=0, lty=3)
legend("topright", legend = c("Malignant", "Benign"), col=c("blue", "red"), pch=1, cex=0.75)

library(FactoMineR)
PCA(num_breast)
