
library(ggplot2)
library(FactoMineR)
data <- diamonds
head(data)
data <- data[,c(1,5,6,7,8,9,10)]

pca <- PCA(data[1:3000,])
library(devtools)
devtools::load_all()
gpca <- getPCA(data[1:3000,], method="JacobiR")
print(gpca)
plot(gpca, type="individuals")
plot(gpca)

datai <- diamonds[1:3000,]
sube <- datai[gpca$Y[,1]< -2,]
subd <- datai[gpca$Y[,1]> -2,]

prop.table(rbind(table(sube[,2]),table(subd[,2])),1)
prop.table(rbind(table(sube[,3]),table(subd[,3])),1)
prop.table(rbind(table(sube[,4]),table(subd[,4])),1)
