# Example "iris" data set. Jacobi function

head(iris)
iris$Species

i = iris[,-5]
head(i)

gJ <- getPCA(i,method="Jacobi")
plot(gJ)
plot(gJ, type="individuals")

library(FactoMineR)
gp <- PCA(i)

setosa <- iris$Species == "setosa"  
versicolor <- iris$Species == "versicolor"
virginica <- iris$Species == "virginica"
plot(gJ$Y[setosa,c(1,2)], main="Individuals",
     type="p", col="blue",xlab="PC1 (72.96%)", ylab="PC2 (22.85%)", xlim=c(-4,5.5))
points(gJ$Y[versicolor,c(1,2)], col="red")
points(gJ$Y[virginica,c(1,2)], col="green")
legend("bottomright", legend = c("Setosa", "Versicolor", "Virginica"), col=c("blue", "red", "green"), pch=1, cex=0.75)


