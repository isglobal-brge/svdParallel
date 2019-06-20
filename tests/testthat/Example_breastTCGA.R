# Application of the blockSVD. breastTCGA data.
# Variable under study 'breast_carcinoma_estrogen_receptor_status'.
# We are interested in which genes relate to the positive and negative strogen receptor.


library(Biobase)
library(made4)
library(brgedata)
library(MultiDataSet)

load("data/breastMulti.rda")
breastTCGA <- breastMulti[["expression"]]
breastTCGA


library(SummarizedExperiment)
genes <- assay(breastTCGA)
dim(genes)
genes[1:5,1:12]


library(FactoMineR)
res_pca <- PCA(t(genes))
plot(res_pca)

library(made4)
group<-as.factor(breastTCGA$ER.Status)
out <- ord(genes, classvec=group, type = "pca")
plot(out, nlab=3)
plotarrays(out)
par(mfrow=c(1,1))


gpca <- getPCA(t(genes), mc.cores = 5, center=T, scale=T)
plot(gpca, type="individuals")

plot(gpca$Y[group=="Positive",c(1,2)], main="Individuals",
     type="p", col="blue",xlab="PC1 (17.84%)", ylab="PC2 (7.61%)", xlim=c(-60,70), ylim=c(-70,60))
points(gpca$Y[group=="Negative",c(1,2)], col="red")
abline(a=0,b=0, lty=3)
abline(v=0, lty=3)
legend("topright", legend = c("Positve", "Negative"), col=c("blue", "red"), pch=1, cex=0.75)


mb <- microbenchmark::microbenchmark("Block method (k=2)"=getPCA(t(genes), parts = 2, center=F, scale=F),
                                     "Block method (k=5)"=getPCA(t(genes), parts = 5, center=F, scale=F),
                                     "Block method (k=7)"=getPCA(t(genes), parts = 7, center=F, scale=F),
                                     "Block method (k=10)"=getPCA(t(genes), parts = 10, center=F, scale=F),
                                     "PCA0 no method"=getPCA(t(genes), center=F, scale=F),
                                     "PCA function"=PCA(t(genes), graph=FALSE), 
                                     "ord function" = ord(genes, classvec=group, type = "pca"), times=10)

ggplot2::autoplot(mb)
