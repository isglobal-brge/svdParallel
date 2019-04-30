# Application of the blockSVD. breastTCGA data.
# Variable under study 'breast_carcinoma_estrogen_receptor_status'.
# We are interested in which genes relate to the positive and negative strogen receptor.

BiocManager::install("Biobase", version = "3.8")
library(Biobase)
BiocManager::install("made4", version = "3.8")
library(made4)
load("breastTCGA.RData")
breastTCGA

genes <- exprs(breastTCGA)
dim(genes)
genes[1:5,1:12]

pheno <- pData(breastTCGA)
#'breast_carcinoma_estrogen_receptor_status' is variable number 12.
pheno[1:10,12]
head(pheno$breast_carcinoma_estrogen_receptor_status)
head(breastTCGA$breast_carcinoma_estrogen_receptor_status)

library(FactoMineR)
res_pca <- PCA(t(genes))
plot(res_pca)

pb <- pheno$breast_carcinoma_estrogen_receptor_status
pb <- as.factor(pb)
group<-droplevels(pbf)
out <- ord(genes, classvec=group)
plot(out, nlab=3)
res_ord <- ord(genes, type="pca")
plotarrays(res_ord)



gpca <- getPCA(t(genes), mc.cores = 5, center=T, scale=T)
plot(gpca, type="individuals")

microbenchmark::microbenchmark(getPCA(t(genes), mc.cores = 5, center=F, scale=F),
                               getPCA(t(genes), mc.cores = 10, center=F, scale=F),
                               PCA(t(genes), graph=FALSE))


positive <- pheno$breast_carcinoma_estrogen_receptor_status == "Positive"  
negative <- pheno$breast_carcinoma_estrogen_receptor_status == "Negative"
plot(gpca$Y[positive,c(1,2)], main="Individuals",
       type="p", col="blue",xlab="PC1", ylab="PC2")
points(gpca$Y[negative,c(1,2)], col="red")
