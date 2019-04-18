# Application of the blockSVD. breastTCGA data.
# Variable under study 'breast_carcinoma_estrogen_receptor_status'.
# We are interested in which genes relate to the positive and negative strogen receptor.

BiocManager::install("Biobase", version = "3.8")
library(Biobase)
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
PCA(t(genes))
gpca <- getPCA(t(genes), mc.cores = 5)

