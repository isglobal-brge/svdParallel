# Application of the blockSVD. breastTCGA data.
# Variable under study 'breast_carcinoma_estrogen_receptor_status'.
# We are interested in which genes relate to the positive and negative strogen receptor.

require(snpStats)
snps <- read.plink("data/obesity") # there are three files obesity.fam, obesity.bim, obesity.bed
names(snps)