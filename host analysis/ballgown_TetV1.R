#set working directory first

#load library
library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)
#library(RskittleBrewer) not available for R v 4.0.2, this is use for setting up colors

#create the phenotype data for each sample
pheno_data = read.csv("TetV1_pheno_edited.csv")

#load ballgown data structures for each sample
bg_TetV1 = ballgown(dataDir = "ballgown", samplePattern = "inf", pData=pheno_data)

#filter to remove low-abundance genes
bg_TetV1_filt = subset(bg_TetV1,"rowVars(texpr(bg_TetV1)) >1",genomesubset=TRUE)

#export the transcript FPKM in a csv file
TetV1_trans_FPKM <- bg_TetV1_filt@expr[["trans"]]
write.csv(TetV1_trans_FPKM, 'TetV1_trans_FPKM.csv')
