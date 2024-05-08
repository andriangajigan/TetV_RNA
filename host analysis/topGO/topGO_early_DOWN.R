# load topGO
library("topGO")
sessionInfo()

# read in the 'gene universe' file 
# which in this case the trinotate_tet_host_GO a gene to GO mapping output of Trinotate
geneID2GO <- readMappings(file = "trinotate_tet_host_GO.txt")
geneUniverse <- names(geneID2GO)

# read in the genes of interest from DE analysis
genesOfInterest <- read.table("DE_early_DOWN.txt",header=FALSE)
genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

# build the GOdata object in topGO
myGOdata_BP <- new("topGOdata", description="Early Stage DOWN BP", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
myGOdata_BP

# run the Fisher's exact tests
resultClassic_BP <- runTest(myGOdata_BP, algorithm="classic", statistic="fisher")
resultClassic_BP
resultElim_BP <- runTest(myGOdata_BP, algorithm="elim", statistic="fisher")
resultElim_BP
resultTopgo_BP <- runTest(myGOdata_BP, algorithm="weight01", statistic="fisher")
resultTopgo_BP
resultParentchild_BP <- runTest(myGOdata_BP, algorithm="parentchild", statistic="fisher")
resultParentchild_BP

# see how many results we get where weight01 gives a P-value <= 0.05
mysummary_BP <- summary(attributes(resultTopgo_BP)$score <= 0.05)
numsignif_BP <- as.integer(mysummary_BP[[3]]) # how many terms is it true that P <= 0.05

# print out results:
allRes_BP <- GenTable(myGOdata_BP, classicFisher = resultClassic_BP, elimFisher = resultElim_BP, topgoFisher = resultTopgo_BP, parentchildFisher = resultParentchild_BP, orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif_BP)
write.table(allRes_BP, file="topGO_early_DOWN_BP.txt", append=F, quote=F, sep="\t", eol="\n")

# for CC
myGOdata_CC <- new("topGOdata", description="Early Stage DOWN CC", ontology="CC", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
myGOdata_CC
resultClassic_CC <- runTest(myGOdata_CC, algorithm="classic", statistic="fisher")
resultClassic_CC
resultElim_CC <- runTest(myGOdata_CC, algorithm="elim", statistic="fisher")
resultElim_CC 
resultTopgo_CC <- runTest(myGOdata_CC, algorithm="weight01", statistic="fisher")
resultTopgo_CC
resultParentchild_CC <- runTest(myGOdata_CC, algorithm="parentchild", statistic="fisher")
resultParentchild_CC
mysummary_CC <- summary(attributes(resultTopgo_CC)$score <= 0.05)
numsignif_CC <- as.integer(mysummary_CC[[3]]) # how many terms is it true that P <= 0.05
allRes_CC <- GenTable(myGOdata_CC, classicFisher = resultClassic_CC, elimFisher = resultElim_CC, topgoFisher = resultTopgo_CC, parentchildFisher = resultParentchild_CC, orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif_CC)
write.table(allRes_CC, file="topGO_early_DOWN_CC.txt", append=F, quote=F, sep="\t", eol="\n")

#for MF
myGOdata_MF <- new("topGOdata", description="Early Stage DOWN MF", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
myGOdata_MF
resultClassic_MF <- runTest(myGOdata_MF, algorithm="classic", statistic="fisher")
resultClassic_MF
resultElim_MF <- runTest(myGOdata_MF, algorithm="elim", statistic="fisher")
resultElim_MF
resultTopgo_MF <- runTest(myGOdata_MF, algorithm="weight01", statistic="fisher")
resultTopgo_MF
resultParentchild_MF <- runTest(myGOdata_MF, algorithm="parentchild", statistic="fisher")
resultParentchild_MF
mysummary_MF <- summary(attributes(resultTopgo_MF)$score <= 0.05)
numsignif_MF <- as.integer(mysummary_MF[[3]]) # how many terms is it true that P <= 0.05
allRes_MF <- GenTable(myGOdata_MF, classicFisher = resultClassic_MF, elimFisher = resultElim_MF, topgoFisher = resultTopgo_MF, parentchildFisher = resultParentchild_MF, orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif_MF)
write.table(allRes_MF, file="topGO_early_DOWN_MF.txt", append=F, quote=F, sep="\t", eol="\n")
