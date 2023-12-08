## ---------------------------
##
## This script was created with R version 4.2.1 (2022-06-23)
##
##
## Title: Personalized identification and characterization of genome-wide gene expression differences 
##        between patient-matched intracranial and extracranial melanoma metastasis pairs
##
## Purpose of script: Script to rebuild analysis of manuscript
##
## Author: Theresa Kraft
##
## Date created: 06th July 2023
## Date last modified: 07th December 2023
##
## R-code and data usage license: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)
##
## Copyright (c) Theresa Kraft, 2023
## Email: theresa.kraft@tu-dresden.de
##
## ---------------------------
##
## Notes: - Please set your directories before you start. 
##        - All required R-functions can be found in the R-script HelperFunctions.R.
##
## ---------------------------


################ 
## DIRECTORIES 
################ 

# workingDirectory
workingDirectory <- "/home/kraftth/ExpressionDirectory/"

# This directory is where data for figures are collected and where lateron, figures are created separately 
figureDirectory <- paste0(workingDirectory, "Figures/")

# This directory is where the expression data is stored
# please download data from TODO
dataDirectory <- paste0(workingDirectory, "Data/")

# This directory is where all annotation files are stored
annotationDirectory <- paste0(workingDirectory, "Annotations/")

# This directory is where you can find the ARHMM.jar file 
hmmDirectory <- paste0(workingDirectory, "Programs/")


setwd(workingDirectory)


################ 
## LIBRARIES 
################ 

library(stringr)
library(pvclust)
library(pbapply)
library(plyr)
library(clusterProfiler)
library("org.Hs.eg.db")

################ 
## START 
################ 


# read in expression data file 
expressionSamples <- read.csv(file=paste0(dataDirectory,
              "MelBrainSys_geneExpression_1822_metastases.tsv"),
              sep = "\t", check.names = T) 

################ 
## FIGURE 1
# stability dendrogram
################ 

# create data for dendrogram 
expressionValues <- expressionSamples[,9:45]


# use pvclust with 10000 bootstraping repetitions for stability analysis
# parallel = T or define parallel = as.integer(numCores) if parallel computing is desired
# took approx. 14 minutes in a parallel mode with 100 cores

t1 <- Sys.time()
valuesManhD2PV <- pvclust(expressionValues, method.hclust = "ward.D2", method.dist = "manhattan", parallel = as.integer(100), nboot = 10000)

print(difftime(Sys.time(), t1) )
timePassed <- difftime(Sys.time(), t1)

outputFileFig1 <- paste0(figureDirectory, "Figure1-Stability-Dendrogram/valuesManhD2PV.RData" ) 
save(valuesManhD2PV, file = outputFileFig1)


################ 
## Metastases pairs 
## creation 
################ 

# create metastases pairs: Intracranial expression value minus extacranial expression value 
expresssionPairs <- data.frame("ID" = expressionSamples$Ensembl.Gene.ID, "Gene" = expressionSamples$Associated.Gene.Name, "Chr" = expressionSamples$Chr, 
                               "Start" = expressionSamples$Start, "End" = expressionSamples$End ,  "Strand" = expressionSamples$Strand,  
                               "P101_BLiv" = expressionSamples$P101_B - expressionSamples$P101_Liv, "P106_BLym" = expressionSamples$P106_B - expressionSamples$P106_Lym, 
                               "P107_BLun" = expressionSamples$P107_B - expressionSamples$P107_Lun, "P108_BLym" = expressionSamples$P108_B - expressionSamples$P108_Lym, 
                               "P111_BLym" = expressionSamples$P111_B - expressionSamples$P111_Lym, "P13_BLym" = expressionSamples$P13_B - expressionSamples$P13_Lym, 
                               "P16_BLun" = expressionSamples$P16_B - expressionSamples$P16_Lun, "P18_BLun-1" = expressionSamples$P18_B - expressionSamples$P18_Lun_1, 
                               "P18_BLun-2" = expressionSamples$P18_B - expressionSamples$P18_Lun_2, "P3_BLun" = expressionSamples$P39_B - expressionSamples$P3_Lun, 
                               "P39_BLun" = expressionSamples$P39_B - expressionSamples$P39_Lun, "P4_BSki-1" = expressionSamples$P4_B - expressionSamples$P4_Ski_1, 
                               "P4_BSki-2" = expressionSamples$P4_B - expressionSamples$P4_Ski_2, "P42_BLym-1" = expressionSamples$P42_B - expressionSamples$P42_Lym_1,
                               "P42_BLym-2" = expressionSamples$P42_B - expressionSamples$P42_Lym_2, "P74_BLym" = expressionSamples$P74_B - expressionSamples$P74_Lym, 
                               "P77_Blym" = expressionSamples$P77_B - expressionSamples$P77_Lym, "P78_BSmi" = expressionSamples$P78_B - expressionSamples$P78_Smi, 
                               "P8_BSof-1" = expressionSamples$P8_B - expressionSamples$P8_Sof_1, "P8_BSof-2" = expressionSamples$P8_B - expressionSamples$P8_Sof_2, 
                               "P8_BSof-3" = expressionSamples$P8_B - expressionSamples$P8_Sof_3, check.names = F )

expresssionPairs <- expresssionPairs[order(expresssionPairs$Chr,expresssionPairs$Start),]


################ 
## FIGURE 2 
# autocorrelation
################ 

allAcfBiol <- create_acf_df_weighted_biol(expresssionPairs)

set.seed(1)
rows <- sample(nrow(expresssionPairs))
permutatedDF <- expresssionPairs[rows, ]

allAcfPermutated <- create_acf_df_weighted(permutatedDF)


for(i in 2:1000){
  set.seed(i)
  rows <- sample(nrow(expresssionPairs))
  permutatedDF <- expresssionPairs[rows, ]
  
  allAcfPermutated2 <- create_acf_df_weighted(permutatedDF)
  
  allAcfPermutated <- cbind(allAcfPermutated, allAcfPermutated2 )
}

allAcfPermutatedPlot <- data.frame("lag" = 0:20, "acfMean" = apply(allAcfPermutated, MARGIN =1, FUN = mean), "acfMedian" =apply(allAcfPermutated, MARGIN =1, FUN = median), 
                                "acfSd" = apply(allAcfPermutated, MARGIN =1, FUN =sd),"acfBiolMean" = allAcfBiol$acfMean, "acfBiolMedian" = allAcfBiol$acfMedian, 
                                "sdBiol"  = allAcfBiol$acfSd)

allAcfPermutatedPlot <- allAcfPermutatedPlot[-1,]

outputFileFig2 <- paste0(figureDirectory, "Figure2-Autocorrelation/allAcfPermutatedPlot.RData" ) 
save(allAcfPermutatedPlot, file = outputFileFig2)


################ 
## HMM training 
################ 


# create training data in the desired format for ARHMM
hmmTrainingData <- data.frame("GeneID" = expresssionPairs$ID, "chr" = expresssionPairs$Chr, 
                              "start" = expresssionPairs$Start, expresssionPairs[,7:27],check.names = F)


outputFileHMMTrain <- paste0(workingDirectory, "Expression_TrainingData.txt")
write.table( hmmTrainingData, file = outputFileHMMTrain, row.names = FALSE,
             col.names = TRUE, dec = ".", sep = "\t", quote = FALSE )

# bash script for HMM training using optimal parameters
# will create the StatePosteriorDecoding file that contains the most probable state for each CpG in each metastases pair
# ARHMM_Trainer.jar from Seifert et. al, 2014 (PMID: 24955771) 
bash <- paste0("java -jar ",hmmDirectory,"ARHMM_Trainer.jar -hmmOrder 1 -arOrder 0 -initialStateDist 0.1 0.8 0.1 -means -3 0 3 -sds 0.5 1 0.5 -scaleMeans 2500 1000 2500 -shapeSds 5000 10 5000 -scaleSds 1E-4 1E-4 1E-4 -statePosterior -modelBasics -model model1 -dataSet Expression_TrainingData.txt")
system(bash)


expressionClassification <- read.delim(paste0(workingDirectory, 
                                               "model1_Expression_TrainingData.txt_StatePosteriorDecoding.txt"),check.names = F)



################ 
## Create big table including HMM predictions  
################ 

# merge expression classification and expression values to one dataframe
expressionPairsRun <- merge(expressionClassification, expresssionPairs, by.x = "Genes",  by.y = "ID" )

colnames(expressionPairsRun)[1] <- "ID"



################ 
## annotate table  
################ 


pathwayPath <- paste0(annotationDirectory, "Relevant_Pathways_KEGG_103-0_ENSEMBL_106_GenomeV_GrCh38-p13_new.txt")
pathways <- read.table( file=pathwayPath, header=T, check.names = FALSE, sep = "\t")




pathwayArraySignal <- c()
for(i in 1:nrow(pathways)){
  gene <- as.character(pathways[i,1])
  path <- ""
  for(j in 8:36){
    if(as.character(pathways[i,j]) == "1"){
      path <- paste0(path, colnames(pathways)[j],";")
    }
  }
  pathwayArraySignal <- rbind(pathwayArraySignal, c(gene, path))
}
pathwayArraySignal <- as.data.frame(pathwayArraySignal)
colnames(pathwayArraySignal) <- c("EnsID", "pathways")

expressionPairs.annotated <- merge(expressionPairsRun, pathwayArraySignal, by.x = "ID", by.y = "EnsID")





pathwayArrayImmune <- c()
for(i in 1:nrow(pathways)){
  gene <- as.character(pathways[i,1])
  path <- ""
  for(j in 37:49){
    if(as.character(pathways[i,j]) == "1"){
      path <- paste0(path, colnames(pathways)[j],";")
    }
  }
  pathwayArrayImmune <- rbind(pathwayArrayImmune, c(gene, path))
}
pathwayArrayImmune <- as.data.frame(pathwayArrayImmune)
colnames(pathwayArrayImmune) <- c("EnsID", "immunePathways")

expressionPairs.annotated <- merge(expressionPairs.annotated, pathwayArrayImmune, by.x = "ID", by.y = "EnsID")





# change colnames, metastases pairs *_p = expression prediction
colnames(expressionPairs.annotated) <-  c("ID","P101_BLiv_p","P106_BLym_p","P107_BLun_p","P108_BLym_p","P111_BLym_p","P13_BLym_p","P16_BLun_p","P18_BLun-1_p","P18_BLun-2_p",
                                          "P3_BLun_p","P39_BLun_p","P4_BSki-1_p","P4_BSki-2_p","P42_BLym-1_p","P42_BLym-2_p","P74_BLym_p","P77_BLym_p",
                                          "P78_BSmi_p","P8_BSof-1_p","P8_BSof-2_p","P8_BSof-3_p","Gene","Chr","Start","End","Strand","P101_BLiv","P106_BLym","P107_BLun",
                                          "P108_BLym","P111_BLym","P13_BLym","P16_BLun","P18_BLun-1","P18_BLun-2","P3_BLun","P39_BLun","P4_BSki-1","P4_BSki-2",
                                          "P42_BLym-1","P42_BLym-2","P74_BLym","P77_BLym","P78_BSmi","P8_BSof-1","P8_BSof-2","P8_BSof-3",
                                          "pathways","immunePathways")




PatientID <- c()
Expression <- c()
value <- c()
overviewDF <- NULL
for(i in 28:48){
  sampleValues <- as.character(colnames(expressionPairs.annotated)[i])
  samplePredictions <- paste0(sampleValues,"_p")
  
  
  highExpr <- sum(expressionPairs.annotated[,samplePredictions] == "+")
  lowExpr <- sum(expressionPairs.annotated[,samplePredictions] == "-")
  
  PatientID <- c(PatientID, sampleValues)
  Expression <- c(Expression, "Decreased Expression")
  value <- c(value, lowExpr)
  
  PatientID <- c(PatientID, sampleValues )
  Expression <- c(Expression, "Increased Expression")
  value <- c(value, highExpr)
}


ranking <- c()
for(i in 1:length(PatientID)){
  if(str_detect(PatientID[i], "Lu")){
    ranking <- c(ranking,5)
  }
  if(str_detect(PatientID[i], "Lym")){
    ranking <- c(ranking,6)
  }
  if(str_detect(PatientID[i], "Ski")){
    ranking <- c(ranking,3)
  }
  if(str_detect(PatientID[i], "Liv")){
    ranking <- c(ranking,2)
  }
  if(str_detect(PatientID[i], "Sof")){
    ranking <- c(ranking,4)
  }
  if(str_detect(PatientID[i], "Smi")){
    ranking <- c(ranking,1)
  }
}

barplotDF <- data.frame(PatientID, Expression, value, ranking)

for(i in 1:length(barplotDF[,1])){
  if(as.character(barplotDF[i,2]) == "Decreased Expression"){
    barplotDF[i,3] <- - barplotDF[i,3]
  }
}
barplotOverviewDF <-  barplotDF


outputFileFig3 <- paste0(figureDirectory, "Figure3-Genomic-Regions/barplotOverviewDF.RData" ) 
save(barplotOverviewDF, file = outputFileFig3)


mean(barplotOverviewDF$value[barplotOverviewDF$Expression == "Decreased Expression"])
sd(barplotOverviewDF$value[barplotOverviewDF$Expression == "Decreased Expression"])

mean(barplotOverviewDF$value[barplotOverviewDF$Expression == "Increased Expression"])
sd(barplotOverviewDF$value[barplotOverviewDF$Expression == "Increased Expression"])


pparDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "ppar"),]
mapkDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "mapk"),]
erbbDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "erbb"),]
cytokineDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "cytokine"),]
cell_cycleDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "cell_cycle"),]
p53DF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "p53"),]
mtorDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "mtor"),]
pi3k_aktDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "pi3k_akt"),]
apoptosisDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "apoptosis"),]
wntDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "wnt"),]
tgfbDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "tgfb"),]
vegfDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "vegf"),]
focal_adhDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "focal_adh"),]
ECMDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "ECM"),]
adh_junDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "adh_jun"),]
jak_statDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "jak_stat"),]
notchDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "notch"),]
hedgehogDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "hedgehog"),]
replicationDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "replication"),]
BERDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "BER"),]
NERDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "NER"),]
HRDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "HR"),]
NHEJDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "NHEJ"),]
mismatchDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "mismatch"),]
telomereDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "telomere"),]
calciumDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "calcium"),]
cAMPDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "cAMP"),]
HIF1DF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "HIF1"),]
estrogenDF <-  expressionPairs.annotated[str_detect(expressionPairs.annotated$pathways, "estrogen"),]

listofDFs <- list(pparDF, mapkDF , erbbDF , cytokineDF , cell_cycleDF , p53DF , mtorDF , pi3k_aktDF , apoptosisDF , wntDF , tgfbDF , vegfDF , focal_adhDF ,
                  ECMDF , adh_junDF , jak_statDF , notchDF , hedgehogDF , replicationDF , BERDF , NERDF , HRDF , NHEJDF , 
                  mismatchDF , telomereDF,calciumDF,cAMPDF,HIF1DF,estrogenDF)
listofAnnotations <- list("ppar","mapk","erbb","cytokine","cell_cycle","p53","mtor","pi3k_akt","apoptosis","wnt","tgfb",	
                          "vegf","focal_adh","ECM","adh_jun","jak_stat","notch","hedgehog","replication","BER","NER",
                          "HR","NHEJ","mismatch","telomere","calcium","cAMP","HIF1","estrogen")


pathwayGenePlotDf <- create_annotation_analysis_df(listOfAnnotations = listofAnnotations, listOfDFs = listofDFs)


outputFileFig3 <- paste0(figureDirectory, "Figure3-Genomic-Regions/pathwayGenePlotDf.RData" ) 
save(pathwayGenePlotDf, file = outputFileFig3)







MelanomaDF <-   expressionPairs.annotated[str_detect(expressionPairs.annotated$immunePathways, "melanoma"),]
PDL1exprDF <-   expressionPairs.annotated[str_detect(expressionPairs.annotated$immunePathways, "pd1"),]
NKillerDF <-   expressionPairs.annotated[str_detect(expressionPairs.annotated$immunePathways, "Nkcyt"),]
RIGIlikeDF <-   expressionPairs.annotated[str_detect(expressionPairs.annotated$immunePathways, "RIGI"),]
NODlikeDF <-   expressionPairs.annotated[str_detect(expressionPairs.annotated$immunePathways, "NOD"),]
TolllikeDF <-   expressionPairs.annotated[str_detect(expressionPairs.annotated$immunePathways, "Toll"),]
TcellReceptorDF <-   expressionPairs.annotated[str_detect(expressionPairs.annotated$immunePathways, "TcellR"),]
AntigenDF <-   expressionPairs.annotated[str_detect(expressionPairs.annotated$immunePathways, "antigen"),]
BcellReceptorDF <-   expressionPairs.annotated[str_detect(expressionPairs.annotated$immunePathways, "BcellR"),]
Th1Th2DF <-   expressionPairs.annotated[str_detect(expressionPairs.annotated$immunePathways, "Th1Th2"),]
Th17DF <-   expressionPairs.annotated[str_detect(expressionPairs.annotated$immunePathways, "Th17"),]
Il17DF <-   expressionPairs.annotated[str_detect(expressionPairs.annotated$immunePathways, "IL17"),]
LeukocyteDF <-   expressionPairs.annotated[str_detect(expressionPairs.annotated$immunePathways, "ltm"),]

listofDFs <- list(MelanomaDF, PDL1exprDF, NKillerDF, RIGIlikeDF, NODlikeDF, TolllikeDF, TcellReceptorDF, AntigenDF, BcellReceptorDF, Th1Th2DF,
                  Th17DF, Il17DF, LeukocyteDF)
listofAnnotations <- list("Melanoma","PD-L1 expression and PD-1 checkpoint pathway","Natural killer cell mediated cytotoxicity",
                          "RIG-I-like receptor","NOD-like receptor","Toll-like receptor","T cell receptor","Antigen processing and presentation",
                          "B cell receptor","Th1 and Th2","Th17",	
                          "IL-17","Leukocyte transendothelial migration")


immunePlotDf <- create_annotation_analysis_df(listOfAnnotations = listofAnnotations, listOfDFs = listofDFs)



outputFileFig3 <- paste0(figureDirectory, "Figure3-Genomic-Regions/immunePlotDf.RData" ) 
save(immunePlotDf, file = outputFileFig3)




################ 
## similarities of samples from the same patient  
################ 

expressionValuesPairs <- expressionPairs.annotated[,28:48]


valuesManhD2PVPairs <- pvclust(expressionValuesPairs, method.hclust = "ward.D2", method.dist = "manhattan", parallel = as.integer(100), nboot = 10000)


save(valuesManhD2PVPairs, file = paste0(figureDirectory, "Figure4-histologicalRegions/valuesManhD2PVPairs.RData" ) )




diffsBetweenPairsTable <- countDiffsBetweenPairs("P18_BLun-1","P18_BLun-2")
diffsBetweenPairsTable <- rbind(diffsBetweenPairsTable, countDiffsBetweenPairs("P4_BSki-1","P4_BSki-2"))
diffsBetweenPairsTable <- rbind(diffsBetweenPairsTable, countDiffsBetweenPairs("P42_BLym-1","P42_BLym-2"))
diffsBetweenPairsTable <- rbind(diffsBetweenPairsTable, countDiffsBetweenPairs("P8_BSof-1","P8_BSof-2"))
diffsBetweenPairsTable <- rbind(diffsBetweenPairsTable, countDiffsBetweenPairs("P8_BSof-1","P8_BSof-3"))
diffsBetweenPairsTable <- rbind(diffsBetweenPairsTable, countDiffsBetweenPairs("P8_BSof-2","P8_BSof-3"))

diffsBetweenPairsTable <- as.data.frame(diffsBetweenPairsTable)

colnames(diffsBetweenPairsTable) <- c("comparison", "++", "--", "==", "+=", "=+", "-=", "=-","+-","-+")

allPatients <- colnames(expressionPairs.annotated)[28:48]

diffsBetweenSamplesTable <- c()

for(i in 1:length(allPatients)){
  for(j in 1:length(allPatients)){
    if(i<j){
      diffsBetweenSamplesTable <- rbind(diffsBetweenSamplesTable,countDiffsBetweenPairs(allPatients[i],allPatients[j]))
    }
  }
}

colnames(diffsBetweenSamplesTable) <- c("comparison", "++", "--", "==", "+=", "=+", "-=", "=-","+-","-+")


diffsBetweenSamplesTable <- as.data.frame(diffsBetweenSamplesTable)
diffsBetweenSamplesTable$`++` <- as.numeric(as.character(diffsBetweenSamplesTable$`++`))
diffsBetweenSamplesTable$`--` <- as.numeric(as.character(diffsBetweenSamplesTable$`--`))
diffsBetweenSamplesTable$`==` <- as.numeric(as.character(diffsBetweenSamplesTable$`==`))
diffsBetweenSamplesTable$`+=` <- as.numeric(as.character(diffsBetweenSamplesTable$`+=`))
diffsBetweenSamplesTable$`=+` <- as.numeric(as.character(diffsBetweenSamplesTable$`=+`))
diffsBetweenSamplesTable$`-=` <- as.numeric(as.character(diffsBetweenSamplesTable$`-=`))
diffsBetweenSamplesTable$`=-` <- as.numeric(as.character(diffsBetweenSamplesTable$`=-`))
diffsBetweenSamplesTable$`+-` <- as.numeric(as.character(diffsBetweenSamplesTable$`+-`))
diffsBetweenSamplesTable$`-+` <- as.numeric(as.character(diffsBetweenSamplesTable$`-+`))

diffsBetweenSamplesTable$OverlapPerc <- (( diffsBetweenSamplesTable$`++` + diffsBetweenSamplesTable$`--` + diffsBetweenSamplesTable$`==` ) / sum(diffsBetweenSamplesTable[1,2:10]) ) *100

diffsBetweenSamplesTable$PairComparison <- "different patients"
diffsBetweenSamplesTable$PairComparison[diffsBetweenSamplesTable$comparison %in% diffsBetweenPairsTable$comparison ] <- "same patient"

diffsBetweenSamplesTable$PairComparison <- as.factor(diffsBetweenSamplesTable$PairComparison)


save(diffsBetweenSamplesTable, file = paste0(figureDirectory, "Figure4-histologicalRegions/diffsBetweenSamplesTable.RData" ) )


################ 
## mean data from multiple patients  
################ 

runPatData <- data.frame("ID" = expressionPairs.annotated$ID, "P101_BLiv" = expressionPairs.annotated$P101_BLiv_p, 
                         "P106_BLym" = expressionPairs.annotated$P106_BLym_p, 
                         "P107_BLun" = expressionPairs.annotated$P107_BLun_p, 
                         "P108_BLym" = expressionPairs.annotated$P108_BLym_p, 
                         "P111_BLym"= expressionPairs.annotated$P111_BLym_p, 
                         "P13_BLym" = expressionPairs.annotated$P13_BLym_p, 
                         "P16_BLun" = expressionPairs.annotated$P16_BLun_p ,
                         "P18_BLun" =  create_mean_from_runDataCols(expressionPairs.annotated$`P18_BLun-1_p`, expressionPairs.annotated$`P18_BLun-2_p`), 
                         "P3_BLun" = expressionPairs.annotated$P3_BLun_p, 
                         "P39_BLun" = expressionPairs.annotated$P39_BLun_p, 
                         "P4_BSki" = create_mean_from_runDataCols(expressionPairs.annotated$`P4_BSki-1_p`, expressionPairs.annotated$`P4_BSki-2_p`),   
                         "P42_BLym" = create_mean_from_runDataCols(expressionPairs.annotated$`P42_BLym-1_p`, expressionPairs.annotated$`P42_BLym-2_p`), 
                         "P74_BLym" = expressionPairs.annotated$P74_BLym_p, 
                         "P77_BLym" = expressionPairs.annotated$P77_BLym_p, 
                         "P78_BSmi" = expressionPairs.annotated$P78_BSmi_p, 
                         "P8_BSof" = create_mean_from_runDataCols(expressionPairs.annotated$`P8_BSof-1_p`, expressionPairs.annotated$`P8_BSof-2_p`, expressionPairs.annotated$`P8_BSof-3_p`)   )


################ 
## Rank genes according to their number of diff. expression in patients   
################ 


expressionPairsCount <- data.frame( "ID" = runPatData$ID, "countNeg" = apply(runPatData, MARGIN = 1, function(x) sum(x == "-")), "countPos" = apply(runPatData, MARGIN = 1, function(x) sum(x == "+")) )
expressionPairs.annotated_count <- merge(expressionPairsCount, expressionPairs.annotated, by.x = "ID",  by.y = "ID" )



numCandGens <- get_number_of_candidates_on_cutoff(expressionPairs.annotated_count)

save(numCandGens, file = paste0(figureDirectory, "Figure5-GOenrichment/numCandGens.RData" ) )



################ 
## Create GO analysis for genes altered in >= 50% of all patients 
################ 

candGenesPlus <- expressionPairs.annotated_count$ID[expressionPairs.annotated_count$countPos >= 8]
candGenesMinus <- expressionPairs.annotated_count$ID[expressionPairs.annotated_count$countNeg >= 8]


egoPlus <- enrichGO(gene          = candGenesPlus,
                    universe      = expressionPairs.annotated_count$ID,
                    keyType = "ENSEMBL",
                    OrgDb         = org.Hs.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "fdr",
                    pvalueCutoff  = 0.05,
                    readable      = TRUE)

egoMinus <- enrichGO(gene          = candGenesMinus,
                     universe      = expressionPairs.annotated_count$ID,
                     keyType = "ENSEMBL",
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "fdr",
                     pvalueCutoff  = 0.05,
                     readable      = TRUE)


resultPlus <- egoPlus@result[order(egoPlus@result$qvalue),]
resultMinus <- egoMinus@result[order(egoMinus@result$qvalue),]

resultPlus$Expression.direction <- "Increased"
resultMinus$Expression.direction <- "Decreased"
resultGO <- rbind(resultMinus,resultPlus)


save(resultGO, file = paste0(figureDirectory, "Figure5-GOenrichment/resultGO.RData" ) )



################ 
## Gene candidates: differentially expressed in 11 or more patients  
################ 



geneCandidates <- expressionPairs.annotated_count[expressionPairs.annotated_count$countPos >= 11 | expressionPairs.annotated_count$countNeg >= 11,c(25,30:50)]
rownames(geneCandidates) <- geneCandidates$Gene

geneCandidates$Gene <- NULL

save(geneCandidates, file = paste0(figureDirectory, "Figure6-candHeatmap/geneCandidates.RData" ) )


# overlap of gene candidates with other studies 

Fischer_candidate_table <- read.csv(paste0(annotationDirectory, "Fischer2019_differential_expression_intra_vs_extra_melanomaMetastases.csv"))

Biermann2022__all_table <- read.csv(paste0(annotationDirectory, "Biermann2022_singleCell_all_genes.csv") ) 
Biermann2022__cand_table <- Biermann2022__all_table[Biermann2022__all_table$p_val_adj < 0.05,]

Chen2014_cand_table <- read.csv(paste0(annotationDirectory, "/Chen2014_differential_expression_intra_vs_extra_melanomaMetastases.csv"))


candGenes_biolFunction_summary <- read.csv(paste0(annotationDirectory, "candGenes_biolFunction_summary.csv"))

candGenes_biolFunction_summary$Fischer <- "not listed"
candGenes_biolFunction_summary$Fischer[candGenes_biolFunction_summary$Gene %in% Fischer_candidate_table$ID] <- "listed"
candGenes_biolFunction_summary$Fischer <- as.character(candGenes_biolFunction_summary$Fischer)

candGenes_biolFunction_summary$Chen <- "not listed"
candGenes_biolFunction_summary$Chen[candGenes_biolFunction_summary$Gene %in% Chen2014_cand_table$Gene.Symbol] <- "listed"
candGenes_biolFunction_summary$Chen <- as.character(candGenes_biolFunction_summary$Chen)


candGenes_biolFunction_summary$Biermann <- "not listed"
candGenes_biolFunction_summary$Biermann[candGenes_biolFunction_summary$Gene %in% Biermann2022__cand_table$gene] <- "listed"
candGenes_biolFunction_summary$Biermann <- as.character(candGenes_biolFunction_summary$Biermann)

for(g in 1:length(candGenes_biolFunction_summary$Gene)){
  gene <- as.character(candGenes_biolFunction_summary$Gene[g])
  if(candGenes_biolFunction_summary$Fischer[g] == "listed"){
    if(Fischer_candidate_table$logFC[Fischer_candidate_table$ID == gene] > 0){
      candGenes_biolFunction_summary$Fischer[g] <- "increased"
    } else{
      candGenes_biolFunction_summary$Fischer[g] <- "decreased"
    }
  }
  
  if(candGenes_biolFunction_summary$Chen[g] == "listed"){
    if(Chen2014_cand_table$Linear.ratio.BM.EM.[Chen2014_cand_table$Gene.Symbol == gene][1] > 1){
      candGenes_biolFunction_summary$Chen[g] <- "increased"
    } else{
      candGenes_biolFunction_summary$Chen[g] <- "decreased"
    }
  }
  
  
  if(candGenes_biolFunction_summary$Biermann[g] == "listed"){
    if(Biermann2022__cand_table$avg_log2FC[Biermann2022__cand_table$gene == gene] > 0){
      candGenes_biolFunction_summary$Biermann[g] <- "increased"
    } else{
      candGenes_biolFunction_summary$Biermann[g] <- "decreased"
    }
  }
  
}

write.csv(candGenes_biolFunction_summary, file =paste0(figureDirectory,"Figure6-candHeatmap/candGenes_biolFunction_summary_validation.csv"), row.names = F)


### overlap perc and pval 

# Biermann 

1- ( sum(candGenes_biolFunction_summary$Biermann == "not listed") / length (candGenes_biolFunction_summary$Biermann) )

inBoth <- length(candGenes_biolFunction_summary$Gene) - sum(candGenes_biolFunction_summary$Biermann == "not listed")
BiermannOnly <- length(Biermann2022__cand_table$gene[Biermann2022__cand_table$gene %in% expressionPairs.annotated_count$Gene]) - inBoth
OursOnly <- length(candGenes_biolFunction_summary$Gene) - inBoth
inNone <- length(expressionPairs.annotated_count$Gene) - inBoth - BiermannOnly - OursOnly 

fisher.test(matrix(c(inBoth,BiermannOnly,OursOnly,inNone),ncol = 2))

1- sum(candGenes_biolFunction_summary$Chen == "not listed") / length (candGenes_biolFunction_summary$Chen)


inBoth <- length(candGenes_biolFunction_summary$Gene) - sum(candGenes_biolFunction_summary$Chen == "not listed")
ChenOnly <- length(Chen2014_cand_table$Gene.Symbol[Chen2014_cand_table$Gene.Symbol %in% expressionPairs.annotated_count$Gene]) - inBoth
OursOnly <- length(candGenes_biolFunction_summary$Gene) - inBoth
inNone <- length(expressionPairs.annotated_count$Gene) - inBoth - ChenOnly - OursOnly 

fisher.test(matrix(c(inBoth,ChenOnly,OursOnly,inNone),ncol = 2))



1- sum(candGenes_biolFunction_summary$Fischer == "not listed") / length (candGenes_biolFunction_summary$Fischer)


inBoth <- length(candGenes_biolFunction_summary$Gene) - sum(candGenes_biolFunction_summary$Fischer == "not listed")
FischerOnly <- length(Fischer_candidate_table$ID[Fischer_candidate_table$ID %in% expressionPairs.annotated_count$Gene]) - inBoth
OursOnly <- length(candGenes_biolFunction_summary$Gene) - inBoth
inNone <- length(expressionPairs.annotated_count$Gene) - inBoth - FischerOnly - OursOnly 

fisher.test(matrix(c(inBoth,FischerOnly,OursOnly,inNone),ncol = 2))







################ 
## Association of gene candidates and survival using TCGA data   
################ 


# read in TCGA public expression and patient data
TCGA_expression_all <- read.csv(file=paste0(dataDirectory,
                                            "tcga-raw-expression-annotated-usedSamples_FiltLowExpress_Normalized_cyclic_logCPM.txt"),
                                sep = "\t", check.names = F) 
TCGA_patientInfos <- read.csv(file=paste0(dataDirectory,
                                          "Patient_information.csv"),check.names = F)


candidateGenesTable <- rownames(geneCandidates)

candidateGenes_TCGAdata <- TCGA_expression_all[TCGA_expression_all$GeneID %in% candidateGenesTable,]





# create a data frame in a format for Kaplan Meier analysis

TCGA_survival_df <- TCGA_patientInfos[,c(1,105,107)]
TCGA_survival_df[,2] <- as.numeric(as.character(TCGA_survival_df[,2]))
TCGA_survival_df[,3] <-   as.numeric(as.character(TCGA_survival_df[,3]))
TCGA_survival_df <- TCGA_survival_df[complete.cases(TCGA_survival_df),]
TCGA_survival_df$Group <- "no Group"


candidateGenes_TCGAdata.Survival <- candidateGenes_TCGAdata[,colnames(candidateGenes_TCGAdata) %in% TCGA_survival_df$Name | colnames(candidateGenes_TCGAdata) == "GeneID"]



supplTable <- data.frame(stringsAsFactors = F)


pvalueDF <- c()


# create survival analysis for each single gene and track pvalue

fitList <- list()
for(g in 1:length(candidateGenes_TCGAdata.Survival$GeneID)){
  gene <- as.character(candidateGenes_TCGAdata.Survival$GeneID[g])
  
  geneDF <- candidateGenes_TCGAdata.Survival[candidateGenes_TCGAdata.Survival$GeneID == gene,-1]
  FirstQuartile <- quantile(as.numeric(geneDF[1,]),0.25)
  ThirdQuartile <- quantile(as.numeric(geneDF[1,]),0.75)
  
  
  highExpressionID <- colnames(geneDF)[geneDF[1,] > ThirdQuartile]
  lowExpressionID <- colnames(geneDF)[geneDF[1,] < FirstQuartile]
  
  
  for(h in 1:length(highExpressionID)){
    supplTable <- rbind(supplTable, data.frame("Gene" = gene,"Patient" = highExpressionID[h], "days_to_death" =  TCGA_survival_df[TCGA_survival_df$Name ==highExpressionID[h],2],
                                               "status" = TCGA_survival_df[TCGA_survival_df$Name ==highExpressionID[h],3], 
                                               "group" = "high Expression"))
  }
  
  
  for(l in 1:length(lowExpressionID)){
    supplTable <- rbind(supplTable, data.frame("Gene" = gene, "Patient" = lowExpressionID[l],  "days_to_death" =  TCGA_survival_df[TCGA_survival_df$Name ==lowExpressionID[l],2],
                                               "status" = TCGA_survival_df[TCGA_survival_df$Name ==lowExpressionID[l],3], 
                                               "group" = "low Expression"))
  }
  
  survAnalysisList <- perform_survival_analysis_multipleTesting(lowExpressionID,highExpressionID,gene)
  
  
  fitList[[g]] <- survAnalysisList[[1]]
  pvalueDF <- rbind(pvalueDF, data.frame("Gene" = gene, "pval" = survAnalysisList[[2]]))
}

pvalueDF$adj.pval <- p.adjust(pvalueDF$pval, method = "fdr")

write.table(supplTable, file =paste0(figureDirectory,"Figure7-survival/TableS8_Survival_Patient_Information.tsv"), row.names = F,sep = "\t")


# create plot for each significant gene
legend <- NULL
plotList <- list()
counter <- 1
for(f in 1:length(fitList)){
  fit <- fitList[[f]]
  gene <- as.character(candidateGenes_TCGAdata.Survival$GeneID[f])
  
  label <-  paste0("q = ", round(pvalueDF$adj.pval[pvalueDF$Gene == gene],3))
  plot <- ggsurvplot(fit, title = gene, xlim = c(0,11000))$plot + theme_gray() + geom_label(x = 2174, y = 0.1, label = label, inherit.aes = F) + 
    scale_x_continuous(breaks = c(0,2500,5000,7500,10000)) + theme(legend.position = "none")
  
  
  if(pvalueDF$adj.pval[pvalueDF$Gene == gene] < 0.1){
    plotList[[counter]] <- plot
    counter <- counter + 1
  }
  legend <- get_legend(ggsurvplot(fit, title = gene)$plot + theme_gray() + geom_label(x = 2174, y = 0.1, label = label, inherit.aes = F))
}



plotList[[12]] <- legend 

save(plotList, file =paste0(figureDirectory,"Figure7-survival/plotListSurvival.RData"))




