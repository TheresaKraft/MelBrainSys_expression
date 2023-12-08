## ---------------------------
##
## This script was created with R version 4.2.1 (2022-06-23)
##
##
## Title: Personalized identification and characterization of genome-wide gene expression differences 
##        between patient-matched intracranial and extracranial melanoma metastasis pairs
##
## Purpose of script: Helper functions for script MainScript.R
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
## Notes: This is a collections of all helper functions used in the MainScript.R
##   
##
## ---------------------------



#############################################################################################################
##
### function calculates autocorrelation of any dataframe in a specific order
### measurementDF: Dataframe containing measurement values in specific order
##
############################################################################################################



create_acf_df_weighted <- function(rnaDataF ){
  sampleAcfAll <- data.frame("lag" = 0:20)
  for(sampleId in 7:27){
    rnaDataM = as.matrix(rnaDataF[,c(sampleId)])
    
    chrIdx = sapply(unique(rnaDataF$Chr), function(chr) which(rnaDataF$Chr == chr))
    
    neededChr = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)
    
    
    acf <- acf(rnaDataM[chrIdx[[1]],1], lag.max = 20,plot = F)
    chrLength <- length(rnaDataM[chrIdx[[1]],1])
    acf$acf <- chrLength * acf$acf
    for(i in 2:length(neededChr)){
      chrLength <- length(rnaDataM[chrIdx[[i]],1])
      acf$acf <- acf$acf + (chrLength * acf(rnaDataM[chrIdx[[i]],1], lag.max = 20,plot = F)$acf)
    }
    acf$acf <- acf$acf / length(rnaDataM[,1])
    
    sampleAcf <- data.frame("lag" = 0:20, "acf" = acf$acf)
    
    sampleAcfAll <- merge(sampleAcfAll, sampleAcf, by.x = "lag", by.y = "lag")
  }
  
  sampleAcfAll <- data.frame("acf" = apply(sampleAcfAll[,-1], 1, mean) )
  return(sampleAcfAll)
}



create_acf_df_weighted_biol <- function(rnaDataF ){
  sampleAcfAll <- data.frame("lag" = 0:20)
  for(sampleId in 7:27){
    rnaDataM = as.matrix(rnaDataF[,c(sampleId)])
    
    chrIdx = sapply(unique(rnaDataF$Chr), function(chr) which(rnaDataF$Chr == chr))
    
    neededChr = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)
    
    
    acf <- acf(rnaDataM[chrIdx[[1]],1], lag.max = 20,plot = F)
    chrLength <- length(rnaDataM[chrIdx[[1]],1])
    acf$acf <- chrLength * acf$acf
    for(i in 2:length(neededChr)){
      chrLength <- length(rnaDataM[chrIdx[[i]],1])
      acf$acf <- acf$acf + (chrLength * acf(rnaDataM[chrIdx[[i]],1], lag.max = 20,plot = F)$acf)
    }
    acf$acf <- acf$acf / length(rnaDataM[,1])
    
    sampleAcf <- data.frame("lag" = 0:20, "acf" = acf$acf)
    
    sampleAcfAll <- merge(sampleAcfAll, sampleAcf, by.x = "lag", by.y = "lag")
  }
  
  sampleAcfAll <- data.frame("acfMean" = apply(sampleAcfAll[,-1], 1, mean), "acfMedian" = apply(sampleAcfAll[,-1], 1, median), "acfSd" = apply(sampleAcfAll[,-1], 1, sd) )
  return(sampleAcfAll)
}




#############################################################################################################
##
### function that extracts pathways
###
##
############################################################################################################



create_annotation_analysis_df <- function(listOfAnnotations, listOfDFs){
  
  if(length(listOfAnnotations) != length(listofDFs)){
    stop("Length of lists are not equal! Returning without calculations.")
  }
  ##
  # calculate the percentage of CpGs in + and - dependent on their rnaDataRunDf 
  ##
  
  # reset all values
  PatientID <- c()
  association <- c() 
  valuePositive <- c() 
  valueNegative <- c() 
  
  # calculate percentage
  for(i in 28:48){
    sampleValues <- as.character(colnames(expressionPairs.annotated)[i])
    samplePredictions <- paste0(sampleValues,"_p")
    
    for(l in 1:length(listOfAnnotations)){
      PatientID <- c(PatientID, sampleValues )
      
      association <- c(association, listOfAnnotations[[l]])
      
      singleDF <- listOfDFs[[l]]
      highMethyl<- nrow(singleDF[singleDF[,samplePredictions] == "+",])
      lowMethyl <- nrow(singleDF[singleDF[,samplePredictions] == "-",])
      valuePositive <- c(valuePositive, (highMethyl / nrow(singleDF) ) ) 
      valueNegative <- c(valueNegative, (lowMethyl/ nrow(singleDF)) )
    }
    
  }
  
  # form barplot dataframe in right format
  barplotDFPos <- data.frame(PatientID, association, value = valuePositive, expression = "Increased")
  barplotDFNeg <- data.frame(PatientID, association, value =  -valueNegative, expression = "Decreased" )
  barplotDF <- rbind(barplotDFPos,barplotDFNeg)
  
  ##
  # rank the sample names according to their tissue type 
  ##
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
  barplotDF <- cbind(barplotDF,ranking)
  
  ##
  # add significance according to fisher test to each row 
  ##
  
  fisherTestDF <- fisher_test_expression(listOfAnnotations, listOfDFs)
  fisherTestDF[,3] <- p.adjust(fisherTestDF[,3] , method = "fdr")
  
  
  significancePos <- c() 
  significanceNeg <- c() 
  pval <- c()
  for(i in 1:length(barplotDF[,1])){
    association <- as.character(barplotDF[i,2])
    
    # Increased Expression 
    if(barplotDF[i,4]  == "Increased"){
      value <- as.numeric(as.character(fisherTestDF[fisherTestDF[,1] == as.character(barplotDF[i,1]) & 
                                                      fisherTestDF[,4] == association & 
                                                      fisherTestDF[,5] == "Increased",][2]))
      pv <- as.numeric(as.character(fisherTestDF[fisherTestDF[,1] == as.character(barplotDF[i,1]) & 
                                                   fisherTestDF[,4] == association & 
                                                   fisherTestDF[,5] == "Increased",][3]))
      
      
      if(pv < 0.05){
        significancePos <- c(significancePos,"X")
      } else {
        significancePos <- c(significancePos, "")
      }
      significanceNeg <- c(significanceNeg,"")
      
    }     
    if(barplotDF[i,4]  == "Decreased"){
      
      pv <- as.numeric(as.character(fisherTestDF[fisherTestDF[,1] == as.character(barplotDF[i,1]) & 
                                                   fisherTestDF[,4] == association & 
                                                   fisherTestDF[,5] == "Decreased",][3]))
      value <- as.numeric(as.character(fisherTestDF[fisherTestDF[,1] == as.character(barplotDF[i,1]) & 
                                                      fisherTestDF[,4] == association & 
                                                      fisherTestDF[,5] == "Decreased",][2]))
      
      
      if(pv < 0.05){
        significanceNeg <- c(significanceNeg,"X")
      } else {
        significanceNeg <- c(significanceNeg, "")
      }
      significancePos <- c(significancePos,"")
      
    } 
    
    
    pval <- c(pval, pv )

  }
  
  barplotDF <- cbind(barplotDF,significancePos,significanceNeg, pval)
  
  
  return(barplotDF)
  
}

fisher_test_expression <- function(listOfAnnotations, listOfDFs){
  fisherTestAll <- c(samples = colnames(expressionPairs.annotated)[28:48])
  
  for(l in 1:length(listOfDFs)){
    assocDdf <- listOfDFs[[l]]
    
    fisherTestAll <- cbind(fisherTestAll, create_fisher_all_samples(assocDdf, T,F), 
                           create_fisher_all_samples(assocDdf, T,T), 
                           create_fisher_all_samples(assocDdf, F,F), 
                           create_fisher_all_samples(assocDdf, F,T))
  }
  
  fishertestOR <- NULL
  
  for(i in 1:length(fisherTestAll[,1])){
    a <- 1
    for(j in seq(2,length(fisherTestAll[1,]),by = 4)){
      fishertestOR <- rbind(fishertestOR, cbind(samples = as.character(fisherTestAll[i,1]), 
                                                valueOR = fisherTestAll[i,j], valuePV = fisherTestAll[i,j + 1],
                                                Association = listOfAnnotations[[a]] , Expression = "Increased"))
      fishertestOR <- rbind(fishertestOR, cbind(samples = as.character(fisherTestAll[i,1]), 
                                                valueOR = fisherTestAll[i,j + 2], valuePV = fisherTestAll[i,j+3],
                                                Association = listOfAnnotations[[a]] , Expression = "Decreased"))
      a <- a + 1
    }
  }
  
  fishertestOR <- as.data.frame(fishertestOR)
  return(fishertestOR)
}

create_fisher_all_samples <- function(associationDF, positive, pvalue){
  
  if(pvalue){
    return (c(fisher_test_by_sample(2,associationDF,positive)$p.value,fisher_test_by_sample(3,associationDF,positive)$p.value,
              fisher_test_by_sample(4,associationDF,positive)$p.value,fisher_test_by_sample(5,associationDF,positive)$p.value,
              fisher_test_by_sample(6,associationDF,positive)$p.value,fisher_test_by_sample(7,associationDF,positive)$p.value,
              fisher_test_by_sample(8,associationDF,positive)$p.value,fisher_test_by_sample(9,associationDF,positive)$p.value,
              fisher_test_by_sample(10,associationDF,positive)$p.value,fisher_test_by_sample(11,associationDF,positive)$p.value,
              fisher_test_by_sample(12,associationDF,positive)$p.value,fisher_test_by_sample(13,associationDF,positive)$p.value,                                              
              fisher_test_by_sample(14,associationDF,positive)$p.value,fisher_test_by_sample(15,associationDF,positive)$p.value,
              fisher_test_by_sample(16,associationDF,positive)$p.value,fisher_test_by_sample(17,associationDF,positive)$p.value,
              fisher_test_by_sample(18,associationDF,positive)$p.value,fisher_test_by_sample(19,associationDF,positive)$p.value,
              fisher_test_by_sample(20,associationDF,positive)$p.value,fisher_test_by_sample(21,associationDF,positive)$p.value,
              fisher_test_by_sample(22,associationDF,positive)$p.value) )
  }
  if(!pvalue){
    return (c(fisher_test_by_sample(2,associationDF,positive)$estimate,fisher_test_by_sample(3,associationDF,positive)$estimate,
              fisher_test_by_sample(4,associationDF,positive)$estimate,fisher_test_by_sample(5,associationDF,positive)$estimate,
              fisher_test_by_sample(6,associationDF,positive)$estimate,fisher_test_by_sample(7,associationDF,positive)$estimate,
              fisher_test_by_sample(8,associationDF,positive)$estimate,fisher_test_by_sample(9,associationDF,positive)$estimate,
              fisher_test_by_sample(10,associationDF,positive)$estimate,fisher_test_by_sample(11,associationDF,positive)$estimate,
              fisher_test_by_sample(12,associationDF,positive)$estimate,fisher_test_by_sample(13,associationDF,positive)$estimate,                                              
              fisher_test_by_sample(14,associationDF,positive)$estimate,fisher_test_by_sample(15,associationDF,positive)$estimate,
              fisher_test_by_sample(16,associationDF,positive)$estimate,fisher_test_by_sample(17,associationDF,positive)$estimate,
              fisher_test_by_sample(18,associationDF,positive)$estimate,fisher_test_by_sample(19,associationDF,positive)$estimate,
              fisher_test_by_sample(20,associationDF,positive)$estimate,fisher_test_by_sample(21,associationDF,positive)$estimate,
              fisher_test_by_sample(22,associationDF,positive)$estimate) )
  }
  
}



fisher_test_by_sample <- function(sampleNr, associationDF, positive){
  
  # number of all, all that are positive and all that are negative in this sample 
  numAll <- sum(expressionPairs.annotated[,sampleNr] == "+") + sum(expressionPairs.annotated[,sampleNr] == "-") + sum(expressionPairs.annotated[,sampleNr] == "=")  
  numAllPos <- sum(expressionPairs.annotated[,sampleNr]== "+")  
  numAllNeg <- sum(expressionPairs.annotated[,sampleNr]== "-")
  
  # number of all associated and all positive and negative associated 
  numAssoc <- sum(associationDF[,sampleNr] == "+") +
    sum(associationDF[,sampleNr] == "-") + 
    sum(associationDF[,sampleNr] == "=")
  numAssocPos <- sum(associationDF[,sampleNr] == "+")
  numAssocNeg <- sum(associationDF[,sampleNr] == "-")
  
  ##
  # Negative expression
  ##
  if(! positive){
    # associated with region; in - state
    obenLinks <- numAssocNeg
    
    # not associated with region: in - state
    obenRechts <- numAllNeg - numAssocNeg
    
    # associated with region; not in - state
    untenLinks <- numAssoc - numAssocNeg
    
    # not associated with region; not in - state
    untenRechts <- (numAll - numAllNeg) - (numAssoc - numAssocNeg)
  }
  
  ##
  # Positive expression
  ##
  if(positive){
    # associated with region; in + state
    obenLinks <- numAssocPos
    
    # not associated with region: in - state
    obenRechts <- numAllPos - numAssocPos
    
    # associated with region; not in - state
    untenLinks <- numAssoc - numAssocPos
    
    # not associated with region; not in - state
    untenRechts <- (numAll - numAllPos) - (numAssoc - numAssocPos)
  }
  
  fisherTestMatrix <- matrix(c(obenLinks,untenLinks,obenRechts, untenRechts), ncol = 2 , nrow = 2)
  return(fisher.test(fisherTestMatrix, alternative = "greater"))
  
  
}



#############################################################################################################
##
### function that counts HMM prediction differences between two pairs 
###
##
############################################################################################################




countDiffsBetweenPairs <- function(pat1, pat2){
  predictionPat1 <- paste0(pat1,"_p")
  predictionPat2 <- paste0(pat2,"_p")
  
  col1 <- expressionPairs.annotated[,predictionPat1]
  col2 <- expressionPairs.annotated[,predictionPat2]
  comparison <- paste0(pat1, " vs ", pat2)
  
  return(c(comparison, 
           sum((col1 == "+" & col2 == "+")),
           sum((col1 == "-" & col2 == "-")),
           sum((col1 == "=" & col2 == "=")),
           sum((col1 == "+" & col2 == "=")),
           sum((col1 == "=" & col2 == "+")),
           sum((col1 == "-" & col2 == "=")),
           sum((col1 == "=" & col2 == "-")),
           sum((col1 == "+" & col2 == "-")),
           sum((col1 == "-" & col2 == "+"))))
}



#############################################################################################################
##
### function that means HMM predictions from metastases of the same patient 
###
##
############################################################################################################



create_mean_from_runDataCols <- function(firstRow, secondRow, thirdRow = NULL, fourthRow = NULL){
  patientRun <- c()
  firstRow <- convert_runData_to_numeric(firstRow)
  secondRow <- convert_runData_to_numeric(secondRow)
  allRows <- data.frame(firstRow,secondRow)
  
  if(!is.null(thirdRow)){
    thirdRow <- convert_runData_to_numeric(thirdRow)
    allRows <- data.frame(firstRow,secondRow, thirdRow)
  }
  if(!is.null(fourthRow)){
    fourthRow <- convert_runData_to_numeric(fourthRow)
    allRows <- data.frame(firstRow,secondRow, thirdRow, fourthRow)
  } 
  
  
  patientRNr <- pbapply::pbapply(allRows, FUN = mean, MARGIN = 1)
  
  patientR <- patientRNr
  patientR[patientRNr > 0 ] <- "+"
  patientR[patientRNr < 0 ] <- "-"
  patientR[patientRNr == 0 ] <- "="
  
  return(patientR)
}



convert_runData_to_numeric <- function(data){
  data <- as.character(data)
  data[data == "+"] <- 1
  data[data == "-"] <- -1
  data[data == "="] <- 0
  data <- as.numeric(data)
  return(data)
}




#############################################################################################################
##
### function that calculated the number of candidate genes dependent on cut-off 
###
##
############################################################################################################


get_number_of_candidates_on_cutoff <- function(genesPatientDF){
  numCand <- NULL
  
  
  for(i in 0:16){
    numCand <- rbind(numCand, c("cutOff" = i, "genesCount" = length( genesPatientDF[genesPatientDF$countNeg >= i,1] ), "Expression" = "Decreased"))
    numCand <- rbind(numCand, c("cutOff" = i, "genesCount" = length( genesPatientDF[genesPatientDF$countPos >= i,1] ), "Expression" = "Increased"))
  }
  numCand <- as.data.frame(numCand)
  numCand$genesCount <- as.numeric(as.character(numCand$genesCount))
  numCand$cutOff<- as.numeric(as.character(numCand$cutOff))
  
  
  return(numCand)
}



#############################################################################################################
##
### function that does Kaplan Meier analysis and prepare for multiple testing adjustment 
###
##
############################################################################################################





perform_survival_analysis_multipleTesting <- function(lowGroupIds, highGroupIds,geneID ){
  TCGA_survival_df$Group <- "no Group"
  TCGA_survival_df$Group[TCGA_survival_df$Name %in% lowGroupIds] <- "low expression"
  TCGA_survival_df$Group[TCGA_survival_df$Name %in% highGroupIds] <- "high expression"
  
  TCGA_survival_df <- TCGA_survival_df[TCGA_survival_df$Group != "no Group",]
  fit <- survfit(Surv(time = TCGA_survival_df[,2], 
                      event = TCGA_survival_df[,3]) ~ Group, 
                 data = TCGA_survival_df)
  text <- survdiff(Surv(TCGA_survival_df[,2], TCGA_survival_df[,3]) ~ Group)$pvalue
  return(list( fit, text) )
}


