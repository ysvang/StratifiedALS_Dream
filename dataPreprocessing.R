#' dataPreprocessing.R
#' 
#' Reads in a data frame datastructure and performs data cleaning and NA imputation
#' 
#' Author: Yeeleng Scott Vang (ysvang@uci.edu)
#================================================================================================================


library(stats)

#' Function Name: imputationCleanup
#' input: (None)
#' output: (none)
#' ********************************************************************
#' This function takes in a csv with raw data and imputes data into cells missing data.
#' For categorical features, the most frequent feature is used to fill in missing cells.
#' For real-value features, the average across all patients is used to fill in missing cells.
#' Note: file automatically reads "trainingDF.csv" file generated from textProcessingMain.R,
#'       and outputs two additional CSV files used by other scripts down the pipeline.
#' ********************************************************************
imputationCleanup <- function(){
  tDF <- read.csv("trainingDF.csv")  # reads in data structure
  
  trainingY <- tDF$ALSFRS_slope  # save target columns into own variable
  tDF$ALSFRS_slope <- NULL
  
  # count the number of NA in each features
  numNAPerFeature <- rep(NA, ncol(tDF))
  for(i in 1:length(numNAPerFeature)){
    numNAPerFeature[i] <- (sum(is.na(tDF[i])))
  }
  
  # if more than 1000 NA in feature, delete the feature out of design matrix
  for (i in length(numNAPerFeature):1){
    if (numNAPerFeature[i] > 1000) {
      tDF <- tDF[,-i]
    }
  }
  
  
  # numNAPerFeature <- rep(NA, ncol(tDF))
  # for(i in 1:length(numNAPerFeature)){
  #   numNAPerFeature[i] <- (sum(is.na(tDF[i])))
  # }
  
  
  # Map Gender feature into binary values
  #table(tDF$Gender) # M: 1720, F: 982
  tDF$Gender[is.na(tDF$Gender)] <- names(which.max(table(tDF$Gender)))
  tDF$Gender <- sub("F", 0, tDF$Gender)
  tDF$Gender <- sub("M", 1, tDF$Gender)
  
  # replace Age with average if missing data
  avg_Age <- mean(tDF$Age, na.rm = TRUE) # find average of Age
  tDF$Age <- ifelse(is.na(tDF$Age), avg_Age, tDF$Age) # replace 'na' with average
  
  # replace NA values in onset_delta with the average
  tDF$onset_delta <- as.numeric(as.character(tDF$onset_delta)) # change datatype from factor to numeric...
  avg_onset_delta <- mean(tDF$onset_delta, na.rm = TRUE) # find average of feature
  tDF$onset_delta <- ifelse(is.na(tDF$onset_delta), avg_onset_delta, tDF$onset_delta) # replace 'na' with average
  
  # replace NA values of if_use_Riluzole with the most frequent value 'Yes' and map feature to binary values
  #table(tDF$if_use_Riluzole) # Yes: 1351, No: 709
  tDF$if_use_Riluzole[is.na(tDF$if_use_Riluzole)] <- names(which.max(table(tDF$if_use_Riluzole)))
  tDF$if_use_Riluzole <- sub("No", 0, tDF$if_use_Riluzole)
  tDF$if_use_Riluzole <- sub("Yes", 1, tDF$if_use_Riluzole)
  
  # replace NA values of treatment_group with the most frequent value 'Active' and map feature to binary values
  table(tDF$treatment_group) # Active: 1672, Placebo: 549
  tDF$treatment_group[is.na(tDF$treatment_group)] <- names(which.max(table(tDF$treatment_group)))
  tDF$treatment_group <- sub("Placebo", 0, tDF$treatment_group)
  tDF$treatment_group <- sub("Active", 1, tDF$treatment_group)
  
  # for remaining real-valued features, replace missing NA cell with average of that feature
  for(i in 9:ncol(tDF)){
    tDF[,i][is.infinite(tDF[,i])] <- NA
    tDF[,i] <- as.numeric(as.character(tDF[,i])) # change datatype from factor to numeric...
    avg_i <- mean(tDF[,i], na.rm = TRUE) # find average of column
    tDF[,i] <- ifelse(is.na(tDF[,i]), avg_i, tDF[,i]) # replace 'na' with average
  }
  
  # check to see there's no NA cell left in design matrix
  # numNAPerFeature <- rep(NA, ncol(tDF))
  # for(i in 1:length(numNAPerFeature)){
  #   numNAPerFeature[i] <- (sum(is.na(tDF[i])))
  # }
  
  
  tDF$id <- NULL # delete id column
  
  # use one-hot encoding of categorical feature Race
  table(tDF$Race)
  tDF$Race[is.na(tDF$Race)] <- names(which.max(table(tDF$Race)))
  tDF$Race <- factor(tDF$Race)
  Race_oneShot <- model.matrix(~tDF$Race+0)
  tDF$Race <- NULL
  tDF <- data.frame(tDF, Race_oneShot)
  
  # use one-hot encoding of categorical feature onset_site
  table(tDF$onset_site)
  tDF$onset_site[is.na(tDF$onset_site)] <- names(which.max(table(tDF$onset_site)))
  tDF$onset_site <- factor(tDF$onset_site)
  onset_site_oneShot <- model.matrix(~tDF$onset_site+0)
  tDF$onset_site <- NULL
  tDF <- data.frame(tDF, onset_site_oneShot)
  
  # write out support files for downstream scripts
  write.csv(tDF, file="processedTrainingDF.csv", row.names=FALSE)
  write.csv(trainingY, file="processedTrainingTargetSlope.csv", row.names=FALSE)
}


#' Function Name: newPatientImputation
#' input: tDF (data frame) 
#' output: vector
#' ********************************************************************
#' This function takes in raw tDF data frame and imputes data into cells missing data.
#' For categorical features, the most frequent feature is used to fill in missing cells.
#' For real-value features, the average across all patients is used to fill in missing cells.
#' ********************************************************************
newPatientImputation <- function(tDF){
  
  tDF$ALSFRS_slope <- NULL
  
  # count the number of NA in each features
  numNAPerFeature <- rep(NA, ncol(tDF))
  for(i in 1:length(numNAPerFeature)){
    numNAPerFeature[i] <- (sum(is.na(tDF[i])))
    
  }
  
  # if more than 1000 NA in feature, delete the feature out of design matrix
  for (i in length(numNAPerFeature):1){
    if (numNAPerFeature[i] > 1000) {
      tDF <- tDF[,-i]
    }
  }
  
  # Map Gender feature into binary values
  #table(tDF$Gender) # M: 1720, F: 982
  tDF$Gender[is.na(tDF$Gender)] <- names(which.max(table(tDF$Gender)))
  tDF$Gender <- sub("F", 0, tDF$Gender)
  tDF$Gender <- sub("M", 1, tDF$Gender)
  
  # replace Age with average if missing data
  avg_Age <- mean(tDF$Age, na.rm = TRUE) # find average of Age
  tDF$Age <- ifelse(is.na(tDF$Age), avg_Age, tDF$Age) # replace 'na' with average
  
  # replace NA values in onset_delta with the average
  tDF$onset_delta <- as.numeric(as.character(tDF$onset_delta)) # change datatype from factor to numeric...
  avg_onset_delta <- mean(tDF$onset_delta, na.rm = TRUE) # find average of feature
  tDF$onset_delta <- ifelse(is.na(tDF$onset_delta), avg_onset_delta, tDF$onset_delta) # replace 'na' with average
  
  # replace NA values of if_use_Riluzole with the most frequent value 'Yes' and map feature to binary values
  table(tDF$if_use_Riluzole) # Yes: 1351, No: 709
  tDF$if_use_Riluzole[is.na(tDF$if_use_Riluzole)] <- names(which.max(table(tDF$if_use_Riluzole)))
  tDF$if_use_Riluzole <- sub("No", 0, tDF$if_use_Riluzole)
  tDF$if_use_Riluzole <- sub("Yes", 1, tDF$if_use_Riluzole)
  
  # replace NA values of treatment_group with the most frequent value 'Active' and map feature to binary values
  table(tDF$treatment_group) # Active: 1672, Placebo: 549
  tDF$treatment_group[is.na(tDF$treatment_group)] <- names(which.max(table(tDF$treatment_group)))
  tDF$treatment_group <- sub("Placebo", 0, tDF$treatment_group)
  tDF$treatment_group <- sub("Active", 1, tDF$treatment_group)
  
  # for remaining real-valued features, replace missing NA cell with average of that feature
  for(i in 9:ncol(tDF)){
    tDF[,i][is.infinite(tDF[,i])] <- NA
    tDF[,i] <- as.numeric(as.character(tDF[,i])) # change datatype from factor to numeric...
    avg_i <- mean(tDF[,i], na.rm = TRUE) # find average of column
    tDF[,i] <- ifelse(is.na(tDF[,i]), avg_i, tDF[,i]) # replace 'na' with average
  }
  
  tDF$id <- NULL # delete id column
  
  # use one-hot encoding of categorical feature Race
  #table(tDF$Race)
  tDF$Race[is.na(tDF$Race)] <- names(which.max(table(tDF$Race))) # replace NA cell in column with most frequent value
  tDF$Race <- factor(tDF$Race)
  Race_oneShot <- model.matrix(~tDF$Race+0)
  tDF$Race <- NULL
  tDF <- data.frame(tDF, Race_oneShot)
  
  # use one-hot encoding of categorical feature onset_site
  table(tDF$onset_site)
  tDF$onset_site[is.na(tDF$onset_site)] <- names(which.max(table(tDF$onset_site))) # replace NA cell in column with most frequent value
  tDF$onset_site <- factor(tDF$onset_site)
  onset_site_oneShot <- model.matrix(~tDF$onset_site+0)
  tDF$onset_site <- NULL
  tDF <- data.frame(tDF, onset_site_oneShot)
  
  return (tail(tDF, n=1)) # return just the new patient back
}