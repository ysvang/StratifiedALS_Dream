#' textProcessingMain.R
#' 
#' Script performs data cleansing and generate all the necessary files needed to run selector.R and predictor.R
#'
#' Author: Yeeleng Scott Vang (ysvang@uci.edu)
#================================================================================================================


rm(list=ls()) # clear all variables from environment
source("textProcessingFunctions.R")
source("dataPreprocessing.R")
source("clusteringScript.R")

set.seed(50) # set seed arbitrary for reproducibility

# Reads in and process PROACT training file
inputFile <- "corrected\\all_forms_PROACT.txt"
trainingFile <- read.table(inputFile, sep="|", header=T, comment.char="", quote="")
dataMatrix_PROACT <- processPROACT_Text(trainingFile)

# Read in additional validation and leaderboard training files
inputFile2 <- "corrected\\all_forms_validate_spike.txt"
inputFile3 <- "leaderboard_data_for_releaseV2\\all_forms_validate_leader.txt"
trainingFile2 <- read.table(inputFile2, sep="|", header=T, comment.char="", quote="")
colnames(trainingFile2)[1] <- "id"
trainingFile3 <- read.table(inputFile3, sep="|", header=T, comment.char="", quote="")
colnames(trainingFile3)[1] <- "id"
validationLeaderData <- rbind(trainingFile2, trainingFile3) # combine these two training files before feading through preprocessor
dataMatrix_VL <- processValidateText(validationLeaderData, dataMatrix_PROACT)
dataMatrix <- rbind(dataMatrix_PROACT, dataMatrix_VL) # merge all training patients into one feature matrix

# Reads in target values of PROACT and additional training datasets
inputYFile <- c("corrected\\ALSFRS_slope_PROACT.txt", "corrected\\ALSFRS_slope_validate_spike.txt", 
                "leaderboard_data_for_releaseV2\\ALSFRS_slope_validate_leader2.txt")
for (i in 1:length(inputYFile)){
  slopeDF <- read.table(inputYFile[i], sep="|", header=T, comment.char="", quote="")
  colnames(slopeDF)[1] <- "id"
  #memory.limit(size=4000)
  if (i == 1){
    targetY <- slopeDF
  } else {
    targetY  <- rbind(targetY , slopeDF)
  }
}

# combine feature matrix and target matrix into one training dataframe
trainingDF <- merge(dataMatrix, targetY, all=TRUE)
trainingDF <- trainingDF[!is.na(trainingDF$ALSFRS_slope),] # remove rows with no ALSFRS slope value
# write out support files for downstream scripts
write.csv(trainingDF, file="trainingDF.csv", row.names=FALSE)

# imputes all NA cells in training dataframe
imputationCleanup()

# performs the stratification of training data frame
clusterFunc(16) # argument is number of clusters




