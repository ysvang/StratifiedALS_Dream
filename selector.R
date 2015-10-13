#' selector .R
#' 
#' finds cluster to which new patient belongs to and selects 6 features for that patient
#' 
#' Author: Yeeleng Scott Vang (ysvang@uci.edu)
#================================================================================================================


rm(list=ls()) # clear all variables from environment

library(stats)
library(caret)
source("textProcessingFunctions.R")
source("dataPreprocessing.R")

args <- commandArgs(trailingOnly = TRUE)

set.seed(50) # set seed arbitrary for reproducibility

k = 16 # number of clusters
data <- read.csv("trainingDF.csv")
tDF <- read.csv("processedTrainingDF.csv")
targetY <- read.csv("processedTrainingTargetSlope.csv")
clusterCentroids <- read.csv("clusterCentroids.csv")
sizeFeatureNames <- read.csv("featureRank.csv")
featureReducedModel <- read.csv("featureReducedModel.csv")

patient <- read.table(args[1], sep="|", header=T, comment.char="", quote="")
#patient <- read.table("testPatient.txt", sep="|", header=T, comment.char="", quote="")

# center and normalize trainingDF
preProc  <- preProcess(tDF, method=c("center", "scale")) # calculate mean and variance of training set
tDF <- predict(preProc, tDF) # apply centering/scaling to training set

# initialize row for new patient
numRowData <- nrow(data)
tempRow <- matrix(c(rep.int(NA,length(data))), nrow=1, ncol=length(data))
newRow <- data.frame(tempRow)
colnames(newRow) <- colnames(data)

# process new patient data into a useable row in dataframe
newRow <- processValidateText(patient, data)

# fill in appropriate NA cells for new Patient
data1 <- rbind(data,newRow)
newPatient <- newPatientImputation(data1)

# convert newPatient non-numeric into numeric dataframe
x <- as.numeric(as.character(newPatient))
tempX <- matrix(x, nrow=1, ncol=length(newPatient))
newX <- data.frame(tempX)
colnames(newX) <- colnames(newPatient)
newPatientNorm <- predict(preProc, newX) # apply centering/scaling to new patient

# extract only relevant features
newPatientNorm <- newPatientNorm[,colnames(featureReducedModel)]

# replicate new patient k-times for matrix arithmetic
newPatientDF <- newPatientNorm[rep(seq_len(nrow(newPatientNorm)), each=k),]
# determines cluster via L1 norm
clusterDiff <- clusterCentroids - newPatientDF
clusterDiffAbs <- abs(clusterDiff)
L1norm <- rowSums(clusterDiffAbs)
minClusterIndex <- which.min(L1norm)

clusterFeature <- sizeFeatureNames[,minClusterIndex]
for (j in nrow(patient):1){
  if (!(toupper(patient$feature_name[j]) %in% toupper(clusterFeature))){
    patient <- patient[-j,]
  }
}



# write out patient file
clusterLine <- paste("cluster: ", minClusterIndex)
write(clusterLine, file=args[2])
#write(clusterLine, file="patientOutput.txt")
for(i in 1:nrow(patient)){
  write(paste(patient[i,1],patient[i,2],patient[i,3],patient[i,4],patient[i,5],patient[i,6],sep="|"), file=args[2],append=T)
  #write(paste(patient[i,1],patient[i,2],patient[i,3],patient[i,4],patient[i,5],patient[i,6],sep="|"), file="patientOutput.txt",append=T)
}


