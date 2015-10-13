#' predictor.R
#' 
#' predicts the ALSFRS slope for new patient using ensemble method of Gradient Boosting and Bayesian Additive Regression Tree (BART)
#' 
#' Author: Yeeleng Scott Vang (ysvang@uci.edu)
#================================================================================================================


rm(list=ls()) # clear all variables from environment

library(stats)
library(xgboost)
library(caret)
library(BayesTree)
source("dataPreprocessing.R")
source("textProcessingFunctions.R")

args <- commandArgs(trailingOnly=TRUE)

set.seed(50) # set seed arbitrary for reproducibility

con <- file(args[1]) 
#con <- file('patientOutput.txt') 
open(con);
clusterLine <- readLines(con, n = 1, warn = FALSE)
close(con)

s <- strsplit(clusterLine, "[^[:digit:]]")
solution <- as.numeric(unlist(s))
clusterNum <- unique(solution[!is.na(solution)])

# read in file without first row
#patient <- read.table("patientOutput.txt", sep="|", header=F, skip=1,comment.char="", quote="")
patient <- read.table(args[1], sep="|", header=F, skip=1,comment.char="", quote="")
colnames(patient) <- c("SubjectID", "form_name", "feature_name", "feature_value", "feature_unit", "feature_delta")

data <- read.csv("trainingDF.csv")
tDF <- read.csv("processedTrainingDF.csv")
targetY <- read.csv("featureReducedTarget.csv")
clusterCentroids <- read.csv("clusterCentroids.csv")
sizeFeatureNames <- read.csv("featureRank.csv")
featureReducedModel <- read.csv("featureReducedModel.csv")
tClusterIndices <- read.csv("tClusterIndices.csv")
tClusterIndices <- tClusterIndices[,1]

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
validFeatureSet <- colnames(newRow[,which(!is.na(newRow))])


# fill in appropriate NA cells for new Patient
data1 <- rbind(data,newRow)
newPatient <- newPatientImputation(data1)

# convert newPatient non-numeric into numeric dataframe
x <- as.numeric(as.character(newPatient))
tempX <- matrix(x, nrow=1, ncol=length(newPatient))
newX <- data.frame(tempX)
colnames(newX) <- colnames(newPatient)
newPatientNorm <- predict(preProc, newX) # apply centering/scaling to new patient


clusterModels <- split(featureReducedModel, f=tClusterIndices) # split training Dataset according to cluster index
targetModels <- split(targetY, f=tClusterIndices)  # split targets according to cluster index

# select subset of clustered data only
clusterData <- clusterModels[[clusterNum]]
targetData <- targetModels[[clusterNum]]


# extract only relevant features corresponding to feature selected model
newPatientNorm <- newPatientNorm[,colnames(featureReducedModel)]

# select only features available from 6 original feature and those available from reduced featured data model
featureSetIntersection <- intersect(validFeatureSet, colnames(newPatientNorm))
clusterData <- clusterData[, featureSetIntersection]
newPatientNorm <- newPatientNorm[, featureSetIntersection]

#Grid search for best gradient boosting parameters
etaOption <- c(.7, .6, .5, .4, .3, .2)
depthOption <-c(3,4,5)
bestParam <- c(1,1) # saves best parameters
old_test_error <- 10
for (i in 1:length(etaOption)){
  for (j in 1:length(depthOption)){
    param <- list("objective" = "reg:linear",
                  "eval_metric" = "rmse",
                  eta = etaOption[i], max.depth = depthOption[j])
    
    cv.nround <- 10
    
    bst.cv = xgb.cv(param=param, data = as.matrix(clusterData), label = targetData[,1],
                    nfold = 10, nrounds = cv.nround, verbose = F)
    if (tail(bst.cv$test.rmse.mean, n=1) < old_test_error){
      bestParam[1] <- etaOption[i]
      bestParam[2] <- depthOption[j]
      old_test_error <- tail(bst.cv$test.rmse.mean, n=1)
    }
  }
}

# gradient boosting predictor using best parameters from CV
bst <- xgboost(data = as.matrix(clusterData), label = targetData[,1], max.depth = bestParam[2], eta = bestParam[1], nround = 6,
               objective = "reg:linear", verbose = F)
bstPred <- predict(bst, as.matrix(newPatientNorm))

# Bayesian Additive Regression Tree (BART) predictor
temp <- rbind(newPatientNorm, newPatientNorm) # quick fix around dimension issue when supplying only one test patient row
# Averages 10 run of BART training to reduce variance
bartSum <- 0
for (b in 1:10){
  bartFit <- bart(clusterData, targetData[,1], temp, ndpost = 10000, ntree = 100, sigquant = 0.5, verbose = F)
  pred <- bartFit$yhat.test.mean
  bartSum <- bartSum + pred[1]
}
bartPred <- bartSum / 10

# ensemble predictor weighing gradient boosting and BART equally
ensemblePred <- 0.5*bstPred + 0.5*bartPred

predLine <- paste(ensemblePred, "|1")
#write(predLine, file="predictionOutput.txt")
write(predLine, file=args[2])




