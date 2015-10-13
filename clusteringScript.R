#' clusteringScript .R
#' 
#' Reads in a preprocesssed data frame and stratifies it into k clusters
#' 
#' Author: Yeeleng Scott Vang (ysvang@uci.edu)
#================================================================================================================


library(stats)
library(caret)
library(MASS)

# step-AIC feature selection option
aicFeatureSelection <- function(model){
  modelTemp <- stepAIC(model, trace=F, direction="both") # step AIC
  return (model[, attr(modelTemp$term, "term.labels")])
}

# sbf feature selection option
sbfFeatureSelection <- function(model, target){
  fit1 <- sbf(model, target, sbfControl =
                    sbfControl(functions = lmSBF, method = "repeatedcv", repeats = 10))
  return(model[, fit1$optVariables]) 
  
}

#' Function Name: clusterFunc
#' input: k
#' output: (none)
#' ********************************************************************
#' This function takes in the desired number of clusters, automatically reads in processed data files
#' generated upstream, and stratifies the training patients via k-means
#' ********************************************************************
clusterFunc <- function(k){
  # reads in processed data files
  tDF <- read.csv("processedTrainingDF.csv")
  targetY <- read.csv("processedTrainingTargetSlope.csv")
  
  preProc  <- preProcess(tDF, method=c("center", "scale")) # calculate mean and variance of training set
  tDF_Norm <- predict(preProc, tDF) # normalize training set
  
  #clusterCentroids <- matrix(c(rep.int(NA,length(tDF))), nrow=k, ncol=length(tDF)) # centroids of clusters
  #tClusterIndices <- rep(NA, nrow(tDF)) # centroid indices of each patients
  
  targetY <- targetY[,1] # transform target dataframe into a vector
  
  # initial feature selection to reduce feature set to avoid curse of dimensionality
  featureReducedModel <- sbfFeatureSelection(tDF_Norm, targetY)
  # perform k-means clustering
  c1 <- kmeans(featureReducedModel, k, nstart = 25, iter.max=200)
  
  # re-run k-means if one cluster consists of a single patient only
  while (1 %in% c1$size){
    oneDatapointCluster <- which(c1$size == 1)
    rowNumbers <- which(c1$cluster == oneDatapointCluster)
    featureReducedModel <- featureReducedModel[-rowNumbers,]
    targetY <- targetY[-rowNumbers]
    c1 <- kmeans(featureReducedModel, k, nstart = 25, iter.max=200)
  }
  
  clusterCentroids <- c1$centers # centroids of clusters
  tClusterIndices <- c1$cluster # index of each patients' centroid membership
  clusterModels <- split(featureReducedModel, f=tClusterIndices) # split training Dataset according to patient centroid membership
  targetModels <- split(targetY, f=tClusterIndices) # split target vector according to patient centroid membership
  
  # initialize a dataframe to track feature ranking of each cluster
  featureRank <- data.frame(matrix("", ncol=k, nrow=ncol(tDF_Norm)), stringsAsFactors=F)
  
  # feature ranking for individual clusters
  for (i in 1:k){
    if (nrow(clusterModels[[i]]) < 2) {
      # associate default features to cluster with less than 2 patients
      featureRank[1:6, i] <- c("Gender", "Age", "onset_delta", "if_use_Riluzole", "Q1_Speech", "Q2_Salivation") 
    } else {
    
      # rank remaining features 
      model <- train(as.matrix(targetModels[[i]]) ~ ., data=clusterModels[[i]], method="lm")
      varRanking <- varImp(model, scale=F)
      var <- varRanking[[1]]
      varRanked <- var[with(var, order(-Overall)), , drop=F]
      featureNames <- rownames(varRanked)
      
      # save ranked features
      featureNames <- rownames(varRanked)
      featureRank[1:length(featureNames),i] <- featureNames
    }
  }
  
  # process each feature name to obtain original feature name (relevant to one-hot encoded categorical features)
  originalFeatureNames <- data.frame(matrix("", ncol=k, nrow=ncol(tDF_Norm)), stringsAsFactors=F)
  for (i in 1:nrow(featureRank)){
    for (j in 1:ncol(featureRank)){
      if (!is.null(featureRank[i,j])){
        strings <- strsplit(featureRank[i,j], "[.]")
        if (length(strings[[1]]) > 1) {
          if (strings[[1]][1] != "tDF"){
            originalFeatureNames[i,j] <- strings[[1]][1]
          } else {
            originalFeatureNames[i,j] <- strings[[1]][2]
          }
        } else {
          originalFeatureNames[i,j] <- strings[[1]][1]
        }
      }
    }
  }
  
  # returns the 6 top feature for each clusters
  sixFeatureNames <- data.frame(matrix("", ncol=k, nrow=6), stringsAsFactors=F)
  for (i in 1:k){
    colFeatureSet = c()
    for (j in 1:nrow(originalFeatureNames)){
      if ((!is.na(originalFeatureNames[j,i])) && (!(originalFeatureNames[j,i] %in% colFeatureSet))){
        name <- originalFeatureNames[j,i]
        if (toupper(substr(name, 1, 4)) == "RACE"){
          name <- "Race"
        } else if (toupper(substr(name, 1, 10)) == "ONSET_SITE"){
          name <- "onset_site"
        } else {
          name <- name
        }
        colFeatureSet <- union(colFeatureSet, name)
      }
      if (length(colFeatureSet) >= 6){
        break
      }
    }
    
    for(w in 1:length(colFeatureSet)){
      sixFeatureNames[w,i] <- colFeatureSet[w]
    }
  }
  
  # write out support files for downstream scripts
  write.csv(clusterCentroids, file="clusterCentroids.csv", row.names=F)
  write.csv(sixFeatureNames, file="featureRank.csv", row.names=F)
  write.csv(tClusterIndices, file="tClusterIndices.csv", row.names=F)
  write.csv(featureReducedModel, file="featureReducedModel.csv", row.names=F)
  write.csv(targetY, file="featureReducedTarget.csv", row.names=F)
}



