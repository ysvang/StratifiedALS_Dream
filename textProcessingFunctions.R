#' textProcessingFunctions.R
#' 
#' functions to perform filtering of the unstructured text files and returns usable dataframe
#'
#' Author: Yeeleng Scott Vang (ysvang@uci.edu)
#================================================================================================================

# static features to use 
staticFeatures <- c("AGE", "GENDER", "RACE", "ONSET_DELTA", "DIAG_DELTA", "ONSET_SITE", "IF_USE_RILUZOLE", 
                    "TREATMENT_GROUP", "FAMILY_ALS_HIST")

# Form names of time-dependent features to use
timeDependentFormNames <- c("ALSFRS", "FVC", "SVC", "VITALS")


#' Function Name: processPROACT_Text
#' input: trainingFile (table)
#' output: data frame
#' ********************************************************************
#' This function takes in the PROACT table and filters the it into a 
#' data frame with each rows corresponding to a patient and each columns
#' corresponding to a feature
#' ********************************************************************
processPROACT_Text <- function(trainingFile) {
  featureSet <- c()
  patientSet <- c()
  dataMatrix <- data.frame()
  featureCol <- data.frame()
  numRowsTrainingFile <- nrow(trainingFile)
  i <- 0 # increments for every new feature encounter
  j <- 1 # tracks current row
  k <- j+1 # tracks all row that pertains to time series data of patient j
  
  while (j <= numRowsTrainingFile) {
    
    # static features processing
    if (toupper(trainingFile$feature_name[j]) %in% staticFeatures){
      
      # for age feature, take the floor
      if (toupper(trainingFile$feature_name[j]) == "AGE"){
        trainingFile$feature_value[j] <- floor(as.numeric(as.character(trainingFile$feature_value[j])))
      }
      
      # if new feature, then add previous feature column to data matrix
      if (!(trainingFile$feature_name[j] %in% featureSet)){ 
        if (i == 1){ 
          dataMatrix <- featureCol
          x <- featureCol
        } 
        if (i > 1) {
          dataMatrix <- merge(dataMatrix, featureCol, all=TRUE, row.names=FALSE) # add previous column to dataMatrix
        }
        featureName <- trainingFile$feature_name[j]  
        featureCol <- data.frame(id=c(trainingFile$SubjectID[j]), value=trainingFile$feature_value[j])  
        names(featureCol) <- c("id", toString(featureName))
        featureSet <- union(featureSet, featureName) # keep set of feature names
        i <- i + 1
      } else { # else add new patient to active feature column
        featureName <- trainingFile$feature_name[j]  
        newPatient <- data.frame(id=c(trainingFile$SubjectID[j]), value=trainingFile$feature_value[j])
        names(newPatient) <- c("id", toString(featureName))
        featureCol <- rbind(featureCol, newPatient)  # add new patient (row) to feature data (column) 
      }
      oldFeatureName <- trainingFile$feature_name[j]
    }
    
    # time series features processing
    if (toupper(trainingFile$form_name[j]) %in% timeDependentFormNames){ 
      # find last row of patient's data for given time-series feature (and guard against EOF)
      while((trainingFile$SubjectID[j] == trainingFile$SubjectID[k]) && (k <= numRowsTrainingFile)){
        k <- k + 1
      }
      
      # grab chuck of time series data to do statistical reduction
      if (toupper(trainingFile$feature_name[j]) == "ALSFRS_TOTAL" || toupper(trainingFile$feature_name[j]) == "ALSFRS_R_TOTAL"){
        # do nothing with these two feature
      } else {
        # grab chuck of time series data corresponding to patient j
        timeSeriesChunk <- trainingFile[j:(k-1),]
        # process chuck to remove irrelevant rows
        timeSeriesChunk <- timeSeriesChunk[!is.na(timeSeriesChunk$feature_delta),] # removes row with no delta value
        timeSeriesChunk <- timeSeriesChunk[!is.na(timeSeriesChunk$feature_value),] # removes row with no feature value
        timeSeriesChunk <- timeSeriesChunk[with(timeSeriesChunk, order(feature_delta)),] # order rows by ascending delta value
        timeSeriesChunk  <- timeSeriesChunk[as.numeric(as.character(timeSeriesChunk$feature_delta)) <= 91,] # keep rows with data less than 4 months
        
        # perform some statistic calculations
        if (nrow(timeSeriesChunk) > 0){
          featureValues<- as.numeric(as.character(timeSeriesChunk$feature_value)) # convert factor into numeric
          featureDelta <- as.numeric(as.character(timeSeriesChunk$feature_delta))
          
          # simple statistics
          maxValue <- max(featureValues)
          minValue <- min(featureValues)
          meanValue <- mean(featureValues)
          lastValue <- tail(featureValues, n=1)
          
          derivativeArray <- rep(NA, (length(featureDelta)-1)) # preallocate slope array
          # derivative time series
          if (length(featureDelta) == 1) {
            derivativeArray <- NA
            slope <- NA
          } else {
            for(z in 1:(length(featureDelta)-1)){
              derivativeArray[z] <- (featureValues[z+1] - featureValues[z])/(featureDelta[z+1] - featureDelta[z])
            }
            slope <- (tail(featureValues, n=1) - featureValues[1]) / (tail(featureDelta, n=1) - featureDelta[1])
          }
          
          maxSlope <- max(derivativeArray)
          minSlope <- min(derivativeArray)
          meanSlope <- mean(derivativeArray)
          lastSlope <- tail(derivativeArray, n=1)
          
          
          # if new feature, then add previous feature column to data matrix
          if (!(trainingFile$feature_name[j] %in% featureSet)){ 
            if (i > 1) {
              dataMatrix <- merge(dataMatrix, featureCol, all=TRUE) # add previous column to dataMatrix
            }
            featureName <- trainingFile$feature_name[j]
            featureSet <- union(featureSet, featureName) # keep set of feature names
            featureCol <- data.frame(id=c(trainingFile$SubjectID[j]), maxValue, minValue, meanValue, lastValue, maxSlope,
                                     minSlope, meanSlope, lastSlope, slope)
            names(featureCol) <- c("id", paste(featureName, "maxValue", sep="."), paste(featureName, "minValue", sep="."),
                                   paste(featureName, "meanValue", sep="."), paste(featureName, "lastValue", sep="."),
                                   paste(featureName, "maxSlope", sep="."), paste(featureName, "minSlope", sep="."),
                                   paste(featureName, "meanSlope", sep="."), paste(featureName, "lastSlope", sep="."),
                                   paste(featureName, "slope", sep="."))
            i <- i + 1
          } else { # else add patient to active feature columns
            newPatient <- data.frame(id=c(trainingFile$SubjectID[j]), maxValue, minValue, meanValue, lastValue, maxSlope,
                                     minSlope, meanSlope, lastSlope, slope)
            names(newPatient) <- c("id", paste(featureName, "maxValue", sep="."), paste(featureName, "minValue", sep="."),
                                   paste(featureName, "meanValue", sep="."), paste(featureName, "lastValue", sep="."),
                                   paste(featureName, "maxSlope", sep="."), paste(featureName, "minSlope", sep="."),
                                   paste(featureName, "meanSlope", sep="."), paste(featureName, "lastSlope", sep="."),
                                   paste(featureName, "slope", sep="."))
            featureCol <- rbind(featureCol, newPatient)  # add new patient (row) to feature data (column)  
          }
          oldFeatureName <- trainingFile$feature_name[j] 
        }
      } 
    }
    patientSet <- union(patientSet, as.integer(as.character(trainingFile$SubjectID[j]))) # set of patientIDs
    k <- k + 1
    j <- k - 1  
  }
  dataMatrix <- merge(dataMatrix, featureCol, all=TRUE) # merge last processed feature to training matrix
  return (dataMatrix)
}


#' Function Name: processValidateText
#' input: trainingFile (table), dataMatrixOld (data frame) 
#' output: data frame
#' ********************************************************************
#' This function takes in the validation and leaderboard training files and 
#' filters them into a data frame with each rows corresponding to a patient and each columns
#' corresponding to a feature, with the features matching those of dataMatrixOld
#' ********************************************************************
processValidateText <- function(trainingFile, dataMatrixOld){
  colnames(trainingFile)[1] <- "id"

  # extract table subset corresponding to static features
  staticRows <- trainingFile[toupper(as.character(trainingFile$feature_name)) %in% staticFeatures, ]
  # extract table subset corresponding to time dependent features
  timeDependentRows <- trainingFile[toupper(as.character(trainingFile$form_name)) %in% timeDependentFormNames, ]
  # further process training file to remove irrelevant rows
  trainingFile <- rbind(staticRows, timeDependentRows) 
  trainingFile <- trainingFile[!is.na(trainingFile$feature_delta),] # removes row with no delta value
  trainingFile <- trainingFile[!is.na(trainingFile$feature_value),] # removes row with no feature value
  trainingFile <- trainingFile[as.numeric(as.character(trainingFile$feature_delta)) <= 91,] # keep rows with data less than 4 months only
  
  numRowData <- nrow(trainingFile)
  newPatients <- unique(trainingFile$id)
  # initialze new data frame for validation and leaderboard patients
  tempRows <- matrix(c(rep.int(NA,length(dataMatrixOld))), nrow=length(newPatients), ncol=length(dataMatrixOld))
  dataMatrixNew <- data.frame(tempRows)
  colnames(dataMatrixNew ) <- colnames(dataMatrixOld)
  
  for (j in 1:length(newPatients)) {
    # grab only data rows relevant to patient j
    currentPatient <- trainingFile[as.numeric(as.character(trainingFile$id)) == newPatients[j],]
    # set of static feature for which patient j has data for
    currentPatientStaticFeatureSet <- intersect(toupper(unique(currentPatient$feature_name)), staticFeatures)
    # set of form_name for which patient j has data for
    currentPatientUniqueFormSet <-  intersect(toupper(unique(currentPatient$form_name)), timeDependentFormNames)
    
    # process static feature for patient j
    k <- 1
    while (k <= length(currentPatientStaticFeatureSet)){  
      currentPatientCluster <- currentPatient[toupper(as.character(currentPatient$feature_name)) == currentPatientStaticFeatureSet[k], ]
      
      # for age feature, take the floor
      if (toupper(currentPatientStaticFeatureSet[k]) == "AGE"){
        currentPatientCluster$feature_value <- floor(as.numeric(as.character(currentPatientCluster$feature_value)))
      }
      
      # fill in dataMatrixNew with with data at relevant columns
      for (m in 2:ncol(dataMatrixNew)){
        if (colnames(dataMatrixNew[m]) == currentPatientCluster$feature_name) {
          #print(j)
          dataMatrixNew[j,m] <- as.character(currentPatientCluster$feature_value[1])
        }
      }      
      k <- k + 1
    }

    
    # process time dependent series feature for patient j
    k <- 1
    while (k <= length(currentPatientUniqueFormSet)){  
      currentPatientCluster <- currentPatient[toupper(as.character(currentPatient$form_name)) == currentPatientUniqueFormSet[k], ]
      # find last row of patient's data for given time-series feature (and guard against EOF)
      uniqueTimeDependentFeatureSet <- unique(currentPatientCluster$feature_name)
      
      for (w in 1:length(uniqueTimeDependentFeatureSet)){
        currentPatientFeatureSubcluster <- currentPatientCluster[as.character(currentPatientCluster$feature_name) == uniqueTimeDependentFeatureSet[w], ]
        if (toupper(currentPatientFeatureSubcluster$feature_name[1]) == "ALSFRS_TOTAL" || toupper(currentPatientFeatureSubcluster$feature_name[1]) == "ALSFRS_R_TOTAL"){
          # Do nothing for these two features
        } else {
          # grab chuck of time series data to do statistical reduction
          currentPatientFeatureSubcluster <- currentPatientFeatureSubcluster[with(currentPatientFeatureSubcluster, order(feature_delta)),] # order rows by ascending delta value
         
          if (nrow(currentPatientFeatureSubcluster) > 0){ 
            featureValues<- as.numeric(as.character(currentPatientFeatureSubcluster$feature_value)) # convert factors into numeric
            featureDelta <- as.numeric(as.character(currentPatientFeatureSubcluster$feature_delta))
            
            # simple statistics
            maxValue <- max(featureValues)
            minValue <- min(featureValues)
            meanValue <- mean(featureValues)
            lastValue <- tail(featureValues, n=1)
            
            derivativeArray <- rep(NA, (length(featureDelta)-1)) # preallocate slope array
            if (length(featureDelta) == 1) {
              derivativeArray <- NA
              slope <- NA
            } else {
              for(z in 1:(length(featureDelta)-1)){
                derivativeArray[z] <- (featureValues[z+1] - featureValues[z])/(featureDelta[z+1] - featureDelta[z])
              }
              slope <- (tail(featureValues, n=1) - featureValues[1]) / (tail(featureDelta, n=1) - featureDelta[1])
            }
            
            maxSlope <- max(derivativeArray)
            minSlope <- min(derivativeArray)
            meanSlope <- mean(derivativeArray)
            lastSlope <- tail(derivativeArray, n=1)
            featureName <- currentPatientFeatureSubcluster$feature_name[1]
                 
            for (m in 1:ncol(dataMatrixNew )){
              if (colnames(dataMatrixNew[m]) == paste(featureName, "maxValue", sep=".")) {
                dataMatrixNew [j,m] <- maxValue
              }
              if (colnames(dataMatrixNew[m]) == paste(featureName, "minValue", sep=".")) {
                dataMatrixNew[j,m] <- minValue
              }
              if (colnames(dataMatrixNew[m]) == paste(featureName, "meanValue", sep=".")) {
                dataMatrixNew[j,m] <- meanValue
              }
              if (colnames(dataMatrixNew[m]) == paste(featureName, "lastValue", sep=".")) {
                dataMatrixNew[j,m] <- lastValue
              }
              if (colnames(dataMatrixNew[m]) == paste(featureName, "maxSlope", sep=".")) {
                dataMatrixNew[j,m] <- maxSlope
              }
              if (colnames(dataMatrixNew[m]) == paste(featureName, "minSlope", sep=".")) {
                dataMatrixNew[j,m] <- minSlope
              }
              if (colnames(dataMatrixNew[m]) == paste(featureName, "meanSlope", sep=".")) {
                dataMatrixNew[j,m] <- meanSlope
              }
              if (colnames(dataMatrixNew[m]) == paste(featureName, "lastSlope", sep=".")) {
                dataMatrixNew[j,m] <- lastSlope
              }
              if (colnames(dataMatrixNew[m]) == paste(featureName, "slope", sep=".")) {
                dataMatrixNew[j,m] <- slope
              }
            }
          }
        }
      } 
      k <- k + 1
    }
    dataMatrixNew[j,1] <- as.numeric(as.character(currentPatient$id[1]))
  }
  dataMatrixNew$Age <- as.numeric(as.character(dataMatrixNew$Age))
  dataMatrixNew$onset_delta <- as.numeric(as.character(dataMatrixNew$onset_delta))
  dataMatrixNew$diag_delta <- as.numeric(as.character(dataMatrixNew$diag_delta))
  
  return (dataMatrixNew)
}
