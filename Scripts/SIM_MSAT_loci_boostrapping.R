# this script combines the 
rm(list = ls())
library(abind)
sim.wd <- 'C:/Users/gsalas/Documents/resampling_CIs/Code/'
setwd(sim.wd)
source('Morton_SSRvSNP_Simulations/RScripts/functions_SSRvSNP_Sim.R')
readGeninds_MSAT(paste0(sim.wd,'Datasets/MSAT_N1200/'))
## Building resampling array ##
gm_resamp_array_function <- function(insert_genind, num_loci, num_reps){
  # Create an empty array named 'resamp_category5loc' to store results
  resamp_category <- array(dim = c(nrow(insert_genind@tab),4,num_reps))
  # loop 25 times for sets (5 loci in one set) of randomly selected loci
  for (i in 1:num_reps) {
    # Randomly sample an amount of loci from the genind object based on user input
    samp_loc <- sample(locNames(insert_genind), size = num_loci, replace = FALSE)
    # Subset the genind object to include only the columns corresponding to the sampled loci
    gm.Wild.genind <- insert_genind[, loc = samp_loc]
    # declare objects
    # access the genind matrix that shows the type of alleles and quantity present among wild individuals
    wildSamp <- gm.Wild.genind@tab
    # calculate the sum of each column in the matrix, ignoring NA values.
    # identify the indices where the sum is not equal to zero, this indicates columns with 
    # variation in allele counts. Subset the original matrix by selecting only the columns identified
    # in the previous step, removing columns with no variation in allele counts
    wildSamp <- wildSamp[, which(colSums(wildSamp, na.rm = TRUE) != 0)]
    # calculating a vector of allele frequencies (sum of allele counts divided by number of haplotypes, or individuals * 2)
    wildComplete <- colSums(wildSamp, na.rm = TRUE) / (nrow(wildSamp) * 2)
    # Subset 'wildComplete' to include only the alleles with non-zero frequency
    wildSubset <- wildComplete[wildComplete > 0]
    # initialize vectors to store results 
    total <- vector(length = nrow(wildSamp))
    common <- vector(length = nrow(wildSamp))
    lowfreq <- vector(length = nrow(wildSamp))
    rare <- vector(length = nrow(wildSamp))
    # Loop for each sample size, ranging from 1 tto the number of loci in the wild population
    for (j in 1:nrow(wildSamp)) {
      # randomly select a subset of rows from the matrix, without replacement
      samp <- sample(nrow(wildSamp), size = j, replace = FALSE)
      # subset the original matrix, to inlcude only the rows randomly selected
      samp <- wildSamp[samp,]
      # Measure the proportion of allelic representation in that sample
      # Check if its the first iteration of the inner loop
      if (j == 1) {
        # Identify the names of alleles in the sample that are also present in 'wildSubset'. 
        # Calculate the proportion of alleles in the sample that are also present in 'wildSubset'.
        total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(samp > 0)))]) / length(names(wildSubset))
        # Identify the names of the alleles in the sample that are also present in the 'wildSubset' for alleles with frequency >0.1
        # Calculate the proportion of common alleles in the sample compared to the 'wildSubset' for alleles with frequency >0.1
        common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.1))
        # Identify the names of alleles in the sample that are also present in the 'wildSubset' for alleles with frequency >0.01 and frequency <0.1
        # Calculate the proportion of low frequency alleles in the sample compared to the 'wildSubset' for allels with frequency >0.01 and frequency <0.1
        lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
        # Identify the names of alleles in the sample that are also present in the 'wildSubset' for allleles with frequency <0.02.
        # Calcualte the proportion of rare alleles in the sample compared to the 'wildSubset' for alleles with frequency <0.01.
        rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(samp > 0)))) / length(which(wildSubset < 0.01))
      } else {
        # Calculate the proportion of alleles in the sample that are also present in the 'wildSubset' for the case where j is not 1
        total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))]) / length(names(wildSubset))
        # Calculate the proportion of alleles of common alleles in the sample compared toe the 'wildSubset' where j is not 1
        common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.1))
        # Calculate the proportion of alleles of low frequency alleles in the sample compared to the 'wildSubset' for the case where j is not 1
        lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
        # Calculate the proportion of rare alleles in the sample compared ot the 'wildSubset' for the case where j is not 1
        rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset < 0.01))
      }
      # Combine the results (proportions) for each sample size into a matrix named 'categorymat_10loc'.
      categorymat <- cbind(total, common, lowfreq, rare)
      # extract the column names (categories) from the matrix
      category <- colnames(categorymat)
      # set row and column names for 'resamp_category10loc'
      dimnames(resamp_category) <- list(paste0("sample ", 1:nrow(categorymat)), category, paste0("replicate ", 1:num_reps))
      # Store the results for the current replicate in resamp_category5loc
      resamp_category[, , i] <- categorymat
    }
  }
  return(resamp_category)
}

## Linear model ##
# pass all arrays to a dataframe using the function. define the function that takes an input 'data_array'
analyze_resampling_array <- function(data_array) {
  # linear model of resampling array. Extract the column named 'total' from the array.
  # concatenate the extracted column values into the vector.
  totalsVector <- c(data_array[,1,]) 
  
  # Specify sample numbers column.
  # Create a vector of sample numbers from 1 to the number of rows in the 'total' column 
  gm_sampleNumbers <- 1:(nrow(data_array[,1,]))
  # Repeat the sample numbers vector for the number of replicates
  gm_sampleNumbers <- rep(gm_sampleNumbers, dim(data_array)[[3]])
  
  # Create data frame from resampling array values
  gm_DF <- data.frame(sampleNumbers=gm_sampleNumbers, totalValues=totalsVector)
  
  # Build and analyze linear models
  gm_Model <- lm(sampleNumbers ~ I((totalValues)^3), data = gm_DF)
  # Create a new data fram 'gm_newData with a single column 'totalValues' containing the value 0.95
  gm_newData <- data.frame(totalValues=0.95)
  # Use the linear model to predict the response for the new data frame. Specify 'interval = prediction to obtain a prediction interval
  gm_95MSSEprediction <- predict(gm_Model, gm_newData, interval = "prediction")
  
  # Pass the gm_95MSSEprediction to the object storing our results 
  # Store the predicted values and predictino interval in the object named 'result' 
  result <- gm_95MSSEprediction
  # Calculate the width of the prediction interval by substracting the lower limit from the upper limit
  piWidth <- gm_95MSSEprediction[3] - gm_95MSSEprediction[2]
  # Return a list containing the predicted values and the width of the prediction interval
  return(list(result = result, piWidth = piWidth))
}

# ## Filling in matrix ##
# # Iterate through the arrays and store results in the matrix
# # Initiate loop that iterates over the indices of 'array_list'
# build_matrix_func <- function(array_list, input_matrix){
#   for (i in 1:length(array_list)) {
#     # Store results and piWidth values in the ith row of the matrix
#     input_matrix[i, ] <- c(array_list[[i]]$result, 
#                            array_list[[i]]$piWidth)
#   }
#   return(input_matrix)
# }
# 
# supremeArray <- array(dim = c(1200, 5, 25))
# 
# MSAT_levels <- c(5, 10, 15, 20, 25)
# 
# QUAC_array_list = list(length(MSAT_levels))
# 
# QUAC_predict_results <- list(length(MSAT_levels))
# 
# for (i in 1:length(MSAT_levels)) {
#   for (j in 1:length(MSAT_01pop_migHigh.genList)) {
#     QUAC_array_list <- gm_resamp_array_function(MSAT_01pop_migHigh.genList[[j]],MSAT_levels[i],5)
#   }
#   supremeArray[,,i] <- abind(QUAC_array_list[[i]])
#   # store the results into a list 
#   QUAC_predict_results[[i]] <- analyze_resampling_array(supremeArray)
#   print(QUAC_predict_results)
# }

##############
# 03/26/2024 #
##############
# This script generates a list "QUAC_array_list" of 5 arrays. Each array stores 
# representation values across replicates from each scenario based on a number of loci
# (5-25 loci are bootstrapped at intervals of 5)
rm(list = ls())
library(abind)
sim.wd <- 'C:/Users/gsalas/Documents/resampling_CIs/Code/'
setwd(sim.wd)
source('Morton_SSRvSNP_Simulations/RScripts/functions_SSRvSNP_Sim.R')
# readGeninds_MSAT(paste0(sim.wd,'Datasets/MSAT_N1200/'))
## Building resampling array ##
gm_resamp_array_function <- function(insert_genind, num_loci, num_reps){
  # Create an empty array named 'resamp_category5loc' to store results
  resamp_category <- array(dim = c(nrow(insert_genind@tab),4,num_reps))
  # loop 25 times for sets (5 loci in one set) of randomly selected loci
  for (i in 1:num_reps) {
    # Randomly sample an amount of loci from the genind object based on user input
    samp_loc <- sample(locNames(insert_genind), size = num_loci, replace = FALSE)
    # Subset the genind object to include only the columns corresponding to the sampled loci
    gm.Wild.genind <- insert_genind[, loc = samp_loc]
    # declare objects
    # access the genind matrix that shows the type of alleles and quantity present among wild individuals
    wildSamp <- gm.Wild.genind@tab
    # calculate the sum of each column in the matrix, ignoring NA values.
    # identify the indices where the sum is not equal to zero, this indicates columns with 
    # variation in allele counts. Subset the original matrix by selecting only the columns identified
    # in the previous step, removing columns with no variation in allele counts
    wildSamp <- wildSamp[, which(colSums(wildSamp, na.rm = TRUE) != 0)]
    # calculating a vector of allele frequencies (sum of allele counts divided by number of haplotypes, or individuals * 2)
    wildComplete <- colSums(wildSamp, na.rm = TRUE) / (nrow(wildSamp) * 2)
    # Subset 'wildComplete' to include only the alleles with non-zero frequency
    wildSubset <- wildComplete[wildComplete > 0]
    # initialize vectors to store results 
    total <- vector(length = nrow(wildSamp))
    common <- vector(length = nrow(wildSamp))
    lowfreq <- vector(length = nrow(wildSamp))
    rare <- vector(length = nrow(wildSamp))
    # Loop for each sample size, ranging from 1 tto the number of loci in the wild population
    for (j in 1:nrow(wildSamp)) {
      # randomly select a subset of rows from the matrix, without replacement
      samp <- sample(nrow(wildSamp), size = j, replace = FALSE)
      # subset the original matrix, to inlcude only the rows randomly selected
      samp <- wildSamp[samp,]
      # Measure the proportion of allelic representation in that sample
      # Check if its the first iteration of the inner loop
      if (j == 1) {
        # Identify the names of alleles in the sample that are also present in 'wildSubset'. 
        # Calculate the proportion of alleles in the sample that are also present in 'wildSubset'.
        total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(samp > 0)))]) / length(names(wildSubset))
        # Identify the names of the alleles in the sample that are also present in the 'wildSubset' for alleles with frequency >0.1
        # Calculate the proportion of common alleles in the sample compared to the 'wildSubset' for alleles with frequency >0.1
        common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.1))
        # Identify the names of alleles in the sample that are also present in the 'wildSubset' for alleles with frequency >0.01 and frequency <0.1
        # Calculate the proportion of low frequency alleles in the sample compared to the 'wildSubset' for allels with frequency >0.01 and frequency <0.1
        lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
        # Identify the names of alleles in the sample that are also present in the 'wildSubset' for allleles with frequency <0.02.
        # Calcualte the proportion of rare alleles in the sample compared to the 'wildSubset' for alleles with frequency <0.01.
        rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(samp > 0)))) / length(which(wildSubset < 0.01))
      } else {
        # Calculate the proportion of alleles in the sample that are also present in the 'wildSubset' for the case where j is not 1
        total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))]) / length(names(wildSubset))
        # Calculate the proportion of alleles of common alleles in the sample compared toe the 'wildSubset' where j is not 1
        common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.1))
        # Calculate the proportion of alleles of low frequency alleles in the sample compared to the 'wildSubset' for the case where j is not 1
        lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
        # Calculate the proportion of rare alleles in the sample compared ot the 'wildSubset' for the case where j is not 1
        rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset < 0.01))
      }
      # Combine the results (proportions) for each sample size into a matrix named 'categorymat_10loc'.
      categorymat <- cbind(total, common, lowfreq, rare)
      # extract the column names (categories) from the matrix
      category <- colnames(categorymat)
      # set row and column names for 'resamp_category10loc'
      dimnames(resamp_category) <- list(paste0("sample ", 1:nrow(categorymat)), category, paste0("replicate ", 1:num_reps))
      # Store the results for the current replicate in resamp_category5loc
      resamp_category[, , i] <- categorymat
    }
  }
  return(resamp_category)
}

## Linear model ##
# pass all arrays to a dataframe using the function. define the function that takes an input 'data_array'
analyze_resampling_array <- function(data_array) {
  # linear model of resampling array. Extract the column named 'total' from the array.
  # concatenate the extracted column values into the vector.
  totalsVector <- c(data_array[,1,]) 
  
  # Specify sample numbers column.
  # Create a vector of sample numbers from 1 to the number of rows in the 'total' column 
  gm_sampleNumbers <- 1:(nrow(data_array[,1,]))
  # Repeat the sample numbers vector for the number of replicates
  gm_sampleNumbers <- rep(gm_sampleNumbers, dim(data_array)[[3]])
  
  # Create data frame from resampling array values
  gm_DF <- data.frame(sampleNumbers=gm_sampleNumbers, totalValues=totalsVector)
  
  # Build and analyze linear models
  gm_Model <- lm(sampleNumbers ~ I((totalValues)^3), data = gm_DF)
  # Create a new data fram 'gm_newData with a single column 'totalValues' containing the value 0.95
  gm_newData <- data.frame(totalValues=0.95)
  # Use the linear model to predict the response for the new data frame. Specify 'interval = prediction to obtain a prediction interval
  gm_95MSSEprediction <- predict(gm_Model, gm_newData, interval = "prediction")
  
  # Pass the gm_95MSSEprediction to the object storing our results 
  # Store the predicted values and predictino interval in the object named 'result' 
  result <- gm_95MSSEprediction
  # Calculate the width of the prediction interval by substracting the lower limit from the upper limit
  piWidth <- gm_95MSSEprediction[3] - gm_95MSSEprediction[2]
  # Return a list containing the predicted values and the width of the prediction interval
  return(list(result = result, piWidth = piWidth))
}

## Filling in matrix ##
# Iterate through the arrays and store results in the matrix
# Initiate loop that iterates over the indices of 'array_list'
build_matrix_func <- function(array_list, input_matrix){
  for (i in 1:length(array_list)) {
    # Store results and piWidth values in the ith row of the matrix
    input_matrix[i, ] <- c(array_list[[i]]$result, 
                           array_list[[i]]$piWidth)
  }
  return(input_matrix)
}
# declare all objects for the loop
# this array stores a resampling array from a simulated genind object of a scenario, of which there are 5 simulated genind objects.
replicatesArray <- array(dim = c(1200, 4, 5))
# this array stores the 
scenarioArray <- array(dim = c(1200, 4, 0))
# this array 
combinedArray <- array(dim = c(1200, 4, 0))
MSAT_levels <- c(5, 10, 15, 20, 25)
QUAC_array_list = list(length(MSAT_levels))
# reading in the names of each of the scenarios to be processed
MSATscenarios <- list.files(path = 'Datasets/MSAT_N1200/', pattern = "genind.MSAT_", full.names = TRUE)
# a list that stores all the simulation scenarios after being added to the environment
MSATscenariosList = list(length(MSATscenarios))

for (i in 1:length(MSAT_levels)) {
  for (j in 1:length(MSATscenarios)) {
    # browser()
    # MSATscenariosList[[j]] <- readRDS(MSATscenarios[j])
    currentScenario <- readRDS(MSATscenarios[[j]])
    for (k in 1:length(currentScenario)){
      replicatesArray <- gm_resamp_array_function(currentScenario[[k]],MSAT_levels[[i]],5)
      scenarioArray <- abind(replicatesArray,scenarioArray)
    }
    combinedArray <- abind(scenarioArray,combinedArray)
    # empties out scenario array
    scenarioArray <- array(dim = c(1200, 4, 0))
    # QUAC_predict_results <- analyze_resampling_array(combinedArray)
    # print(QUAC_predict_results)
  }
  QUAC_array_list[[i]] <- combinedArray
  combinedArray <- array(dim = c(1200,4,0))
}