################################################
# GDS 2024/01/10 QUAC QUBO Loci Bootstrapping #
################################################
setwd("C:/Users/gsalas/Documents/resampling_CIs/Code/Datasets")
# get a list of file names in the directory contating "Wild" in their names and return the full file paths
datasets <- list.files(path = "LociBootstrapping_Datasets", pattern = "Wild", full.names = TRUE)
# Read RDS files and store them in a list named 'gm_list'
gm_list <- list(
  readRDS(datasets[1]),
  readRDS(datasets[2]),
  readRDS(datasets[3]),
  readRDS(datasets[4])
)
# Assign each element of 'gm_list' to a seperate object representing different genetic datasets
QUAC.MSAT.Wild.genind <- gm_list[[1]]
QUAC.SNP.Wild.genind <- gm_list[[2]]
QUBO.MSAT.Wild.genind <- gm_list[[3]]
QUBO.SNP.Wild.genind <- gm_list[[4]]
library(adegenet)

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
  totalsVector <- c(data_array[,"total",]) 
  
  # Specify sample numbers column.
  # Create a vector of sample numbers from 1 to the number of rows in the 'total' column 
  gm_sampleNumbers <- 1:(nrow(data_array[,"total",]))
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

## QUAC Dataset ##
QUAC_loci_ranges <- c(3:15,50,100,200,500,1000,5000,10000,nLoc(QUAC.SNP.Wild.genind))
# Create a list equal to the length of the different levels that the arrays will be stored into
QUAC_array_list = list(length(QUAC_loci_ranges))
# Create a list equal to the length of the different levels that the predict values will be stored into
QUAC_predict_results <- list(length(QUAC_loci_ranges))
# this matrix will store the pi values and pi widths
# Create an empty matrix to store the results
QUAC_results_matrix <- matrix(nrow = length(QUAC_loci_ranges), ncol = 4)
# Set column names for 'results_Matrix'
colnames(QUAC_results_matrix) <- c("fit", "lower", "upper", "piWidth")
# Set row names for 'results_matrix'
# declare an object to store the character values of the row names
loci_names <- vector(length = nrow(QUAC_results_matrix))
# Create a loop that iterates through the length of the different levels
for (i in 1:length(QUAC_loci_ranges)) {
  # declare a variable that is going to take the value of the level of loci you want to subset based on the list
  loci_amount <- QUAC_loci_ranges[i]
  # if else statement that will store the arrays of each set of loci into the list from the MSAT genind,
  # once the string reaches to 50, the else statement will iterate through the SNP genind.
  # the gm_resamp_array_function takes three arguments, the genind object, the loci amount, and the amount of replicates.
  if (loci_amount < 50) {
    QUAC_array_list[[i]] <- gm_resamp_array_function(QUAC.MSAT.Wild.genind, loci_amount, 25)
  } else{
    QUAC_array_list[[i]] <- gm_resamp_array_function(QUAC.SNP.Wild.genind, loci_amount, 25)
  }
  QUAC_predict_results[[i]] <- analyze_resampling_array(QUAC_array_list[[i]])
  print(QUAC_predict_results)
  # if else statement that stores the output of the character results for MSAT and SNP levels
  if (QUAC_loci_ranges[i] < 50) {
    loci_names[i] <- paste0("QUAC_MSAT ",QUAC_loci_ranges[i]," loci")
  } else{
    loci_names[i] <- paste0("QUAC_SNP ",QUAC_loci_ranges[i]," loci")
  }
  # fill the row names with the loci_names
  rownames(QUAC_results_matrix) <- loci_names
  # build the matrix using the predict results and the empty results matrix
  QUAC_final_matrix_results <- build_matrix_func(QUAC_predict_results, QUAC_results_matrix)
}
print(QUAC_final_matrix_results)
write.csv(QUAC_final_matrix_results, 
          file = "C:/Users/gsalas/Documents/resampling_CIs/Code/Outputs/QUAC_resamp_loci.csv",
          row.names = TRUE)

## QUBO dataset ##
QUBO_loci_ranges <- c(3:9, 50, 100, 200, 500, 1000, 2000, 4000, nLoc(QUBO.SNP.Wild.genind))
# Create a list equal to the length of the different levels that the arrays will be stored into
QUBO_array_list = list(length(QUBO_loci_ranges))
# Create a list equal to the length of the different levels that the predict values will be stored into
QUBO_predict_results <- list(length(QUBO_loci_ranges))
# this matrix will store the pi values and pi widths
# Create an empty matrix to store the results
QUBO_results_matrix <- matrix(nrow = length(QUBO_loci_ranges), ncol = 4)
# Set column names for 'results_Matrix'
colnames(QUBO_results_matrix) <- c("fit", "lower", "upper", "piWidth")
# Set row names for 'results_matrix'
# declare an object to store the character values of the row names
loci_names <- vector(length = nrow(QUBO_results_matrix))
# Create a loop that iterates through the length of the different levels
for (i in 1:length(QUBO_loci_ranges)) {
  # declare a variable that is going to take the value of the level of loci you want to subset based on the list
  loci_amount <- QUBO_loci_ranges[i]
  # if else statement that will store the arrays of each set of loci into the list from the MSAT genind,
  # once the string reaches to 50, the else statement will iterate through the SNP genind.
  # the gm_resamp_array_function takes three arguments, the genind object, the loci amount, and the amount of replicates.
  if (loci_amount < 50) {
    QUBO_array_list[[i]] <- gm_resamp_array_function(QUBO.MSAT.Wild.genind, loci_amount, 25)
  } else{
    QUBO_array_list[[i]] <- gm_resamp_array_function(QUBO.SNP.Wild.genind, loci_amount, 25)
  }
  # store the results into a list 
  QUBO_predict_results[[i]] <- analyze_resampling_array(QUBO_array_list[[i]])
  print(QUBO_predict_results)
  # if else statement that stores the output of the character results for MSAT and SNP levels
  if (QUBO_loci_ranges[i] < 50) {
    loci_names[i] <- paste0("QUBO_MSAT ",QUBO_loci_ranges[i]," loci")
  } else{
    loci_names[i] <- paste0("QUBO_SNP ",QUBO_loci_ranges[i]," loci")
  }
  # fill the row names with the loci_names
  rownames(QUBO_results_matrix) <- loci_names
  # build the matrix using the predict results and the empty results matrix
  QUBO_final_matrix_results <- build_matrix_func(QUBO_predict_results, QUBO_results_matrix)
}
print(QUBO_final_matrix_results)
write.csv(QUBO_final_matrix_results,
          file = "C:/Users/gsalas/Documents/resampling_CIs/Code/Outputs/QUBO_resamp_loci.csv",
          row.names = TRUE)