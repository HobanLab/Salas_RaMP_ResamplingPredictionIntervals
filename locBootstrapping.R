#####################################
# GDS 2023/12/06 Loci Bootstrapping #
#####################################
setwd("C:/Users/gsalas/Documents/resampling_CIs/Code/Datasets")
library(adegenet)
# This resampling array was created in the resampling walkthrough titled 
# resamplingwalkthru.R using the Q. Acerifolia MSAT genetic marker dataset. 
# resamp_category consists of 5 resampling replicates, each with unique allelic
# representation values based on random sampling amount of individuals across allele frequency categories.
# resamp_category was saved as a .rds file and uploaded to this R script to ensure repeatability.
resamp_category <- readRDS("resamp_category.RDS")

# One of the goals of this R script was to pass resampling arrays (of different loci amounts) through a 
# linear model and calculate the prediction interval values and the prediction 
# interval widths. we copied over the code from the predict.R script and made 
# the changes necessary to do this

# linear model of resamp_category.
totalsVector <- c(resamp_category[,"total",])
# Specify sample numbers column
gm_sampleNumbers <- 1:(nrow(resamp_category[,"total",]))
gm_sampleNumbers <- rep(gm_sampleNumbers, dim(resamp_category)[[3]])
# Create data.frame from resampling array values
gm_DF <- data.frame(sampleNumbers=gm_sampleNumbers, totalValues=totalsVector)
# Build and analyze linear models
gm_Model <- lm(sampleNumbers ~ I((totalValues)^3), data = gm_DF)
# ensure that the 95% value being predicted is typed as a proportion
gm_newData <- data.frame(totalValues=0.95)
gm_95MSSEprediction <- predict(gm_Model, gm_newData, interval = "prediction")
gm_95MSSEprediction
# calculate the prediction interval widths
piWidth <- gm_95MSSEprediction[3] - gm_95MSSEprediction[2]

# We used the Q. Acerifolia MSAT file to create resampling arrays of 
# different loci quantities. We followed all the same steps to read in a genind file
# from the resampling walkthrough, however, we create two genind objects that are
# separated by the garden and wild populations.

# Change the filepath below to the filepath for your particular system
QUAC.MSAT.filepath <- "C:/Users/gsalas/Documents/resampling_CIs/Code/Datasets/QUAC_woK_allpop_clean.gen"

# Number of characters used to code an allele. 
QUAC.MSAT.genind <- read.genepop(QUAC.MSAT.filepath, ncode = 3)

# "G" means garden, and the "W" means wild
QUAC.MSAT.genind@pop
# rename the populations to be "garden" or "wild" based on the G or W in the current string
# Correct garden popNames
levels(QUAC.MSAT.genind@pop)[grep(pattern = "QAc-G-", levels(QUAC.MSAT.genind@pop))] <-
  rep("garden", length(grep(pattern = "QAc-G-", levels(QUAC.MSAT.genind@pop))))
# Correct wild popNames
levels(QUAC.MSAT.genind@pop)[grep(pattern = "QAc-W-", levels(QUAC.MSAT.genind@pop))] <-
  rep("wild", length(grep(pattern = "QAc-W-", levels(QUAC.MSAT.genind@pop))))
# Recall the populations, to make sure they look better
QUAC.MSAT.genind@pop

# Garden and Wild Individuals
gardenRows <- which(QUAC.MSAT.genind@pop == "garden")
wildRows <- which(QUAC.MSAT.genind@pop == "wild")

# SEPARATE BY POPULATION
# Creating a genind object specifically isolating wild individual alleles will 
# ensure we are not including garden individuals when subsampling loci. 
populations <- seppop(QUAC.MSAT.genind)
QUAC.MSAT.WILD.genind <- populations$wild
QUAC.MSAT.GARDEN.genind <- populations$garden

# Below, we complete the process of loci bootstrapping. We begin with 5 loci, 
# then 10 loci, and the total loci amount. The code randomly samples sets of loci 25 times. The skeleton
# for creating the resampling array code is copied from the resampling walkthrough. 

###########
# 5 loci
# Create an empty 3D array with dimensions 164 rows, 4 columns and 25 slices, to store the results.
resamp_category5loc <- array(dim = c(164, 4, 25))
# Initiate a loop 25 times for 25 sets (5 loci in one set) of randomly selected loci.
for (i in 1:25) {
  # Randomly sample 5 loci names from the genind object
  samp_5loc <- sample(locNames(QUAC.MSAT.WILD.genind), size = 5, replace = FALSE)
  # Subset the genind object QUAC.MSAT.WILD.genind to include only the columns corresponding to the sampled loci.
  QUAC.MSAT.5loc.WILD.genind <- QUAC.MSAT.WILD.genind[, loc = samp_5loc]
  # objects
  # access the matrix that shows the type of alleles and quantity of said allele present in wild individuals
  wildSamp_5loc <- QUAC.MSAT.5loc.WILD.genind@tab
  # Calculate the sum of each column in the matrix 'wildSamp_5loc', ignoring NA values. Identify the indices
  # where the sum is not equal to zero, this indicates columns with variation in allele counts. Subset the 
  # original matrix 'wildSamp_5loc' by selecting only the columns identified in the previous step, removing
  # columns with no variation in allele counts
  wildSamp_5loc <- wildSamp_5loc[, which(colSums(wildSamp_5loc, na.rm = TRUE) != 0)]
  # Calculating a vector of allele frequencies (sum of allele counts divided by number of haplotypes, or individuals * 2)
  wildComplete <- colSums(wildSamp_5loc, na.rm = TRUE) / (nrow(wildSamp_5loc) * 2)
  # Subset 'wildComplete' to include only the alleles with non-zero frequency.
  wildSubset <- wildComplete[wildComplete > 0]
  # Initialize vectors to store results
  total <- vector(length = nrow(wildSamp_5loc))
  common <- vector(length = nrow(wildSamp_5loc))
  lowfreq <- vector(length = nrow(wildSamp_5loc))
  rare <- vector(length = nrow(wildSamp_5loc))
  # Loop for each sample size, ranging from 1 to the number of loci in the wild population
  for (j in 1:nrow(wildSamp_5loc)) {
    # Randomly select a subset of rows from the matrix, without replacement
    samp <- sample(nrow(wildSamp_5loc), size = j, replace = FALSE)
    # Subset the original matrix 'wildSamp_5loc' to include only the rows randomly selected in the previous step.
    samp <- wildSamp_5loc[samp,]
    # Measure the proportion of allelic representation in that sample
    # Check if it's the first iteration of the inner loop
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
    # Combine the results (proportions) for each sample size into a matrix named 'categorymat_5loc'.
    categorymat_5loc <- cbind(total, common, lowfreq, rare)
    # extract the column names (categories) from the matrix
    category <- colnames(categorymat_5loc)
    # set row and column names for 'resamp_category5loc'
    dimnames(resamp_category5loc) <- list(paste0("sample ", 1:nrow(categorymat_5loc)), category, paste0("replicate ", 1:25))
    # Store the results for the current replicate in resamp_category5loc
    resamp_category5loc[, , i] <- categorymat_5loc
  }
}
# array containing the results
resamp_category5loc

############
# 10 loci
resamp_category10loc <- array(dim = c(164, 4, 25))
# Loop 25 times for 25 sets (5 loci in one set) of randomly selected loci
for (i in 1:25) {
  # Randomly sample 5 loci names from the genind object
  samp_10loc <- sample(locNames(QUAC.MSAT.WILD.genind), size = 10, replace = FALSE)
  # Subset the genind object QUAC.MSAT.WILD.genind to include only the columns corresponding to the sampled loci.
  QUAC.MSAT.10loc.WILD.genind <- QUAC.MSAT.WILD.genind[, loc = samp_10loc]
  # objects
  # access the matrix that shows the type of alleles and quantity of said allele present in wild individuals
  wildSamp_10loc <- QUAC.MSAT.10loc.WILD.genind@tab
  # Calculate the sum of each column in the matrix 'wildSamp_10loc', ignoring NA values. Identify the indices
  # where the sum is not equal to zero, this indicates columns with variation in allele counts. Subset the 
  # original matrix 'wildSamp_10loc' by selecting only the columns identified in the previous step, removing
  # columns with no variation in allele counts
  wildSamp_10loc <- wildSamp_10loc[, which(colSums(wildSamp_10loc, na.rm = TRUE) != 0)]
  # Calculating a vector of allele frequencies (sum of allele counts divided by number of haplotypes, or individuals * 2)
  wildComplete <- colSums(wildSamp_10loc, na.rm = TRUE) / (nrow(wildSamp_10loc) * 2)
  # Subset 'wildComplete' to include only the alleles with non-zero frequency.
  wildSubset <- wildComplete[wildComplete > 0]
  # Initialize vectors to store results
  total <- vector(length = nrow(wildSamp_10loc))
  common <- vector(length = nrow(wildSamp_10loc))
  lowfreq <- vector(length = nrow(wildSamp_10loc))
  rare <- vector(length = nrow(wildSamp_10loc))
  # Loop for each sample size, ranging from 1 to the number of loci in the wild population
  for (j in 1:nrow(wildSamp_10loc)) {
    # Randomly select a subset of rows from the matrix, without replacement
    samp <- sample(nrow(wildSamp_10loc), size = j, replace = FALSE)
    # Subset the original matrix 'wildSamp_5loc' to include only the rows randomly selected in the previous step.
    samp <- wildSamp_10loc[samp,]
    # Measure the proportion of allelic representation in that sample
    # Check if it's the first iteration of the inner loop
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
    categorymat_10loc <- cbind(total, common, lowfreq, rare)
    # extract the column names (categories) from the matrix
    category <- colnames(categorymat_10loc)
    # set row and column names for 'resamp_category10loc'
    dimnames(resamp_category10loc) <- list(paste0("sample ", 1:nrow(categorymat_10loc)), category, paste0("replicate ", 1:25))
    # Store the results for the current replicate in resamp_category10loc
    resamp_category10loc[, , i] <- categorymat_10loc
  }
}
# array containing the results
resamp_category10loc[,,]

###############
# total loci
resamp_categorytotloc <- array(dim = c(164, 4, 25))
# Loop 25 times for 25 sets (5 loci in one set) of randomly selected loci
for (i in 1:25) {
  # Randomly sample 5 loci names
  samp_totloc <- sample(locNames(QUAC.MSAT.WILD.genind), size = 15, replace = FALSE)
  # Subset the genind object QUAC.MSAT.WILD.genind to include only the columns corresponding to the sampled loci.
  QUAC.MSAT.totloc.WILD.genind <- QUAC.MSAT.WILD.genind[, loc = samp_totloc]
  # objects
  # access the matrix that shows the type of alleles and quantity of said allele present in wild individuals
  wildSamp_totloc <- QUAC.MSAT.totloc.WILD.genind@tab
  # Calculate the sum of each column in the matrix 'wildSamp_10loc', ignoring NA values. Identify the indices
  # where the sum is not equal to zero, this indicates columns with variation in allele counts. Subset the 
  # original matrix 'wildSamp_10loc' by selecting only the columns identified in the previous step, removing
  # columns with no variation in allele counts
  wildSamp_totloc <- wildSamp_totloc[, which(colSums(wildSamp_totloc, na.rm = TRUE) != 0)]
  # Calculating a vector of allele frequencies (sum of allele counts divided by number of haplotypes, or individuals * 2)
  wildComplete <- colSums(wildSamp_totloc, na.rm = TRUE) / (nrow(wildSamp_totloc) * 2)
  # Subset 'wildComplete' to include only the alleles with non-zero frequency.
  wildSubset <- wildComplete[wildComplete > 0]
  # Initialize vectors to store results
  total <- vector(length = nrow(wildSamp_totloc))
  common <- vector(length = nrow(wildSamp_totloc))
  lowfreq <- vector(length = nrow(wildSamp_totloc))
  rare <- vector(length = nrow(wildSamp_totloc))
  # Loop for each sample size, ranging from 1 to the number of loci in the wild population
  for (j in 1:nrow(wildSamp_totloc)) {
    # Randomly select a subset of rows from the matrix, without replacement
    samp <- sample(nrow(wildSamp_totloc), size = j, replace = FALSE)
    # Subset the original matrix 'wildSamp_totloc' to include only the rows randomly selected in the previous step.
    samp <- wildSamp_totloc[samp,]
    # Measure the proportion of allelic representation in that sample
    # Check if it's the first iteration of the inner loop
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
    # Combine the results (proportions) for each sample size into a matrix named 'categorymat_totloc'.
    categorymat_totloc <- cbind(total, common, lowfreq, rare)
    # extract the column names (categories) from the matrix
    category <- colnames(categorymat_totloc)
    # set row and column names for 'resamp_categorytotloc'
    dimnames(resamp_categorytotloc) <- list(paste0("sample ", 1:nrow(categorymat_totloc)), category, paste0("replicate ", 1:25))
    # Store the results for the current replicate in resamp_categorytotloc
    resamp_categorytotloc[, , i] <- categorymat_totloc
  }
}
# array containing the results
resamp_categorytotloc

# Here we have all the resampling arrays with the different amounts of loci
resamp_category5loc
resamp_category10loc
resamp_categorytotloc

# We created a function that calculates the prediction intervals. That function 
# was assigned the object name analyze_resampling_array. the input is the resampling array name.

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

# this is a list of the resampling arrays which will be iterated in a loop to 
# execute the analyze_resampling_array function which will output the prediction interval
# values and prediction interval widths. 
array_list <- list(resamp_categorytotloc, resamp_category10loc, resamp_category5loc)

# this matrix will store the pi values and pi widths
# Create an empty matrix to store the results
results_matrix <- matrix(nrow = length(array_list), ncol = 4)
# Set column names for 'results_Matrix'
colnames(results_matrix) <- c("fit", "lower", "upper", "piWidth")
# Set row names for 'results_matrix'
rownames(results_matrix) <- c("resamp_categorytotloc","resamp_category10loc", "resamp_category5loc")

# Iterate through the arrays and store results in the matrix
# Initiate loop that iterates over the idicies of 'array_list'
for (i in 1:length(array_list)) {
  # Call the 'analyze_resampling_array' function on the ith element of 'array_list'. This function returns a list with 
  # resul and piWidth values.
  analyze_resampling_array(array_list[[i]])
  # Store results and piWidth values in the ith row of the matrix
  results_matrix[i, ] <- c(analyze_resampling_array(array_list[[i]])$result, 
                           analyze_resampling_array(array_list[[i]])$piWidth)
  
}
# Print final 'results_matrix' showing the fit, lower, upper and piWidth values for each resampling array
print(results_matrix)

# the results are saved to a .csv file
write.csv(results_matrix, file = "C:/Users/gsalas/Documents/resampling_CIs/Code/resamp_lociMatrix.csv", 
          row.names = TRUE)