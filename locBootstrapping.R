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
# Creating a genind object specifically isolating wild indivudal alleles will 
# ensure we are not including garden individuals when subsampling loci. 
populations <- seppop(QUAC.MSAT.genind)
QUAC.MSAT.WILD.genind <- populations$wild
QUAC.MSAT.GARDEN.genind <- populations$garden

# Below, we complete the process of loci bootstrapping. We begin with 5 loci, 
# then 10 loci, and the total loci amount. The code randomly samples sets of loci 25 times. The skeleton
# for creating the resampling array code is copied from the resampling walkthrough. 

###########
# 5 loci
resamp_category5loc <- array(dim = c(164, 4, 25))
# Loop 25 times for 25 sets (5 loci in one set) of randomly selected loci
for (i in 1:25) {
  # Randomly sample 5 loci names
  samp_5loc <- sample(locNames(QUAC.MSAT.WILD.genind), size = 5, replace = FALSE)
  # Subset the genind object based on the sampled loci names
  QUAC.MSAT.5loc.WILD.genind <- QUAC.MSAT.WILD.genind[, loc = samp_5loc]
  # objects
  wildSamp_5loc <- QUAC.MSAT.5loc.WILD.genind@tab
  wildSamp_5loc <- wildSamp_5loc[, which(colSums(wildSamp_5loc, na.rm = TRUE) != 0)]
  wildComplete <- colSums(wildSamp_5loc, na.rm = TRUE) / (nrow(wildSamp_5loc) * 2)
  wildSubset <- wildComplete[wildComplete > 0]
  total <- vector(length = nrow(wildSamp_5loc))
  common <- vector(length = nrow(wildSamp_5loc))
  lowfreq <- vector(length = nrow(wildSamp_5loc))
  rare <- vector(length = nrow(wildSamp_5loc))
  # Loop for each sample size
  for (j in 1:nrow(wildSamp_5loc)) {
    samp <- sample(nrow(wildSamp_5loc), size = j, replace = FALSE)
    samp <- wildSamp_5loc[samp,]
    # Measure the proportion of allelic representation in that sample
    if (j == 1) {
      total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(samp > 0)))]) / length(names(wildSubset))
      common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.1))
      lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
      rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(samp > 0)))) / length(which(wildSubset < 0.01))
    } else {
      total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))]) / length(names(wildSubset))
      common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.1))
      lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
      rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset < 0.01))
    }
    categorymat_5loc <- cbind(total, common, lowfreq, rare)
    category <- colnames(categorymat_5loc)
    dimnames(resamp_category5loc) <- list(paste0("sample ", 1:nrow(categorymat_5loc)), category, paste0("replicate ", 1:25))
    # Store the results for the current replicate in resamp_category5loc
    resamp_category5loc[, , i] <- categorymat_5loc
  }
}
resamp_category5loc

############
# 10 loci
resamp_category10loc <- array(dim = c(164, 4, 25))
# Loop 25 times for 25 sets (5 loci in one set) of randomly selected loci
for (i in 1:25) {
  # Randomly sample 5 loci names
  samp_10loc <- sample(locNames(QUAC.MSAT.WILD.genind), size = 10, replace = FALSE)
  
  # Subset the genind object based on the sampled loci names
  QUAC.MSAT.10loc.WILD.genind <- QUAC.MSAT.WILD.genind[, loc = samp_10loc]
  # objects
  wildSamp_10loc <- QUAC.MSAT.10loc.WILD.genind@tab
  wildSamp_10loc <- wildSamp_10loc[, which(colSums(wildSamp_10loc, na.rm = TRUE) != 0)]
  wildComplete <- colSums(wildSamp_10loc, na.rm = TRUE) / (nrow(wildSamp_10loc) * 2)
  wildSubset <- wildComplete[wildComplete > 0]
  total <- vector(length = nrow(wildSamp_10loc))
  common <- vector(length = nrow(wildSamp_10loc))
  lowfreq <- vector(length = nrow(wildSamp_10loc))
  rare <- vector(length = nrow(wildSamp_10loc))
  # Loop for each sample size
  for (j in 1:nrow(wildSamp_10loc)) {
    samp <- sample(nrow(wildSamp_10loc), size = j, replace = FALSE)
    samp <- wildSamp_10loc[samp,]
    # Measure the proportion of allelic representation in that sample
    if (j == 1) {
      total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(samp > 0)))]) / length(names(wildSubset))
      common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.1))
      lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
      rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(samp > 0)))) / length(which(wildSubset < 0.01))
    } else {
      total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))]) / length(names(wildSubset))
      common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.1))
      lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
      rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset < 0.01))
    }
    categorymat_10loc <- cbind(total, common, lowfreq, rare)
    category <- colnames(categorymat_10loc)
    dimnames(resamp_category10loc) <- list(paste0("sample ", 1:nrow(categorymat_10loc)), category, paste0("replicate ", 1:25))
    # Store the results for the current replicate in resamp_category10loc
    resamp_category10loc[, , i] <- categorymat_10loc
  }
}
resamp_category10loc[,,]

###############
# total loci
resamp_categorytotloc <- array(dim = c(164, 4, 25))
# Loop 25 times for 25 sets (5 loci in one set) of randomly selected loci
for (i in 1:25) {
  # Randomly sample 5 loci names
  samp_totloc <- sample(locNames(QUAC.MSAT.WILD.genind), size = 15, replace = FALSE)
  
  # Subset the genind object based on the sampled loci names
  QUAC.MSAT.totloc.WILD.genind <- QUAC.MSAT.WILD.genind[, loc = samp_totloc]
  
  # objects
  wildSamp_totloc <- QUAC.MSAT.totloc.WILD.genind@tab
  wildSamp_totloc <- wildSamp_totloc[, which(colSums(wildSamp_totloc, na.rm = TRUE) != 0)]
  wildComplete <- colSums(wildSamp_totloc, na.rm = TRUE) / (nrow(wildSamp_totloc) * 2)
  wildSubset <- wildComplete[wildComplete > 0]
  total <- vector(length = nrow(wildSamp_totloc))
  common <- vector(length = nrow(wildSamp_totloc))
  lowfreq <- vector(length = nrow(wildSamp_totloc))
  rare <- vector(length = nrow(wildSamp_totloc))
  # Loop for each sample size
  for (j in 1:nrow(wildSamp_totloc)) {
    samp <- sample(nrow(wildSamp_totloc), size = j, replace = FALSE)
    samp <- wildSamp_totloc[samp,]
    # Measure the proportion of allelic representation in that sample
    if (j == 1) {
      total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(samp > 0)))]) / length(names(wildSubset))
      common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.1))
      lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
      rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(samp > 0)))) / length(which(wildSubset < 0.01))
    } else {
      total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))]) / length(names(wildSubset))
      common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.1))
      lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
      rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset < 0.01))
    }
    categorymat_totloc <- cbind(total, common, lowfreq, rare)
    category <- colnames(categorymat_totloc)
    dimnames(resamp_categorytotloc) <- list(paste0("sample ", 1:nrow(categorymat_totloc)), category, paste0("replicate ", 1:25))
    # Store the results for the current replicate in resamp_category10loc
    resamp_categorytotloc[, , i] <- categorymat_totloc
  }
}
resamp_categorytotloc

# Here we have all the resampling arrays with the different amounts of loci
resamp_category5loc
resamp_category10loc
resamp_categorytotloc

# We created a function that calculates the prediction intervals. That function 
# was assigned the object name analyze_resampling_array. the input is the resampling array name.

# pass all arrays to a dataframe using the function. 
analyze_resampling_array <- function(data_array) {
  # linear model of resampling array
  totalsVector <- c(data_array[,"total",]) 
  
  # Specify sample numbers column
  gm_sampleNumbers <- 1:(nrow(data_array[,"total",]))
  gm_sampleNumbers <- rep(gm_sampleNumbers, dim(data_array)[[3]])
  
  # Create data frame from resampling array values
  gm_DF <- data.frame(sampleNumbers=gm_sampleNumbers, totalValues=totalsVector)
  
  # Build and analyze linear models
  gm_Model <- lm(sampleNumbers ~ I((totalValues)^3), data = gm_DF)
  gm_newData <- data.frame(totalValues=0.95)
  gm_95MSSEprediction <- predict(gm_Model, gm_newData, interval = "prediction")
  
  # Pass the gm_95MSSEprediction to the object storing our results 
  result <- gm_95MSSEprediction
  
  piWidth <- gm_95MSSEprediction[3] - gm_95MSSEprediction[2]
  
  return(list(result = result, piWidth = piWidth))
}

# this is a list of the resampling arrays which will be iterated in a loop to 
# execute the analyze_resampling_array function which will output the prediction interval
# values and prediction interval widths. 
array_list <- list(resamp_categorytotloc, resamp_category10loc, resamp_category5loc)

# this matrix will store the pi values and pi widths
# Create an empty matrix to store the results
results_matrix <- matrix(nrow = length(array_list), ncol = 4)
colnames(results_matrix) <- c("fit", "lower", "upper", "piWidth")
rownames(results_matrix) <- c("resamp_categorytotloc","resamp_category10loc", "resamp_category5loc")

# Iterate through the arrays and store results in the matrix
for (i in 1:length(array_list)) {
  analyze_resampling_array(array_list[[i]])
  # Store results in the matrix
  results_matrix[i, ] <- c(analyze_resampling_array(array_list[[i]])$result, 
                           analyze_resampling_array(array_list[[i]])$piWidth)
  
}
print(results_matrix)

# the results are saved to a .csv file
write.csv(results_matrix, file = "C:/Users/gsalas/Documents/resampling_CIs/Code/resamp_lociMatrix.csv", 
          row.names = TRUE)