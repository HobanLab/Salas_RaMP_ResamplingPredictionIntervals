#####################################
# GDS 2023/12/06 Loci Bootstrapping #
#####################################
setwd("C:/Users/gsalas/Documents/resampling_CIs/Code/Datasets")
library(adegenet)

resamp_category <- readRDS("resamp_category.RDS")

# linear model of resampling array
totalsVector <- c(resamp_category[,"total",])
# Specify sample numbers column
gm_sampleNumbers <- 1:(nrow(resamp_category[,"total",]))
gm_sampleNumbers <- rep(gm_sampleNumbers, dim(resamp_category)[[3]])
# Create data.frame from resampling array values
gm_DF <- data.frame(sampleNumbers=gm_sampleNumbers, totalValues=totalsVector)
# Build and analyze linear models
gm_Model <- lm(sampleNumbers ~ I((totalValues)^3), data = gm_DF)
gm_newData <- data.frame(totalValues=0.95)
gm_95MSSEprediction <- predict(gm_Model, gm_newData, interval = "prediction")
# Pass the gm_95MSSEprediction to the object storing our results by iterating
# storing them in the rows index of predict_matrix
gm_95MSSEprediction

ciWidth <- gm_95MSSEprediction[3] - gm_95MSSEprediction[2]

## READING IN THE GENIND FILE ## 
# Change the filepath below to the filepath for your particular system
QUAC.MSAT.filepath <- "C:/Users/gsalas/Documents/resampling_CIs/Code/Datasets/QUAC_woK_allpop_clean.gen"

# QUESTION 1: What does the ncode argument signify? How can you find out?
# Number of characters used to code an allele. Information on ncode can be found on the help tab by typing ?read.genepop
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
populations <- seppop(QUAC.MSAT.genind)
QUAC.MSAT.WILD.genind <- populations$wild
QUAC.MSAT.GARDEN.genind <- populations$garden

# Randomly sample a subset of loci
loc_samp <- sample(locNames(QUAC.MSAT.WILD.genind), size = 3, replace = FALSE)
QUAC.MSAT.3loc.WILD.genind <- QUAC.MSAT.WILD.genind[,loc=loc_samp]
locNames(QUAC.MSAT.3loc.WILD.genind)
wildSamp_3loc <- QUAC.MSAT.3loc.WILD.genind@tab

##########################################
samp_10loc <- sample(locNames(QUAC.MSAT.WILD.genind), size = 10, replace = FALSE)
QUAC.MSAT.10loc.WILD.genind <- QUAC.MSAT.WILD.genind[,loc=samp_10loc]
locNames(QUAC.MSAT.10loc.WILD.genind)
wildSamp_10loc <- QUAC.MSAT.10loc.WILD.genind@tab
# ALL ALLELE CATEGORIES, ALL SAMPLE SIZES #
# QUESTION 31:
# 
wildComplete <- colSums(wildSamp_10loc, na.rm = TRUE)/(nrow(wildSamp_10loc)*2) 
wildSubset <- wildComplete[wildComplete > 0]
total <- vector(length = nrow(wildSamp_10loc))
common <- vector(length = nrow(wildSamp_10loc))
lowfreq <- vector(length = nrow(wildSamp_10loc))
rare <- vector(length = nrow(wildSamp_10loc))
wildSamp_10loc <- wildSamp_10loc[,which(colSums(wildSamp_10loc, na.rm = TRUE)!= 0)]

for(i in 1:nrow(wildSamp_10loc)){
  # browser()
  # Use sample to randomly subsample the matrix of wild individuals
  samp <- sample(nrow(wildSamp_10loc), size = i, replace = FALSE)
  samp <- wildSamp_10loc[samp,]
  # samp <- sample(wildSamples[!is.na(wildSamples)], size = i, replace = FALSE)
  # Now, measure the proportion of allelic representation in that sample
  if (i==1) {
    total[i] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(samp > 0)))])/length(names(wildSubset))
    common[i] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(samp > 0))))/length(which(wildSubset > 0.1))
    lowfreq[i] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(samp > 0))))/length(which(wildSubset > 0.01 & wildSubset < 0.1))
    rare[i] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(samp > 0))))/length(which(wildSubset < 0.01))
  } else{
    total[i] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(colSums(samp, na.rm=TRUE) > 0)))])/length(names(wildSubset))
    common[i] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(colSums(samp, na.rm=TRUE) > 0))))/length(which(wildSubset > 0.1))
    lowfreq[i] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(colSums(samp, na.rm=TRUE) > 0))))/length(which(wildSubset > 0.01 & wildSubset < 0.1))
    rare[i] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(colSums(samp, na.rm=TRUE) > 0))))/length(which(wildSubset < 0.01))
    
  }
  categorymat_10loc <- cbind(total,common,lowfreq,rare)
}
print(categorymat_10loc)

# Resampling in replicate: all allele categories, all sample sizes #
# QUESTION 32: Repeat Question 31, but store the values across 5 resampling replicates
resamp_category10loc <- array(dim = c(164,4,5))
category <- colnames(categorymat_10loc)
dimnames(resamp_category10loc) <- list(paste0("sample ", 1:nrow(categorymat_10loc)), category, paste0("replicate ",1:5))
j <- length((dimnames(resamp_category10loc)[[3]]))
for (q in 1:j){
  for(i in 1:nrow(wildSamp_10loc)){
    # browser()
    samp <- sample(nrow(wildSamp_10loc), size = i, replace = FALSE)
    samp <- wildSamp_10loc[samp,]
    # Use sample to randomly subsample the matrix of wild individuals
    # Now, measure the proportion of allelic representation in that sample
    if (i==1) {
      total[i] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(samp > 0)))])/length(names(wildSubset))
      common[i] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(samp > 0))))/length(which(wildSubset > 0.1))
      lowfreq[i] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(samp > 0))))/length(which(wildSubset > 0.01 & wildSubset < 0.1))
      rare[i] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(samp > 0))))/length(which(wildSubset < 0.01))
    } else{
      total[i] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(colSums(samp, na.rm=TRUE) > 0)))])/length(names(wildSubset))
      common[i] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(colSums(samp, na.rm=TRUE) > 0))))/length(which(wildSubset > 0.1))
      lowfreq[i] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(colSums(samp, na.rm=TRUE) > 0))))/length(which(wildSubset > 0.01 & wildSubset < 0.1))
      rare[i] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(colSums(samp, na.rm=TRUE) > 0))))/length(which(wildSubset < 0.01))
    }
  }
  categorymat_10loc <- cbind(total,common,lowfreq,rare)
  resamp_category10loc[,,q] <- categorymat_10loc
}
resamp_category10loc

########################################

samp_5loc <- sample(locNames(QUAC.MSAT.WILD.genind), size = 5, replace = FALSE)
QUAC.MSAT.5loc.WILD.genind <- QUAC.MSAT.WILD.genind[,loc=samp_5loc]
locNames(QUAC.MSAT.5loc.WILD.genind)
wildSamp_5loc <- QUAC.MSAT.5loc.WILD.genind@tab
# ALL ALLELE CATEGORIES, ALL SAMPLE SIZES #
# QUESTION 31:
# 
wildComplete <- colSums(wildSamp_5loc, na.rm = TRUE)/(nrow(wildSamp_5loc)*2) 
wildSubset <- wildComplete[wildComplete > 0]
total <- vector(length = nrow(wildSamp_5loc))
common <- vector(length = nrow(wildSamp_5loc))
lowfreq <- vector(length = nrow(wildSamp_5loc))
rare <- vector(length = nrow(wildSamp_5loc))
wildSamp_5loc <- wildSamp_5loc[,which(colSums(wildSamp_5loc, na.rm = TRUE)!= 0)]

for(i in 1:nrow(wildSamp_5loc)){
  # browser()
  # Use sample to randomly subsample the matrix of wild individuals
  samp <- sample(nrow(wildSamp_5loc), size = i, replace = FALSE)
  samp <- wildSamp_5loc[samp,]
  # samp <- sample(wildSamples[!is.na(wildSamples)], size = i, replace = FALSE)
  # Now, measure the proportion of allelic representation in that sample
  if (i==1) {
    total[i] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(samp > 0)))])/length(names(wildSubset))
    common[i] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(samp > 0))))/length(which(wildSubset > 0.1))
    lowfreq[i] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(samp > 0))))/length(which(wildSubset > 0.01 & wildSubset < 0.1))
    rare[i] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(samp > 0))))/length(which(wildSubset < 0.01))
  } else{
    total[i] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(colSums(samp, na.rm=TRUE) > 0)))])/length(names(wildSubset))
    common[i] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(colSums(samp, na.rm=TRUE) > 0))))/length(which(wildSubset > 0.1))
    lowfreq[i] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(colSums(samp, na.rm=TRUE) > 0))))/length(which(wildSubset > 0.01 & wildSubset < 0.1))
    rare[i] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(colSums(samp, na.rm=TRUE) > 0))))/length(which(wildSubset < 0.01))
    
  }
  categorymat_5loc <- cbind(total,common,lowfreq,rare)
}
print(categorymat_5loc)

# Resampling in replicate: all allele categories, all sample sizes #
# QUESTION 32: Repeat Question 31, but store the values across 5 resampling replicates
resamp_category5loc <- array(dim = c(164,4,5))
category <- colnames(categorymat_5loc)
dimnames(resamp_category5loc) <- list(paste0("sample ", 1:nrow(categorymat_5loc)), category, paste0("replicate ",1:5))
j <- length((dimnames(resamp_category5loc)[[3]]))
for (q in 1:j){
  for(i in 1:nrow(wildSamp_5loc)){
    # browser()
    samp <- sample(nrow(wildSamp_5loc), size = i, replace = FALSE)
    samp <- wildSamp_5loc[samp,]
    # Use sample to randomly subsample the matrix of wild individuals
    # Now, measure the proportion of allelic representation in that sample
    if (i==1) {
      total[i] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(samp > 0)))])/length(names(wildSubset))
      common[i] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(samp > 0))))/length(which(wildSubset > 0.1))
      lowfreq[i] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(samp > 0))))/length(which(wildSubset > 0.01 & wildSubset < 0.1))
      rare[i] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(samp > 0))))/length(which(wildSubset < 0.01))
    } else{
      total[i] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(colSums(samp, na.rm=TRUE) > 0)))])/length(names(wildSubset))
      common[i] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(colSums(samp, na.rm=TRUE) > 0))))/length(which(wildSubset > 0.1))
      lowfreq[i] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(colSums(samp, na.rm=TRUE) > 0))))/length(which(wildSubset > 0.01 & wildSubset < 0.1))
      rare[i] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(colSums(samp, na.rm=TRUE) > 0))))/length(which(wildSubset < 0.01))
    }
  }
  categorymat_5loc <- cbind(total,common,lowfreq,rare)
  resamp_category5loc[,,q] <- categorymat_5loc
}
resamp_category5loc


# three resampling arrays
resamp_category
resamp_category10loc
resamp_category5loc

# test <- vector(length = 25)
# for (i in 1:5) {
#   test[i] <- sample(locNames(QUAC.MSAT.WILD.genind), size = 5, replace = FALSE)
# }
# test
# test <-  function(DATA, x, y){}
# for (i in 1:()) {
#   output[i] <- test(data, x, y)
# }

##########
# 5 sets of loci resampled 5 times and results of calculations stored into a slot the resampling array
genind_list <- list()

# Loop 5 times
for (i in 1:5) {
  # Randomly sample 5 loci names
  samp_5loc <- sample(locNames(QUAC.MSAT.WILD.genind), size = 5, replace = FALSE)
  
  # Subset the genind object based on the sampled loci names
  QUAC.MSAT.5loc.WILD.genind <- QUAC.MSAT.WILD.genind[, loc = samp_5loc]
  
  # Store the subsetted genind object in the list
  genind_list[[i]] <- QUAC.MSAT.5loc.WILD.genind
  wildSamp_5loc <- genind_list[[i]]@tab
  wildComplete <- colSums(wildSamp_5loc, na.rm = TRUE)/(nrow(wildSamp_5loc)*2) 
  wildSubset <- wildComplete[wildComplete > 0]
  total <- vector(length = nrow(wildSamp_5loc))
  common <- vector(length = nrow(wildSamp_5loc))
  lowfreq <- vector(length = nrow(wildSamp_5loc))
  rare <- vector(length = nrow(wildSamp_5loc))
  wildSamp_5loc <- wildSamp_5loc[,which(colSums(wildSamp_5loc, na.rm = TRUE)!= 0)]
  resamp_category5loc <- array(dim = c(164,4,25))
  category <- colnames(categorymat_5loc)
  dimnames(resamp_category5loc) <- list(paste0("sample ", 1:nrow(categorymat_5loc)), category, paste0("replicate ",1:25))
  j <- length((dimnames(resamp_category5loc)[[3]]))
  for (q in 1:j){
    for(i in 1:nrow(wildSamp_5loc)){
      # browser()
      samp <- sample(nrow(wildSamp_5loc), size = i, replace = FALSE)
      samp <- wildSamp_5loc[samp,]
      # Use sample to randomly subsample the matrix of wild individuals
      # Now, measure the proportion of allelic representation in that sample
      if (i==1) {
        total[i] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(samp > 0)))])/length(names(wildSubset))
        common[i] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(samp > 0))))/length(which(wildSubset > 0.1))
        lowfreq[i] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(samp > 0))))/length(which(wildSubset > 0.01 & wildSubset < 0.1))
        rare[i] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(samp > 0))))/length(which(wildSubset < 0.01))
      } else{
        total[i] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(colSums(samp, na.rm=TRUE) > 0)))])/length(names(wildSubset))
        common[i] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(colSums(samp, na.rm=TRUE) > 0))))/length(which(wildSubset > 0.1))
        lowfreq[i] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(colSums(samp, na.rm=TRUE) > 0))))/length(which(wildSubset > 0.01 & wildSubset < 0.1))
        rare[i] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(colSums(samp, na.rm=TRUE) > 0))))/length(which(wildSubset < 0.01))
      }
    }
    categorymat_5loc <- cbind(total,common,lowfreq,rare)
    resamp_category5loc[,,q] <- categorymat_5loc
  }
  resamp_category5loc
}
