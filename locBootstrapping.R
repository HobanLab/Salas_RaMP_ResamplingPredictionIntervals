#####################################
# GDS 2023/12/06 Loci Bootstrapping #
#####################################
setwd("C:/Users/gsalas/Documents/resampling_CIs/Code/Datasets")
library(adegenet)

test <- readRDS("resamp_category.RDS")

# totalallelecatmean <- vector()
# comallelecatmean <- vector()
# lowfreqallelecatmean <- vector()
# rareallelecatmean <- vector()
# totalallelecat95upper <- vector()
# totalallelecat95lower <- vector()
# 
# j <- length(dimnames(test)[[2]])
# 
# for(i in 1:nrow(test)) {
#   totalallelecatmean[i] <- mean(test[i,"total",])
#   totalallelecat95upper[i] <- quantile(test[i,"total",], 0.95)
#   totalallelecat95lower[i] <- quantile(test[i,"total",], 0.05)
# }
# 
# for (i in 1:nrow(test)) {
#     totalallelecatmean[i] <- mean(test[i,"total",])
#     comallelecatmean[i] <- mean(test[i,"common",])
#     lowfreqallelecatmean[i] <- mean(test[i,"lowfreq",])
#     rareallelecatmean[i] <- mean(test[i,"rare",])
# }
# 
# mat <- cbind(totalallelecatmean,comallelecatmean,lowfreqallelecatmean,rareallelecatmean)

# linear model of resampling array
totalsVector <- c(test[,"total",])
# Specify sample numbers column
gm_sampleNumbers <- 1:(nrow(test[,"total",]))
gm_sampleNumbers <- rep(gm_sampleNumbers, dim(test)[[3]])
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

# Randomly sample a subset of loci
# ANSWER:
QUAC.MSAT.Wild.genind <- QUAC.MSAT.genind[1,]
locNames(QUAC.MSAT.genind)
loc_samp <- sample(length(locNames(QUAC.MSAT.genind)), size = 3, replace = FALSE)
loc_samp <- sample(locNames(QUAC.MSAT.genind), size = 3, replace = FALSE)
QUAC.MSAT.3loc.genind <- QUAC.MSAT.genind[,loc=loc_samp]
wildSamp_3loc <- QUAC.MSAT.3loc.genind@tab[wildRows,]
locNames(QUAC.MSAT.3loc.genind)
