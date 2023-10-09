#####################################################################
# 2023/09/28 Predict function over multiple genetic marker datasets #
#####################################################################
setwd("C:/Users/gsalas/Documents/resampling_CIs/Code/Datasets")
datasets <- list.files(path = "Subset", pattern = "Arr.Rdata", full.names =  TRUE)
#gm stands for genetic marker
gm_Array <- list(
  readRDS(datasets[1]),
  readRDS(datasets[2]),
  readRDS(datasets[3])
)

# Declare an object to capture total resampling values
gm_Totals <- vector(length = nrow(gm_Array[[1]]))
# Declare an object to capture the fit, upper, and lower predict values 
predict_matrix = matrix(ncol = 3, nrow = 3)
for(i in 1:length(gm_Array)) {
  # Capture Total resampling values
  gm_Totals <- vector(length = (nrow(gm_Array[[i]]) * (dim(gm_Array[[i]])[[3]])))
  gm_Totals <- c(gm_Array[[i]][,1,])
  # Specify sample numbers column
  gm_sampleNumbers <- 2:(nrow(gm_Array[[i]][,1,])+1)
  gm_sampleNumbers <- rep(gm_sampleNumbers, dim(gm_Array[[i]])[[3]])
  # Create data.frame from resampling array values
  gm_DF <- data.frame(sampleNumbers=gm_sampleNumbers, totalValues=gm_Totals)
  # Build and analyze linear models
  gm_Model <- lm(sampleNumbers ~ I((totalValues)^3), data = gm_DF)
  gm_newData <- data.frame(totalValues=95)
  gm_95MSSEprediction <- predict(gm_Model, gm_newData, interval = "prediction")
  # Pass the gm_95MSSEprediction to the object storing our results by iterating
  # storing them in the rows index of predict_matrix
  predict_matrix[i,] <- gm_95MSSEprediction
  # for (j in 1:3) {
  #   for (z in 1:3) {
  #     predict_matrix[j,z] = gm_95MSSEprediction[j]
  #   }
  # }
}
predict_matrix

# Set work directory by adding path file 
setwd("C:/Users/gsalas/Documents/resampling_CIs/Code/")
# Load dataset into environment
load("Datasets/quercus_final_results_orig.Rdata")
