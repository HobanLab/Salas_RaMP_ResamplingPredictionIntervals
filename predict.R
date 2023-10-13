#####################################################################
# 2023/09/28 Predict function over multiple genetic marker datasets #
#####################################################################
setwd("C:/Users/gsalas/Documents/resampling_CIs/Code/Datasets")
datasets <- list.files(path = "Subset", pattern = "Arr.Rdata", full.names =  TRUE)
#gm stands for genetic marker. 
# read in datasets
gm_listArray <- list(
  readRDS(datasets[1]),
  readRDS(datasets[2]),
  readRDS(datasets[3])
)
# Declare an object to capture total resampling values
gm_Totals <- vector(length = nrow(gm_listArray[[1]]))
# Declare an object to capture the fit, upper, and lower predict values 
genetic_marker_predict_matrix = matrix(ncol = 3, nrow = 3)
for(i in 1:length(gm_listArray)) {
  # Capture Total resampling values
  gm_Totals <- vector(length = (nrow(gm_listArray[[i]]) * (dim(gm_listArray[[i]])[[3]])))
  gm_Totals <- c(gm_listArray[[i]][,1,])
  # Specify sample numbers column
  gm_sampleNumbers <- 2:(nrow(gm_listArray[[i]][,1,])+1)
  gm_sampleNumbers <- rep(gm_sampleNumbers, dim(gm_listArray[[i]])[[3]])
  # Create data.frame from resampling array values
  gm_DF <- data.frame(sampleNumbers=gm_sampleNumbers, totalValues=gm_Totals)
  # Build and analyze linear models
  gm_Model <- lm(sampleNumbers ~ I((totalValues)^3), data = gm_DF)
  gm_newData <- data.frame(totalValues=95)
  gm_95MSSEprediction <- predict(gm_Model, gm_newData, interval = "prediction")
  # Pass the gm_95MSSEprediction to the object storing our results by iterating
  # storing them in the rows index of predict_matrix
  genetic_marker_predict_matrix[i,] <- gm_95MSSEprediction
}
genetic_marker_predict_matrix
colnames(genetic_marker_predict_matrix) <- c("fit", "lwr", "upr")
write.csv(genetic_marker_predict_matrix, file = "C:/Users/gsalas/Documents/resampling_CIs/Code/gm_predict.csv", 
          row.names = FALSE)

# Set work directory by adding path file 
setwd("C:/Users/gsalas/Documents/resampling_CIs/Code/")
# Load dataset into environment
load("Datasets/quercus_final_results_orig.Rdata")
# 
quercus14_genDiv <- vector(length = (nrow(final_quercus_results)))
quercus14_predict_Matrix <- matrix(ncol = 3, nrow = 14)
for (q in 1:(dim(final_quercus_results)[3])) {
  # Capture Total resampling values
  quercus14_genDiv <- vector(length = (nrow(final_quercus_results[,,q])) * (dim(final_quercus_results)[2]))
  quercus14_genDiv <- c(final_quercus_results[,,q])
  # Specify sample numbers column
  quercus14_sampleNumbers <- 2:(nrow(final_quercus_results[,,q])+1)
  quercus14_sampleNumbers <- rep(quercus14_sampleNumbers, dim(final_quercus_results[,,q])[2])
  # Create data.frame from resampling array values
  quercus14_DF <- data.frame(sampleNumbers=quercus14_sampleNumbers, totalValues=quercus14_genDiv)
  # quercus14_DF_list[[q]] <- quercus14_DF
  # Build and analyze linear models
  quercus14_Model <- lm(sampleNumbers ~ I((totalValues)^3), data = quercus14_DF)
  # ensure then value being predicted is 0.95
  quercus14_newData <- data.frame(totalValues=0.95)
  quercus14_95MSSEprediction <- predict(quercus14_Model, quercus14_newData, interval = "prediction")
  # Pass the gm_95MSSEprediction to the object storing our results by iterating
  # storing them in the rows index of predict_matrix
  quercus14_predict_Matrix[q,] <- quercus14_95MSSEprediction
}
quercus14_predict_Matrix
colnames(quercus14_predict_Matrix) <- c("fit", "lwr", "upr")
write.csv(quercus14_predict_Matrix, file = "C:/Users/gsalas/Documents/resampling_CIs/Code/q14_predict.csv", 
          row.names = FALSE)