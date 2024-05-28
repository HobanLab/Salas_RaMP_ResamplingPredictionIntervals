#####################################################################
# 2023/09/28 Predict function over multiple genetic marker datasets #
#####################################################################
setwd("C:/Users/gsalas/Documents/resampling_CIs/Code/Datasets/QUAC_Subset_resampArrs")
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
  ciWidth <- vector()
  for (j in 1:(nrow(genetic_marker_predict_matrix))) {
    ciWidth[j] <- genetic_marker_predict_matrix[j,3] - genetic_marker_predict_matrix[j,2]
  }
}
colnames(genetic_marker_predict_matrix) <- c("fit", "lwr", "upr")
rownames(genetic_marker_predict_matrix) <- c("MSAT", "SNP Subset (R0)", "SNP Subset (R80)")
gmMatrix <- cbind(genetic_marker_predict_matrix, ciWidth)
write.csv(gmMatrix, file = "C:/Users/gsalas/Documents/resampling_CIs/Code/Outputs/QUAC_PI_Values.csv", 
          row.names = TRUE)

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
  # Build and analyze linear models
  quercus14_Model <- lm(sampleNumbers ~ I((totalValues)^3), data = quercus14_DF)
  # ensure then value being predicted is 0.95
  quercus14_newData <- data.frame(totalValues=0.95)
  quercus14_95MSSEprediction <- predict(quercus14_Model, quercus14_newData, interval = "prediction")
  # Pass the gm_95MSSEprediction to the object storing our results by iterating
  # storing them in the rows index of predict_matrix
  quercus14_predict_Matrix[q,] <- quercus14_95MSSEprediction
  piWidth <- vector()
  for (z in 1:(nrow(quercus14_predict_Matrix))) {
    piWidth[z] <- quercus14_predict_Matrix[z,3] - quercus14_predict_Matrix[z,2]
  } 
}
colnames(quercus14_predict_Matrix) <- c("fit", "lwr", "upr")
rownames(quercus14_predict_Matrix) <- c("QUAC","QUAR","QUAU", "QUBO","QUCA","QUCE","QUEN","QUGE","QUGR","QUHA","QUHI","QUOG","QUPA", "QUTO")
q14Matrix <- cbind(quercus14_predict_Matrix, piWidth)
write.csv(q14Matrix, file = "C:/Users/gsalas/Documents/resampling_CIs/Code/Outputs/Quercus14_PI_values.csv", 
          row.names = TRUE)
