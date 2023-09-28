########################
# 2023/09/28 
########################
setwd("C:/Users/gsalas/Documents/resampling_CIs/Code/Datasets")
datasets <- list.files(path = "Subset", pattern = "Arr.Rdata", full.names =  TRUE)
gm_Array <- list(
  readRDS(datasets[1]),
  readRDS(datasets[2]),
  readRDS(datasets[3])
)
for(i in 1:length(gm_Array)) {
  gm_Totals <- c(gm_Array[[i]][,1,])
  gm_sampleNumbers <- 2:(nrow(gm_Array[[i]][,1,])+1)
  gm_sampleNumbers <- rep(gm_sampleNumbers, dim(gm_Array[[i]])[[3]])
  gm_DF <- data.frame(sampleNumbers=gm_sampleNumbers, totalValues=gm_Totals)
  # for (j in 1:length(data)) {
  #   gm_model <- lm(sampleNumbers ~ I((totalValues)^3), data=gm_DF)
  # }
  gm_Model <- lm(sampleNumbers ~ I((totalValues)^3), data = gm_DF)
  gm_newData <- data.frame(totalValues=95)
  gm_95MSSEprediction <- predict(gm_Model, gm_newData, interval = "prediction")
}
gm_95MSSEprediction
str(gm_DF)
##########
# for(i in 1:length(gm_Array)) {
#   gm_Totals<-vector()
#   gm_sampleNumbers<-vector()
#   gm_DF<-vector()
#   gm_Model<-vector()
#   gm_newData<-vector()
#   gm_95MSSEprediction<-vec
#   for (q in 1:nrow(data[[i]][,1,])) {
#     gm_Totals <- c(gm_Array[[i]][,1,])
#     gm_sampleNumbers <- 2:(nrow(gm_Array[[i]][,1,])+1)
#     gm_sampleNumbers <- rep(gm_sampleNumbers, dim(gm_Array[[i]])[[3]])
#     gm_DF <- data.frame(sampleNumbers=gm_sampleNumbers, totalValues=gm_Totals)
#     gm_Model <- lm(sampleNumbers ~ I((totalValues)^3), data = gm_DF)
#     gm_newData <- data.frame(totalValues=95)
#     gm_95MSSEprediction <- predict(gm_Model, gm_newData, interval = "prediction")
#   }
# }