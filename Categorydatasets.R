#####################################################
# 2023/09/12 Exploring Allele Category MSAT Dataset #
#####################################################
# Set work directory by adding path file 
setwd("C:/Users/gsalas/Documents/resampling_CIs/Code/Datasets")
datasets <- list.files(path = "Subset", pattern = "Arr.Rdata", full.names =  TRUE)

data<-list(
  readRDS(datasets[1]),
  readRDS(datasets[2]),
  readRDS(datasets[3])
)
# this represents the allelic representation % for the total allele category, samples 1 through 10, replicate one
data[[3]][1:10,1,1]

QUAC_data_type_name <- c("MSAT", "SNP Subset (R0)", "SNP Subset (R80)")
allele_cat_calcs <- array(dim = c(90,4,3))
dimnames(allele_cat_calcs)<-list(paste0("sample ",1:90),c("meanTotal","upper95","lower95","ciWidth"), QUAC_data_type_name)

# imagesDirectory <- "C:/Users/gsalas/Documents/resampling_CIs/Code/Images/"
# pdf(file = paste0(imagesDirectory, "MSATtotalAllelecatPlot.pdf"))
for (i in 1:length(data)) {
  meanTotalcat<-vector()
  upper95<-vector()
  lower95<-vector()
  ciWidth<-vector()
  meanVerycat <- vector()
  meanComcat <- vector()
  meanLowcat <- vector()
  meanRarecat <- vector()
  for (q in 1:nrow(data[[i]][,1,])) {
    meanTotalcat[q]<-mean(data[[i]][q,1,])
    upper95[q]<-quantile(data[[i]][q,1,],0.95)
    lower95[q]<-quantile(data[[i]][q,1,],0.05)
    ciWidth[q]<-upper95[q]-lower95[q]
    meanVerycat[q]<-mean(data[[i]][q,2,])
    meanComcat[q]<-mean(data[[i]][q,3,])
    meanLowcat[q]<-mean(data[[i]][q,4,])
    meanRarecat[q]<-mean(data[[i]][q,5,])
  } 
  # browser()
  placeholder<-cbind(meanTotalcat,upper95,lower95, ciWidth)
  allele_cat_calcs[,,i]<-placeholder
  meanCI <- vector()
  sample15 <- vector()
  sample30 <- vector()
  sample45 <- vector()
  sample60 <- vector()
  sample75 <- vector()
  sample90 <- vector()
  for (i in 1:3) {
    meanCI[i]<-mean(allele_cat_calcs[,"ciWidth",i])
    sample15[i]<-allele_cat_calcs[15,"ciWidth",i]
    sample30[i]<-allele_cat_calcs[30,"ciWidth",i]
    sample45[i]<-allele_cat_calcs[45,"ciWidth",i]
    sample60[i]<-allele_cat_calcs[60,"ciWidth",i]
    sample75[i]<-allele_cat_calcs[75,"ciWidth",i]
    sample90[i]<-allele_cat_calcs[90,"ciWidth",i]
  }
  plot(meanTotalcat, xlab = "Sample", ylab = "Frequency%", main=QUAC_data_type_name[i],
       col="black", pch = 16)
  lines(upper95, col = "red", lwd = 2, lty = "dashed")
  lines(lower95, col = "blue", lwd = 2, lty = "dashed")
  leg.txt = c("total freq average", "upper 95% confidence interval", "lower 95% confidence interval")
  legend("bottomright", 
         legend = leg.txt, cex = 0.75,
         fill = c("black","red","blue"))
  plot(meanTotalcat,col="black",main=QUAC_data_type_name[i])
  points(meanVerycat,col="blue")
  write.csv(cbind(sample15,sample30,sample45,sample60,sample75,sample90,meanCI),file = "datasetCIwidths.csv")
}

meanCI
write.csv(cbind(sample15,meanCI),file = "datasetCIwidths.csv")
sample15
str(allele_cat_calcs)
