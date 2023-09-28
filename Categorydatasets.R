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

QUAC_data_type_name <- c("Q. acerifoloa genetic diversity: Microsatellites", "SNP Subset (R0)", "SNP Subset (R80)")
allele_cat_calcs <- array(dim = c(90,4,3))
dimnames(allele_cat_calcs)<-list(paste0("sample ",1:90),c("meanTotal","upper95","lower95","ciWidth"), QUAC_data_type_name)

imagesDirectory <- "C:/Users/gsalas/Documents/resampling_CIs/Code/Images/"
pdf(file = paste0(imagesDirectory, "MSATtotalAllelecatPlot.pdf"))
for (i in 1:length(data)) {
  meanTotalcat<-vector()
  upper95<-vector()
  lower95<-vector()
  ciWidth<-vector()
  meanVerycat <- vector()
  meanComcat <- vector()
  meanLowcat <- vector()
  meanRarecat <- vector()
  # nested for loop is soft coded to identify all populations sample sizes 
  # across replicates for each dataset and is assigned to the variable q
  # main calculations will be assigned an object name and iterate across each                   population sample size for each dataset
  # object names with categories are hard coded with their appropriate column number
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
  # cbind will bind the columns together, 
  # however each iteration in the previous for loop will only be captured if they are first placed into the array called allele_cat_calcs
  # objects are also created for the .csv file we eventually want to create for ciWidths at different population sample sizes
  placeholder<-cbind(meanTotalcat,upper95,lower95, ciWidth)
  allele_cat_calcs[,,i]<-placeholder
  meanCI <- vector()
  sample15 <- vector()
  sample30 <- vector()
  sample45 <- vector()
  sample60 <- vector()
  sample75 <- vector()
  sample90 <- vector()
  # the average ci width value and populations sample sizes for each dataset is calculated by subsetting from the array
  for (j in 1:3) {
    meanCI[j]<-mean(allele_cat_calcs[,"ciWidth",j])
    sample15[j]<-allele_cat_calcs[15,"ciWidth",j]
    sample30[j]<-allele_cat_calcs[30,"ciWidth",j]
    sample45[j]<-allele_cat_calcs[45,"ciWidth",j]
    sample60[j]<-allele_cat_calcs[60,"ciWidth",j]
    sample75[j]<-allele_cat_calcs[75,"ciWidth",j]
    sample90[j]<-allele_cat_calcs[90,"ciWidth",j]
  }
  # we can find the minimum point at which the total mean category for each dataset achieves 95% genetic diversity and give it the object name min_95totavg
  # set the file path of where the plots will be saved to and specify the file type
  # standard plotting and formatting 
  min_95totavg<-(min(which(meanTotalcat > 95)))
  png(file = paste0(imagesDirectory, QUAC_data_type_name[i], "total category mean plots.png"),width=930, height=438)
  plot(meanTotalcat, xlab = "Sample", ylab = "Genetic Diversity %", main=QUAC_data_type_name[i],
       col="black", pch = 16)
  lines(upper95, col = "red", lwd = 2, lty = "dashed")
  lines(lower95, col = "blue", lwd = 2, lty = "dashed")
  leg.txt = c("total freq average", "upper 95% confidence interval", "lower 95% confidence interval")
  legend("right",
         legend = leg.txt, cex = 0.75,
         fill = c("black","red","blue"))
  dev.off( )
  png(file = paste0(imagesDirectory, QUAC_data_type_name[i],"category mean plots.png"), width=930, height=438)
  plot(xlab = "Sample", ylab = "Frequency %", meanTotalcat, col="black",main=QUAC_data_type_name[i], pch=16)
  points(meanVerycat,col="blue", pch=16)
  points(meanComcat, col="green", pch = 16)
  points(meanLowcat, col="red", pch = 16)
  points(meanRarecat, col="pink", pch = 16)
  abline(h = 95, lty = "dotted", col = "orange")
  abline(v = min_95totavg, col = "black")
  leg.txt = c("Total allele", "Very common allele", "Common allele", "Low frequency allele", "Rare allele")
  legend("right", 
         legend = leg.txt, cex = 0.75,
         title = (sub = paste("95% point estimation is ", min_95totavg)),
         fill = c("black","blue", "green", "red", "pink", "orange")
  )
  dev.off( )
  # create a .csv and bind the columns together, specify file name 
  write.csv(cbind(QUAC_data_type_name,sample15,sample30,sample45,sample60,sample75,sample90,meanCI),file = "datasetCIwidths.csv")
}
str(allele_cat_calcs)
