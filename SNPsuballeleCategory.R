###########################################################
# 2023/08/31 Exploring Allele Category SNP Subset Dataset #
###########################################################
# Set work directory by adding path file 
setwd("C:/Users/gsalas/Documents/resampling_CIs/Code/")
QUAC_SNP_DN_Subset_resampArr<-readRDS("Datasets/QUAC.SNP.DN.R80.Subset_resampArr.Rdata")
str(QUAC_SNP_DN_Subset_resampArr)

apply(QUAC_SNP_DN_Subset_resampArr[,1,],1,mean)

totalallelecatmean <- vector()
totalallelecat95upper <- vector()
totalallelecat95lower <- vector()

for(i in 1:nrow(QUAC_SNP_DN_Subset_resampArr)) {
  totalallelecatmean[i] <- mean(QUAC_SNP_DN_Subset_resampArr[i,1,])
  totalallelecat95upper[i] <- quantile(QUAC_SNP_DN_Subset_resampArr[i,1,], 0.95)
  totalallelecat95lower[i] <- quantile(QUAC_SNP_DN_Subset_resampArr[i,1,], 0.05)
}

totalallelecatmean

imagesDirectory <- "C:/Users/gsalas/Documents/resampling_CIs/Code/Images/"
pdf(file = paste0(imagesDirectory, "SubsetSNPtotalAllelecatPlot.pdf"))
plot(xlab = "Sample", ylab = "Frequency%", main = "SNP Subset Total Category Confidence Interval",
     totalallelecatmean, col = "black", pch = 16)
lines(totalallelecat95upper, col = "red", lwd = 2, lty = "dashed")
lines(totalallelecat95lower, col = "blue", lwd = 2, lty = "dashed")
leg.txt = c("total freq average", "upper 95% confidence interval", "lower 95% confidence interval")
legend("right", 
       legend = leg.txt, cex = 0.75,
       fill = c("black","red","blue"))
dev.off()

totalallelecatmean <- vector()
verycomallelecatmean <- vector()
comallelecatmean <- vector()
lowfreqallelecatmean <- vector()
rareallelecatmean <- vector()

for (i in 1:nrow(QUAC_SNP_DN_Subset_resampArr)) {
  totalallelecatmean[i] <- mean(QUAC_SNP_DN_Subset_resampArr[i,1,])
  verycomallelecatmean[i] <- mean(QUAC_SNP_DN_Subset_resampArr[i,2,])
  comallelecatmean[i] <- mean(QUAC_SNP_DN_Subset_resampArr[i,3,])
  lowfreqallelecatmean[i] <- mean(QUAC_SNP_DN_Subset_resampArr[i,4,])
  rareallelecatmean[i] <- mean(QUAC_SNP_DN_Subset_resampArr[i,5,])
}

min_95totavg<-(min(which(totalallelecatmean > 95)))
pdf(file = paste0(imagesDirectory, "SubsetSNPallelecatmeanplot.pdf"))
plot(xlab = "Sample", ylab = "Frequency %", main = "SNP Subset Allele Category Averages", totalallelecatmean, pch = 16, ylim = c(0,100))
points(verycomallelecatmean, col="blue", pch = 16)
points(comallelecatmean, col="green", pch = 16)
points(lowfreqallelecatmean, col="red", pch = 16)
points(rareallelecatmean, col="pink", pch = 16)
abline(h = 95, lty = "dotted", col = "orange")
abline(v = min_95totavg, col = "black")
leg.txt = c("Total allele", "Very common allele", "Common allele", "Low frequency allele", "Rare allele")
legend("right", 
       legend = leg.txt, cex = 0.75,
       title = (sub = paste("95% point estimation is ", x)),
       fill = c("black","blue", "green", "red", "pink", "orange"))
dev.off()