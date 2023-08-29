################################################
# 2023/08/28 Exploring Allele Category Dataset #
################################################
# Set work directory by adding path file 
setwd("C:/Users/gsalas/Documents/resampling_CIs/Code/")
QUAC_MSAT_Complete_resampArr <- readRDS("Datasets/QUAC.MSAT.Complete_resampArr.Rdata")
# [sample number, allele category, replicate number]
# this represents the total allelic representation % for samples 1-100 replicate one
QUAC_MSAT_Complete_resampArr[1:100,1,1] 
# this represents the allelic representation % for each allele category of sample 10 replicate one
QUAC_MSAT_Complete_resampArr[10,1:5,1]
# this represents the allelic representation % for each allele category of sample 10 replicate two
QUAC_MSAT_Complete_resampArr[10,1:5,2]

str(QUAC_MSAT_Complete_resampArr)

apply(QUAC_MSAT_Complete_resampArr[,1,],1,mean)

totalallelecatmean <- vector()
totalallelecat95upper <- vector()
totalallelecat95lower <- vector()
for(i in 1:nrow(QUAC_MSAT_Complete_resampArr)) {
  totalallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,1,])
  totalallelecat95upper[i] <- quantile(QUAC_MSAT_Complete_resampArr[i,1,], 0.95)
  totalallelecat95lower[i] <- quantile(QUAC_MSAT_Complete_resampArr[i,1,], 0.05)
}
totalallelecatmean
totalallelecat95upper   
totalallelecat95lower

imagesDirectory <- "C:/Users/gsalas/Documents/resampling_CIs/Code/Images/"
pdf(file = paste0(imagesDirectory, "totalAllelecatPlot.pdf"))
plot(xlab = "", ylab = "", totalallelecatmean)
lines(totalallelecat95upper, col = "red", lwd = 2, lty = "dashed")
lines(totalallelecat95lower, col = "blue", lwd = 2, lty = "dashed")
dev.off()

totalallelecatmean <- vector()
verycomallelecatmean <- vector()
comallelecatmean <- vector()
lowfreqallelecatmean <- vector()
rareallelecatmean <- vector()


for (q in colMeans(QUAC_MSAT_Complete_resampArr)) {
  for (i in 1:nrow(QUAC_MSAT_Complete_resampArr)){
    totalallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,q,])
    verycomallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,q,])
    comallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,q,])
    lowfreqallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,q,])
    rareallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,q,])
  }
}  


  for (i in 1:nrow(QUAC_MSAT_Complete_resampArr)){
    totalallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,1,])
    verycomallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,2,])
    comallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,3,])
    lowfreqallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,4,])
    rareallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,5,])
  }

pdf(file = paste0(imagesDirectory, "allelecatmeanplot.pdf"))
plot(xlab = "Sample", ylab = "Frequency %", main = "Allele Category Averages", totalallelecatmean, pch = 16)
points(verycomallelecatmean, col="blue", pch = 16)
points(comallelecatmean, col="green", pch = 16)
points(lowfreqallelecatmean, col="red", pch = 16)
points(rareallelecatmean, col="brown", pch = 16)
abline(h = 95, lty = "dotted", col = "orange")
abline(v = min(which(totalallelecatmean > 95)), col = "black")
leg.txt = c("Total allele", "Very common allele", "Common allele", "Low frequency allele", "Rare allele")
legend(100, 60, legend = leg.txt,
       fill = c("black","blue", "green", "red", "brown", "orange"))
dev.off()

# allelecatmean <- vector()
# for (q in 1:5) {
#   for (i in 1:nrow(QUAC_MSAT_Complete_resampArr)){
#     allelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,q,])
#   }
# }
# allelecatmean
