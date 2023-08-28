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

totalallelecatmean <- vector()
for (i in 1:nrow(QUAC_MSAT_Complete_resampArr)) {
  totalallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,1,])
  
}
totalallelecatmean

apply(QUAC_MSAT_Complete_resampArr[,1,],1,mean)

totalallelecat95upper <- vector()
totalallelecat95lower <- vector()
for(i in 1:nrow(QUAC_MSAT_Complete_resampArr)) {
  totalallelecat95upper[i] <- quantile(QUAC_MSAT_Complete_resampArr[i,1,], 0.95)
  totalallelecat95lower[i] <- quantile(QUAC_MSAT_Complete_resampArr[i,1,], 0.05)
}

totalallelecat95upper   
totalallelecat95lower

imagesDirectory <- "C:/Users/gsalas/Documents/resampling_CIs/Code/Images/"
pdf(file = paste0(imagesDirectory, "totalAllelecatPlot.pdf"))
plot(xlab = "", ylab = "", totalallelecatmean, xlim = c(0,150))
lines(totalallelecat95upper, col = "red", lwd = 2, lty = "dashed")
lines(totalallelecat95lower, col = "blue", lwd = 2, lty = "dashed")
dev.off()