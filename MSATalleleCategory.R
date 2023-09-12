#####################################################
# 2023/08/28 Exploring Allele Category MSAT Dataset #
#####################################################
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

# this displays the total allele category and its confidence intervals
imagesDirectory <- "C:/Users/gsalas/Documents/resampling_CIs/Code/Images/"
pdf(file = paste0(imagesDirectory, "MSATtotalAllelecatPlot.pdf"))
plot(xlab = "Sample", ylab = "Frequency%", main = "MSAT Total Category Confidence Interval",
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


# for (q in colMeans(QUAC_MSAT_Complete_resampArr)) {
#   for (i in 1:nrow(QUAC_MSAT_Complete_resampArr)){
#     totalallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,q,])
#     verycomallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,q,])
#     comallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,q,])
#     lowfreqallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,q,])
#     rareallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,q,])
#   }
# }  


  for (i in 1:nrow(QUAC_MSAT_Complete_resampArr)){
    totalallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,1,])
    verycomallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,2,])
    comallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,3,])
    lowfreqallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,4,])
    rareallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,5,])
  }

min_95totavg<-(min(which(totalallelecatmean > 95)))
png(file = paste0(imagesDirectory, "MSATallelecatmeanplot.png"))
plot(xlab = "Sample", ylab = "Frequency %", main = "MSAT Allele Category Averages", totalallelecatmean, pch = 16, ylim = c(0,100))
points(verycomallelecatmean, col="blue", pch = 16)
points(comallelecatmean, col="green", pch = 16)
points(lowfreqallelecatmean, col="red", pch = 16)
points(rareallelecatmean, col="pink", pch = 16)
abline(h = 95, lty = "dotted", col = "orange")
abline(v = min_95totavg, col = "black")
leg.txt = c("Total allele", "Very common allele", "Common allele", "Low frequency allele", "Rare allele")
legend("right", 
       legend = leg.txt, cex = 0.75,
       title = (sub = paste("95% point estimation is ", min_95totavg)),
       fill = c("black","blue", "green", "red", "pink", "orange"))
dev.off()

CIwidth <- totalallelecat95upper - totalallelecat95lower
avgCI <- mean(CIwidth)
cbind(CIwidth[15], CIwidth[30], CIwidth[45], CIwidth[60], CIwidth[75], CIwidth[90], CIwidth[105],CIwidth[120], CIwidth[135], CIwidth[150], avgCI)

ciWidthTable <- cbind(CIwidth[15], CIwidth[30], CIwidth[45], CIwidth[60], CIwidth[75], CIwidth[90], CIwidth[105],CIwidth[120], CIwidth[135], CIwidth[150], avgCI)
write.csv2(ciWidthTable, file="ciWidth.csv")
# allelecatmean <- vector()
# for (q in 1:5) {
#   for (i in 1:nrow(QUAC_MSAT_Complete_resampArr)){
#     allelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,q,])
#   }
# }
# allelecatmean