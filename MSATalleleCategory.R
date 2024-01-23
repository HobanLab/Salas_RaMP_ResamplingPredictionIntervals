#####################################################
# 2023/08/28 Exploring Allele Category MSAT Dataset #
#####################################################
# Set work directory by adding path file 
setwd("C:/Users/gsalas/Documents/resampling_CIs/Code/")
QUAC_MSAT_Complete_resampArr <- readRDS("Datasets/QUAC.MSAT.Complete_resampArr.Rdata")
# Display the strucutre of the object 'QUAC_MSAT_Complete_resampArr'. This is information about its type, dimensions, and content
# [sample number, allele category, replicate number]
str(QUAC_MSAT_Complete_resampArr)
# this represents the total allelic representation % for samples 1-100 replicate one
QUAC_MSAT_Complete_resampArr[1:100,1,1] 
# this represents the allelic representation % for each allele category of sample 10 replicate one
QUAC_MSAT_Complete_resampArr[10,1:5,1]
# this represents the allelic representation % for each allele category of sample 10 replicate two
QUAC_MSAT_Complete_resampArr[10,1:5,2]
# Apply the 'mean' function along the sample sizes of individuals of the total allele category, across all replicates
apply(QUAC_MSAT_Complete_resampArr[,1,],1,mean)
# declare objects to store values in vectors
totalallelecatmean <- vector()
totalallelecat95upper <- vector()
totalallelecat95lower <- vector()
# loop through each row of the resampling array
for(i in 1:nrow(QUAC_MSAT_Complete_resampArr)) {
  # calculate the mean of the total allele category along the i-th row across all replicates and store it in 'totalallelecatmean'
  totalallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,1,])
  # calculate the 95th percentile of the total allele category along the i-th row across all replicates and store it in 'totalallelecat95upper'
  totalallelecat95upper[i] <- quantile(QUAC_MSAT_Complete_resampArr[i,1,], 0.95)
  # calculate the 5th percentile of the total allele category along the i-th row across all replicates and store it in 'totalallelecat95lower'
  totalallelecat95lower[i] <- quantile(QUAC_MSAT_Complete_resampArr[i,1,], 0.05)
}
# print the results of each object
print(totalallelecatmean)
print(totalallelecat95upper)   
print(totalallelecat95lower)

# this displays the total allele category and its confidence intervals
# create the file path to upload images
imagesDirectory <- "C:/Users/gsalas/Documents/resampling_CIs/Code/Images/"
# create a PDF file for plotting results
pdf(file = paste0(imagesDirectory, "MSATtotalAllelecatPlot.pdf"))
# create a plot with the x-axis labeled "sample" and the y-axis labaled "Frequency%" and a title "MSAT Total category Confidence Interval
plot(xlab = "Sample", ylab = "Frequency%", main = "MSAT Total Category Confidence Interval",
     totalallelecatmean, col = "black", pch = 16)
# add a dashed red line for the upper 95% confidence interval in red
lines(totalallelecat95upper, col = "red", lwd = 2, lty = "dashed")
# add a dashed blue line for the lower 95% confidence interval in blue
lines(totalallelecat95lower, col = "blue", lwd = 2, lty = "dashed")
leg.txt = c("total freq average", "upper 95% confidence interval", "lower 95% confidence interval")
# add a legend to the right of the plot
legend("right", 
       legend = leg.txt, cex = 0.75,
       fill = c("black","red","blue"))
# complete the PDF file
dev.off()

# declare objects to store values in vectors
totalallelecatmean <- vector()
verycomallelecatmean <- vector()
comallelecatmean <- vector()
lowfreqallelecatmean <- vector()
rareallelecatmean <- vector()
# Loop through each row of the resampling array
  for (i in 1:nrow(QUAC_MSAT_Complete_resampArr)){
    # Calculate the mean of the total allele category along the i-th row across all replicates and store it in 'totalallelecatmean'
    totalallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,1,])
    # calculate the mean of the very common allele category along the i-th row across all replicates and store it in 'verycomallelecatmean'
    verycomallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,2,])
    # calculate the mean of the common allele category along the i-th row across all replicates and store it in 'comallelecatmean'
    comallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,3,])
    # calculate the mean of the low frequency allele category along the i-th row across all replicates and store it in 'lowfreqallelecatmean'
    lowfreqallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,4,])
    # calculate the mean of the rare allele category along the i-th row across all replicates and store it in 'rareallelecatmean'
    rareallelecatmean[i] <- mean(QUAC_MSAT_Complete_resampArr[i,5,])
  }
# find the minimum sample which 'totalallelecatmean' is greater than 95
min_95totavg<-(min(which(totalallelecatmean > 95)))
# Create a PNG file for plotting results
png(file = paste0(imagesDirectory, "MSATallelecatmeanplot.png"))
# Create a plot with the x-axis labeled "Sample", the y-axis labeled "Frequency%" and the title "MSAT Allele Category Averages"
plot(xlab = "Sample", ylab = "Frequency %", main = "MSAT Allele Category Averages", totalallelecatmean, pch = 16, ylim = c(0,100))
# Add points for different allele categories in different colors
points(verycomallelecatmean, col="blue", pch = 16)
points(comallelecatmean, col="green", pch = 16)
points(lowfreqallelecatmean, col="red", pch = 16)
points(rareallelecatmean, col="pink", pch = 16)
# add a horizontal line at 95%
abline(h = 95, lty = "dotted", col = "orange")
# add a vertical line at the minimum sample which 'totalallelecatmean is greater than 95
abline(v = min_95totavg, col = "black")
# add a legend to the plot 
leg.txt = c("Total allele", "Very common allele", "Common allele", "Low frequency allele", "Rare allele")
legend("right", 
       legend = leg.txt, cex = 0.75,
       title = (sub = paste("95% point estimation is ", min_95totavg)),
       fill = c("black","blue", "green", "red", "pink", "orange"))
# Complete the PNG file
dev.off()
# caclulate the width of the confidence interval for each individual in the total allele category
CIwidth <- totalallelecat95upper - totalallelecat95lower
# calculate the average width of the confidence interval for each individual in the total allele category
avgCI <- mean(CIwidth)
# combine the widths of the confidence interval at various sample sizes along with the average confidence interval width into a matrix
cbind(CIwidth[15], CIwidth[30], CIwidth[45], CIwidth[60], CIwidth[75], CIwidth[90], CIwidth[105],CIwidth[120], CIwidth[135], CIwidth[150], avgCI)
# create an object to declare the matrix
ciWidthTable <- cbind(CIwidth[15], CIwidth[30], CIwidth[45], CIwidth[60], CIwidth[75], CIwidth[90], CIwidth[105],CIwidth[120], CIwidth[135], CIwidth[150], avgCI)
# write the matrix to a .csv file 
write.csv2(ciWidthTable, file="ciWidth.csv")