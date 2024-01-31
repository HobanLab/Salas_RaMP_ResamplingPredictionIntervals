# This script contains resampling arrays of allele frequencies based on different categories of rarity for MSAT, complete SNP and subset SNP of Q. acerifolia provided by research assistant Austin Koontz.
# The script calculates means, quantiles, and generates plots for the total and other allele categories. 
#####################################################
# 2023/08/28 Exploring Allele Category MSAT Dataset #
#####################################################
# Set work directory by adding path file 
setwd("C:/Users/gsalas/Documents/resampling_CIs/Code/")
# QUAC.MSAT.Subset_resampArr.Rdata is a resampling array containing allele frequency category data that filters microsatellite loci shared between garden and wild populations of Q. acerifolia
QUAC_MSAT_Subset_resampArr <- readRDS("Datasets/QUAC_Subset_resampArrs/QUAC.MSAT.Subset_resampArr.Rdata")
# Display the strucutre of the object 'QUAC_MSAT_Subset_resampArr'. This is information about its type, dimensions, and content
# [sample number, allele category, replicate number]
str(QUAC_MSAT_Subset_resampArr)
# this represents the total allelic representation % for samples 1-90 replicate one
QUAC_MSAT_Subset_resampArr[1:90,1,1] 
# this represents the allelic representation % for each allele category of sample 10 replicate one
QUAC_MSAT_Subset_resampArr[10,1:5,1]
# this represents the allelic representation % for each allele category of sample 10 replicate two
QUAC_MSAT_Subset_resampArr[10,1:5,2]
# Apply the 'mean' function along the sample sizes of individuals of the total allele category, across all replicates
apply(QUAC_MSAT_Subset_resampArr[,1,],1,mean)
# declare objects to store values in vectors
totalallelecatmean <- vector()
totalallelecat95upper <- vector()
totalallelecat95lower <- vector()
# loop through each row of the resampling array
for(i in 1:nrow(QUAC_MSAT_Subset_resampArr)) {
  # calculate the mean of the total allele category along the i-th row across all replicates and store it in 'totalallelecatmean'
  totalallelecatmean[i] <- mean(QUAC_MSAT_Subset_resampArr[i,1,])
  # calculate the 95th percentile of the total allele category along the i-th row across all replicates and store it in 'totalallelecat95upper'
  totalallelecat95upper[i] <- quantile(QUAC_MSAT_Subset_resampArr[i,1,], 0.95)
  # calculate the 5th percentile of the total allele category along the i-th row across all replicates and store it in 'totalallelecat95lower'
  totalallelecat95lower[i] <- quantile(QUAC_MSAT_Subset_resampArr[i,1,], 0.05)
}
# print the results of each object
print(totalallelecatmean)
print(totalallelecat95upper)   
print(totalallelecat95lower)

# this displays the total allele category and its confidence intervals
# create the file path to upload images
imagesDirectory <- "C:/Users/gsalas/Documents/resampling_CIs/Code/Images/"
# create a PDF file for plotting results
pdf(file = paste0(imagesDirectory, "QUACMSATSubsettotalAllelecatPlot.pdf"))
# create a plot with the x-axis labeled "sample" and the y-axis labaled "Frequency%" and a title "MSAT Total category Confidence Interval"
plot(xlab = "Sample", ylab = "Frequency%", main = "Q. acerifolia MSAT (Subset) Total Allele Category Confidence Interval",
     totalallelecatmean, col = "black", pch = 16)
# add a dashed red line for the upper 95% confidence interval in red
lines(totalallelecat95upper, col = "red", lwd = 2, lty = "dashed")
# add a dashed blue line for the lower 95% confidence interval in blue
lines(totalallelecat95lower, col = "blue", lwd = 2, lty = "dashed")
# add a legend to the right of the plot
leg.txt = c("total freq average", "upper 95% confidence interval", "lower 95% confidence interval")
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
for (i in 1:nrow(QUAC_MSAT_Subset_resampArr)){
  # Calculate the mean of the total allele category along the i-th row across all replicates and store it in 'totalallelecatmean'
  totalallelecatmean[i] <- mean(QUAC_MSAT_Subset_resampArr[i,1,])
  # calculate the mean of the very common allele category along the i-th row across all replicates and store it in 'verycomallelecatmean'
  verycomallelecatmean[i] <- mean(QUAC_MSAT_Subset_resampArr[i,2,])
  # calculate the mean of the common allele category along the i-th row across all replicates and store it in 'comallelecatmean'
  comallelecatmean[i] <- mean(QUAC_MSAT_Subset_resampArr[i,3,])
  # calculate the mean of the low frequency allele category along the i-th row across all replicates and store it in 'lowfreqallelecatmean'
  lowfreqallelecatmean[i] <- mean(QUAC_MSAT_Subset_resampArr[i,4,])
  # calculate the mean of the rare allele category along the i-th row across all replicates and store it in 'rareallelecatmean'
  rareallelecatmean[i] <- mean(QUAC_MSAT_Subset_resampArr[i,5,])
}
# find the minimum sample which 'totalallelecatmean' is greater than 95
min_95totavg<-(min(which(totalallelecatmean > 95)))
# Create a PNG file for plotting results
png(file = paste0(imagesDirectory, "QUACMSATSubsetallelecatmeanplot.png"))
# Create a plot with the x-axis labeled "Sample", the y-axis labeled "Frequency%" and the title "MSAT Allele Category Averages"
plot(xlab = "Sample", ylab = "Frequency %", main = "QUAC MSAT (Subset) Allele Category Averages", totalallelecatmean, pch = 16, ylim = c(0,100))
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
# calculate the width of the confidence interval for each individual in the total allele category
CIwidth <- totalallelecat95upper - totalallelecat95lower
# calculate the average width of the confidence interval for each individaul in the total allele category
avgCI <- mean(CIwidth)
# combine the widths of the confidence interval at various sample sizes along with the average confidence interval width into a vector
c(CIwidth[15], CIwidth[30], CIwidth[45], CIwidth[60], CIwidth[75], CIwidth[90], avgCI)
# create an object to declare the vector
MSAT_Subset <- c(CIwidth[15], CIwidth[30], CIwidth[45], CIwidth[60], CIwidth[75], CIwidth[90], avgCI)


#############################################################
# 2023/08/30 Exploring Allele Category SNP Subset Dataset #
#############################################################
# Set work directory by adding path file 
setwd("C:/Users/gsalas/Documents/resampling_CIs/Code/")
# a singlenucleotide polymorphism resampling array containing allele frequency category data that used the reference processing approach, does not filter the number of samples that share 0% of loci, and does filter loci shared between garden and wild populations of Q. acerifolia
QUAC_SNP_R0_Subset_resampArr<-readRDS("Datasets/QUAC_Subset_resampArrs/QUAC.SNP.REF.R0.Subset_resampArr.Rdata")
# declare objects to store values in vectors
totalallelecatmean <- vector()
totalallelecat95upper <- vector()
totalallelecat95lower <- vector()
# loop through each row of the resampling array
for(i in 1:nrow(QUAC_SNP_R0_Subset_resampArr)) {
  # calculate the mean of the total allele category along the i-th row across all replicates and store it in 'totalallelecatmean'
  totalallelecatmean[i] <- mean(QUAC_SNP_R0_Subset_resampArr[i,1,])
  # calculate the 95th percentile of the total allele category along the i-th row across all replicates and store it in 'totalallelecat95upper'
  totalallelecat95upper[i] <- quantile(QUAC_SNP_R0_Subset_resampArr[i,1,], 0.95)
  # calculate the 5th percentile of the total allele category along the i-th row across all replicates and store it in 'totalallelecat95lower'
  totalallelecat95lower[i] <- quantile(QUAC_SNP_R0_Subset_resampArr[i,1,], 0.05)
}
# print the results of each object
print(totalallelecatmean)
print(totalallelecat95upper)
print(totalallelecat95lower)

# this displays the total allele category and its confidence intervals
# create the file path to upload images
imagesDirectory <- "C:/Users/gsalas/Documents/resampling_CIs/Code/Images/"
# create a PDF file for plotting results
pdf(file = paste0(imagesDirectory, "QUACSNPR0totalAllelecatPlot.pdf"))
# create a plot with the x-axis labeled "sample" and the y-axis labeled "Frequency%" and a title "SNP Subset Total Category Confidence Interval"
plot(xlab = "Sample", ylab = "Frequency%", main = "Q. Acerifolia SNP R0 Subset  Total Category Confidence Interval",
     totalallelecatmean, col = "black", pch = 16)
# add a dashed red line for the upper 95% confidence interval in red
lines(totalallelecat95upper, col = "red", lwd = 2, lty = "dashed")
# add a dashed blue line for the lower 95% confidence interval in blue
lines(totalallelecat95lower, col = "blue", lwd = 2, lty = "dashed")
# add a legend to the right of the plot
leg.txt = c("total freq average", "upper 95% confidence interval", "lower 95% confidence interval")
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
for (i in 1:nrow(QUAC_SNP_R0_Subset_resampArr)) {
  # Calculate the mean of the total allele category along the i-th row across all replicates and store it in 'totalallelecatmean'
  totalallelecatmean[i] <- mean(QUAC_SNP_R0_Subset_resampArr[i,1,])
  # calculate the mean of the very common allele category along the i-th row across all replicates and store it in 'verycomallelecatmean'
  verycomallelecatmean[i] <- mean(QUAC_SNP_R0_Subset_resampArr[i,2,])
  # calculate the mean of the common allele category along the i-th row across all replicates and store it in 'comallelecatmean'
  comallelecatmean[i] <- mean(QUAC_SNP_R0_Subset_resampArr[i,3,])
  # calculate the mean of the low frequency allele category along the i-th row across all replicates and store it in 'lowfreqallelecatmean'
  lowfreqallelecatmean[i] <- mean(QUAC_SNP_R0_Subset_resampArr[i,4,])
  # calculate the mean of the rare allele category along the i-th row across all replicates and store it in 'rareallelecatmean'
  rareallelecatmean[i] <- mean(QUAC_SNP_R0_Subset_resampArr[i,5,])
}
# find the minimum sample which 'totalallelecatmean' is greater than 95
min_95totavg<-(min(which(totalallelecatmean > 95)))
# Create a PNG file for plotting results
png(file = paste0(imagesDirectory, width = 8, height = 11, units="in", "QUACSNPR0Subsetallelecatmeanplot.png"))
# Create a plot with the x-axis labeled "Sample", the y-axis labeled "Frequency%" and the title "SNP Complete Allele Category Averages"
plot(xlab = "Sample", ylab = "Frequency %", main = "SNP Complete Allele Category Averages", totalallelecatmean, pch = 16, ylim = c(0,100))
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
# calculate the average width of the confidence interval for each individaul in the total allele category
avgCIwidth <- mean(CIwidth)
# combine the widths of the confidence interval at various sample sizes along with the average confidence interval width into a vector
SNP_Subset_R0 <- c(CIwidth[15], CIwidth[30], CIwidth[45], CIwidth[60], CIwidth[75], CIwidth[90], avgCIwidth)

###########################################################
# 2023/08/31 Exploring Allele Category SNP Subset Dataset #
###########################################################
# Set work directory by adding path file 
setwd("C:/Users/gsalas/Documents/resampling_CIs/Code/")
QUAC_SNP_R80_Subset_resampArr<-readRDS("Datasets/QUAC_Subset_resampArrs/QUAC.SNP.REF.R80.Subset_resampArr.Rdata")
# Display the strucutre of the object 'QUAC_MSAT_Complete_resampArr'. This is information about its type, dimensions, and content
# [sample number, allele category, replicate number]
str(QUAC_SNP_R80_Subset_resampArr)
# Apply the 'mean' function along the sample sizes of individuals of the total allele category, across all replicates
apply(QUAC_SNP_R80_Subset_resampArr[,1,],1,mean)
# declare objects to store values in vectors
totalallelecatmean <- vector()
totalallelecat95upper <- vector()
totalallelecat95lower <- vector()
# loop through each row of the resampling array
for(i in 1:nrow(QUAC_SNP_R80_Subset_resampArr)) {
  totalallelecatmean[i] <- mean(QUAC_SNP_R80_Subset_resampArr[i,1,])
  totalallelecat95upper[i] <- quantile(QUAC_SNP_R80_Subset_resampArr[i,1,], 0.95)
  totalallelecat95lower[i] <- quantile(QUAC_SNP_R80_Subset_resampArr[i,1,], 0.05)
}
print(totalallelecatmean)
print(totalallelecat95upper)   
print(totalallelecat95lower)

# this displays the total allele category and its confidence intervals
# create the file path to upload images
imagesDirectory <- "C:/Users/gsalas/Documents/resampling_CIs/Code/Images/"
# create a PDF file for plotting results
pdf(file = paste0(imagesDirectory, "QUACSNPR80SubsettotalAllelecatPlot.pdf"))
# create a plot with the x-axis labeled "sample" and the y-axis labaled "Frequency%" and a title "SNP Subset Total Category Confidence Interval"

plot(xlab = "Sample", ylab = "Frequency%", main = "SNP R80 Subset Total Category Confidence Interval",
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
for (i in 1:nrow(QUAC_SNP_R80_Subset_resampArr)) {
  # Calculate the mean of the total allele category along the i-th row across all replicates and store it in 'totalallelecatmean'
  totalallelecatmean[i] <- mean(QUAC_SNP_R80_Subset_resampArr[i,1,])
  # calculate the mean of the very common allele category along the i-th row across all replicates and store it in 'verycomallelecatmean'
  verycomallelecatmean[i] <- mean(QUAC_SNP_R80_Subset_resampArr[i,2,])
  # calculate the mean of the common allele category along the i-th row across all replicates and store it in 'comallelecatmean'
  comallelecatmean[i] <- mean(QUAC_SNP_R80_Subset_resampArr[i,3,])
  # calculate the mean of the low frequency allele category along the i-th row across all replicates and store it in 'lowfreqallelecatmean'
  lowfreqallelecatmean[i] <- mean(QUAC_SNP_R80_Subset_resampArr[i,4,])
  # calculate the mean of the rare allele category along the i-th row across all replicates and store it in 'rareallelecatmean'
  rareallelecatmean[i] <- mean(QUAC_SNP_R80_Subset_resampArr[i,5,])
}

# find the minimum sample which 'totalallelecatmean' is greater than 95
min_95totavg<-(min(which(totalallelecatmean > 95)))
# Create a PDF file for plotting results
pdf(file = paste0(imagesDirectory, "QUACSNPR80allelecatmeanplot.pdf"))
# Create a plot with the x-axis labeled "Sample", the y-axis labeled "Frequency%" and the title "SNP Subset Allele Category Averages"
plot(xlab = "Sample", ylab = "Frequency %", main = "SNP Subset Allele Category Averages", totalallelecatmean, pch = 16, ylim = c(0,100))
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
# Complete the PDF file
dev.off()
# calculate the width of the confidence interval for each individual in the total allele category
CIwidth <- totalallelecat95upper - totalallelecat95lower
# calculate the average width of the confidence interval for each individual in the total allele category
avgCIwidth <- mean(CIwidth)
# combine the widths of the confidence interval at various sample sizes along with the average confidence interval width into a vector
SNP_Subset_R80 <- c(CIwidth[15], CIwidth[30], CIwidth[45], CIwidth[60], CIwidth[75], CIwidth[90], avgCIwidth)
rbind
# Combine the confidence interval width objects by rows
QUACGMCIWidths <- rbind(MSAT_Subset, SNP_Subset_R0,SNP_Subset_R80)
# specify the name of each column by the sample size, and the mean value
colnames(QUACGMCIWidths) <- c("15", "30", "45", "60", "75", "90", "Mean")
# matrix showing the confidence interval widths for each genetic marker dataset of Q. acerifolia
QUACGMCIWidths
# save matrix as .csv file
write.csv(QUACGMCIWidths, file = "C:/Users/gsalas/Documents/resampling_CIs/Code/Outputs/QUAC_CI_values.csv")
