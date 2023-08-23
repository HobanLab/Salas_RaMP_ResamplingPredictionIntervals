##########################################
# 2023/08/01 Exploring Resampling Dataset#
##########################################

# Set work directory by adding path file 
setwd("C:/Users/gsalas/Documents/resampling_CIs/Code/")

# Load dataset into environment
load("Datasets/quercus_final_results_orig.Rdata")

# List objects
ls()

# Examine dimensions of the dataset
dim(final_quercus_results)

# Amount of genetic diversity from a sample size of random individuals ranging between one to five hundred 
# for one thousand replicates of species one
final_quercus_results[,,1]

# This represents the amount of genetic diversity for a sample size of 
# random individuals [ranging between one to one hundred, for one replicate, of species one].
final_quercus_results[1:100,1,1]

# This represents the amount of genetic diversity for a sample size of 
# one random individual replicated one hundred times for species one.
final_quercus_results[1,1:100,1]

# This represents the amount of genetic diversity for a sample size of 
# random individuals ranging between one to one hundred replicated one time 
# for species two
final_quercus_results[1:100,1,2]
 
# This plot displays the genetic diversity of samples the size of 
# one hundred fifty individuals for all the replicates of species two.
plot(final_quercus_results[150,,2])

# This plot displays genetic diversity of all sample sizes for 
# one replicate of species one in blue and all sample sizes for one replicate
# of species two in red.
plot(final_quercus_results[,1,1],col="blue") 
points(final_quercus_results[,1,2],col="red")
# This plot displays the genetic diversity of samples the size of 
# one individual for all the replicates of species two
plot(final_quercus_results[1,,2])

# this gives us the mean genetic diversity values for a range of sample sizes, 
# range of replicates, and species' along the rows
sp <- 1
# we made a for loop alternative for the apply() function
apply(final_quercus_results[,,1],1,mean)
plot(final_quercus_results[,1,sp])
lines(apply(final_quercus_results[,,sp],1,mean),col="red",lwd=2)

# (loop alt) assign a variable to the vector to eventually capture all the means across the replicates
meanRepValues <- vector()

# nrow is soft code that detects how many rows there are in an object
# we use it to get the amount of genetic diversity across replicates 
# we then set an index for e in order to get a mean.
for(i in 1:nrow(final_quercus_results)){
  meanRepValues[i] <- mean(final_quercus_results[i,,1])
  }
# mean genetic diversity across all replicates for species 1
meanRepValues

# graph of meanRepvalues
plot(xlab = "Sample Size", ylab = "Genetic Diversity", meanRepValues)
leg.txt = c("Total Mean Genetic Diversity")
legend(250, 0.7, legend = leg.txt,
       fill = c("black"))

# graph of genetic diversity, replicates 1 - 3 are plotted along with a line showing the tmean gd
plot(xlab = "Sample size", ylab = "Genetic Diversity", final_quercus_results[,1,1], pch = 16)
points(final_quercus_results[,2,1], col="blue", pch = 16)
points(final_quercus_results[,3,1], col="green", pch = 16)
# recall this is the average and the above points are the individual replicates plotted on one graph
lines(meanRepValues, col="red", lwd=2)
leg.txt <- c("Replicate 1", "Replicate 2", "Replicate 3", "Total Mean of Genetic Diversity")
legend(200, 0.7, legend = leg.txt,
       fill = c("black","blue","green","red"))

# get the 95% CI of plot, but first go thru these ideas
# IDEA 1
# this gives us the position in the vector of the min genetic diversity value greater than 0.95 in replicate 1 of species 1
min(which(final_quercus_results[,1,sp]>0.95))

# this gives us the minimum sample size across all replicates
sp<-1; min_samp95<-vector(length = 1000)
for (r in 1:1000) {
  min_samp95[r]<-min(which(final_quercus_results[,r,sp]>0.95))
}

min_samp95

#this gives mean across reps of the first individual to cross 95%
boxplot(min_samp95)
mean(min_samp95)

# distribution of 95% min sample sizes across replicates for species one
# not what we quite want however...
plot(xlab = "Replicate number", ylab = "95% minimum sample size", min_samp95)
leg.txt <- "Sample size"
legend(650, 250, legend = leg.txt,
       fill = c("black"))

# IDEA 2
p<-1
# this shows the genetic diversity value in 95 percentile of the values across the replicates for a sample size of one for species one
quantile(final_quercus_results[p,,1],0.95)

quantile(final_quercus_results[p,,1],.05)

p<-2
quantile(final_quercus_results[p,,1],.95);
quantile(final_quercus_results[p,,1],.05)


# MSSE = minimum sample size estimate
# Using mean(min_samp95) ; incorrect (we think)
# we're plotting the average of the 95% MSSE for each replicate

# Using (min(which(meanRepValues > 0.95))) ; correct (we think)
# we're plotting the 95% MSEE for the *average* representation values across replicates
# we do NOT need to specify range in vector
# also, notice we can the meanRepValues from earlier to the function for efficiency
# declare vectors
meanRepValues <- vector()
upper95 <- vector()
lower95 <- vector()
for (n in 1:nrow(final_quercus_results)) {
  meanRepValues[n] <- mean(final_quercus_results[n,,1])
  upper95[n] <- quantile(final_quercus_results[n,,1],0.95)
  lower95[n] <- quantile(final_quercus_results[n,,1],0.05)
 }

# add the lines to the legend, add an asymptote of 0.95 horizontal, and add the line at which sample size reaches the 0.95 benchmark vertically
plot(xlab = "Sample Size", ylab = "Genetic Diversity", main = "QUAC", meanRepValues, ylim = c(0,1))
leg.txt = c("Total Mean Genetic Diversity", "95% Upper Limit", "95% Lower Limit")
lines(upper95, col = "red",lwd = 2, lty = "dashed")
lines(lower95, col = "green",lwd = 2, lty = "dashed")
abline(h = 0.95, lty = "dotted", col = "orange")
abline(v = min(which(meanRepValues > 0.95)), lty = "dotted", col = "orange")
legend(250, 0.7, legend = leg.txt,
       fill = c("black", "red", "green"))


# now, we want to declare a higher dimension object for the 14 slices (spp.) of the array for the 3 vectors
species_name <- c("QUAC","QUAR","QUAU", "QUBO","QUCA","QUCE","QUEN","QUGE","QUGR","QUHA","QUHI","QUOG","QUPA", "QUTO")
resultsArray <- array(dim = c(500,4,14))
dimnames(resultsArray)<-list(paste0("sample",1:500),c("meangd","upper95","lower95","ciwidth"), species_name)
meanRepValues <- vector()
upper95 <- vector()
lower95 <- vector()
for (q in 1:14) {
  for (i in 1:nrow(final_quercus_results)) {
    # 
    meanRepValues[i] <- mean(final_quercus_results[i,,q])
    upper95[i] <- quantile(final_quercus_results[i,,q],0.95)
    lower95[i] <- quantile(final_quercus_results[i,,q],0.05)
    CIwidth <- upper95 - lower95
  }
  # Bind vectors together into a matrix
  speciesMat <- cbind(meanRepValues, upper95, lower95, CIwidth)
  # Pass the matrix into a slot of the array
  resultsArray[,,q] <- speciesMat
}
str(resultsArray)

# 
imagesDirectory <- "C:/Users/gsalas/Documents/resampling_CIs/Code/Images/"
pdf(file=paste0(imagesDirectory,"14CIPlots.pdf"), width = 8.5, height = 11)
par(mfrow=c(3,2), omi=c(0,0.3,0,1.7))
for (i in 1:14) {
  x <- resultsArray[,,i]
  # plot(xlab = "Sample Size", ylab = "Genetic Diversity",x[,1], ylim = c(0,1))
  plot(xlab = "", ylab = "", main = species_name[i], x[,"meangd"], ylim = c(0,1))
  lines(x[,"upper95"], col = "red",lwd = 2, lty = "dashed")
  lines(x[,"lower95"], col = "green",lwd = 2, lty = "dashed")
  abline(h = 0.95, lty = "dotted", col = "blue")
  abline(v = min(which(x[,1] > 0.95)), lty = "dotted", col = "blue")
  if(i==6){
    legend(550, 4.27, xpd = NA, legend = leg.txt, fill = c("black", "red", "green"))
    # Label for x-axis
    mtext("Number of samples", side = 1, line=3, adj=-6,  cex = 1.5)
    # Label for y-axis
    mtext("Genetic Diversity", side = 2, line=3, adj=7.0, padj = -17.5, cex = 1.5)
  } else{
    if(i==12){
      legend(550, 4.27, xpd = NA, legend = leg.txt, fill = c("black", "red", "green"))
      # Label for x-axis
      mtext("Number of samples", side = 1, line=3, adj=-6, cex = 1.5)
      # Label for y-axis
      mtext("Genetic Diversity", side = 2, line=3, adj=7.0, padj = -17.5, cex = 1.5)
    } else {
      if(i==14){
        legend(550, 0.7, xpd = NA, legend = leg.txt, fill = c("black", "red", "green"))
        # Label for x-axis
        mtext("Number of samples", side = 1, line=3, adj=-6, cex = 1.5)
        # Label for y-axis
        mtext("Genetic Diversity", side = 2, line=3, adj=1.3, padj = -17.5, cex = 1.5)
      }
    }
  }
}
dev.off()


meanCI <- vector()
pdf(file = paste0(imagesDirectory, "14CIWidthplots.pdf"), width = 8.51, height = 7.27)
par(mfrow=c(2,3), omi=c(0,0.3,0,1.7))
for (i in 1:14) {
  x <- resultsArray[,,i]
  plot(xlab = "", ylab = "", main = species_name[i], x[,"ciwidth"], ylim = c(0,0.2), pch=16)
  meanCI[i] <- mean(x[,"ciwidth"])
  if(i==6){
    legend(550, 0.533, xpd = NA, legend = "test", fill = c("black", "red", "green"))
    # Label for x-axis
    mtext("Number of samples", side = 1, line=3, adj=2.75,  cex = 1.5)
    # Label for y-axis
    mtext("Confidence Interval Width", side = 2, line=3, adj=-3.5, padj = -24, cex = 1.5)
  }
}
dev.off()
meanCI

