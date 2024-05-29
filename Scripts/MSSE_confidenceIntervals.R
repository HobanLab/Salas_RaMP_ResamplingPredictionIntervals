###########################################
# 2023/08/01 Exploring Resampling Dataset #
###########################################

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
# abline(v = min(which(upper95 > 0.95)), lty = "dotted", col="red")
# abline(v = min(which(lower95 > 0.95)), lty = "dotted", col="green")
legend(250, 0.7, legend = leg.txt,
       fill = c("black", "red", "green"))


# now, we want to declare a higher dimension object for the 14 slices (spp.) of the array for the 3 vectors
# create a vector of the 14 species names abbreviated
species_name <- c("QUAC","QUAR","QUAU", "QUBO","QUCA","QUCE","QUEN","QUGE","QUGR","QUHA","QUHI","QUOG","QUPA", "QUTO")
# create an array to store the outputs for each species
calcs_for_quercus14 <- array(dim = c(500,4,14))
# dimnames will create descriptions you can use instead of numbers when subsetting calcs_for_quercus14
dimnames(calcs_for_quercus14)<-list(paste0("sample",1:500),c("meanRepValues","upper95","lower95","ciwidth"), species_name)
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
  # Pass the matrices into a slot of the array
  calcs_for_quercus14[,,q] <- speciesMat
}
str(calcs_for_quercus14)
saveRDS(calcs_for_quercus14, "C:/Users/gsalas/Documents/resampling_CIs/Code/Outputs/calcs_for_quercus14.Rdata")
x <- calcs_for_quercus14
# upper95 <- vector()
# lower95 <- vector()
CImatrix <- matrix(nrow = 14, ncol = 4)
rownames(CImatrix) <- dimnames(x)[[3]]
colnames(CImatrix) <- c("MSSE", "Lower", "Upper", "CIWidth")
for (q in 1:14) {
  min_95totavg <-(min(which(x[,,q][,"meanRepValues"]>0.95)))
  lowerbound <- (min(which(x[,,q][,"lower95"] > 0.95)))
  upperbound <- (min(which(x[,,q][,"upper95"] > 0.95)))
  ciwidth <- lowerbound - upperbound
  CImatrix[q,] <- cbind(min_95totavg, upperbound, lowerbound, ciwidth)
}
CImatrix[is.infinite(CImatrix)]<-NA  
getwd()
write.csv(CImatrix, file = "Outputs/Quercus14_CI_values.csv")

# imagesDirectory is an object that is a designated file path you use to paste images of plots
imagesDirectory <- "C:/Users/gsalas/Documents/resampling_CIs/Code/Images/"
# create pdf using pdf(), specify the file path by pasting images directory with the .pdf title of your plot and specify dimensions
pdf(file=paste0(imagesDirectory, "14CIPlots.pdf"), width = 8.5, height = 11)
# create pdf using pdf(), specify the file path by pasting images directory with the .pdf title o
par(mfrow=c(3,2), omi=c(0.8,0.3,0,0))
for (i in 1:14) {
  x <- calcs_for_quercus14[,,i]
  # species 9 and species 11 are passed through an ifelse() function
  # to be subset by 
  if(i==9 | i==11){
    plot(xlab = "", ylab = "", main = species_name[i], x[-which(x[,"meanRepValues"]==0),1], xlim = c(0,525), ylim = c(0,1), pch=16)
    lines(x[,"upper95"][c(1:length(x[-which(x[,"meanRepValues"]==0),1]))], col = "red",lwd = 2, lty = "dashed")
    lines(x[,"lower95"][c(1:length(x[-which(x[,"meanRepValues"]==0),1]))], col = "green",lwd = 2, lty = "dashed")
    abline(h = 0.95, lty = "dotted", col = "blue")
    abline(v = min(which(x[,1] > 0.95)), lty = "dotted", col = "orange")
  } else{
    plot(xlab = "", ylab = "", main = species_name[i], x[,"meanRepValues"], xlim = c(0,525), ylim = c(0,1), pch=16)
    lines(x[,"upper95"], col = "red",lwd = 2, lty = "dashed")
    lines(x[,"lower95"], col = "blue",lwd = 2, lty = "dashed")
    abline(h = 0.95, lty = "dotted", col = "blue")
    abline(v = min(which(x[,1] > 0.95)), lty = "dotted", col = "brown")
    if (i==6|i==12){
      legend(-250, -0.4, xpd = NA, legend = leg.txt, fill = c("black", "red", "green"))
      # Label for x-axis
      mtext("Number of samples", side = 1, line=3, adj=-1.5,  cex = 1.5)
      # Label for y-axis (adj adjusts the vertical position of text)
      mtext("Genetic Diversity", side = 2, line=11.0, adj=11.75, padj = -17.5, cex = 1.5)
    } else{
      if(i==14){
        legend(-225, -0.4, xpd = NA, legend = leg.txt, fill = c("black", "red", "green"))
        # Label for x-axis
        mtext("Number of samples", side = 1, line=3, adj=-1.5, cex = 1.5)
        # Label for y-axis
        mtext("Genetic Diversity", side = 2, line=3, adj=0.4, padj = -22.5, cex = 1.5)
      }
    }
  }
}
dev.off()

# creating a vector to store the average confidence interval width for each species 
meanCI <- vector()
pdf(file = paste0(imagesDirectory, "14CIWidthplots.pdf"), width = 8.51, height = 7.27)
par(mfrow=c(2,3), omi=c(0.8,0.3,0,0))
for (i in 1:14) {
  x <- calcs_for_quercus14[,,i]
  meanCI[i] <- mean(x[,"ciwidth"])
  if(i==9|i==11){
    plot(xlab = "", ylab = "", main = species_name[i], x[,"ciwidth"][c(1:length(x[-which(x[,"ciwidth"]==0),1]))], xlim = c(0,525), ylim = c(0,0.2), pch=16)
  } else{
    plot(xlab = "", ylab = "", main = species_name[i], x[,"ciwidth"], xlim = c(0,525), ylim = c(0,0.2), pch=16)
    if(i==6|i==12){
      legend(-750, -0.09, xpd = NA, legend = "Calculated CI width", fill = c("black", "red", "green"))
      # Label for x-axis
      mtext("Number of samples", side = 1, line=3, adj=11.5,  cex = 1.5)
      # Label for y-axis
      mtext("Confidence Interval Width", side = 2, line=12.0, adj=-1.0, padj = -24, cex = 1.5)
    } else{
      if(i==14){
        legend(-300, -0.08, xpd = NA, legend = "Calculated CI width", fill = c("black", "red", "green"))
        # Label for x-axis
        mtext("Number of samples", side = 1, line=3, adj=11.5,  cex = 1.5)
        # Label for y-axis
        mtext("Confidence Interval Width", side = 2, line=0, adj=0.5, padj = -17, cex = 1.5)
      }
    }
  }
}
dev.off()
meanCI

#table of all the values at given points
sample25 <- calcs_for_quercus14[25,4,1:14]
sample50 <- calcs_for_quercus14[50,4,1:14]
sample100 <- calcs_for_quercus14[100,4,1:14]
write.csv(cbind(sample25, sample50, sample100, meanCI), file = "14CIWidths.csv")


pdf(file = paste0(imagesDirectory, "14CIWidthplotshigh.pdf"), width = 8.51, height = 7.27)
par(mfrow=c(2,3), omi=c(0.8,0.3,0,0))
for (i in 1:14) {
  x <- calcs_for_quercus14[,,i]
  plot(xlab = "", ylab = "", main = species_name[i], x[,"ciwidth"], xlim = c(0,50), ylim = c(0,0.2), pch=16)
  box(col="black",
      which = "figure")
  meanCI[i] <- mean(x[,"ciwidth"])
  if(i==6|i==12){
    legend(-70, -0.08, xpd = NA, legend = "Calculated CI width", fill = c("black", "red", "green"))
    # Label for x-axis
    mtext("Number of samples", side = 1, line=3, adj=11.5,  cex = 1.5)
    # Label for y-axis
    mtext("Confidence Interval Width", side = 2, line=11.0, adj=-1.0, padj = -24, cex = 1.5)
  }else{
    if(i==14){
      # legend(-300, -0.08, xpd = NA, legend = "Calculated CI width", fill = c("black", "red", "green"))
      #   # Label for x-axis
      #   mtext("Number of samples", side = 1, line=3, adj=11.5,  cex = 1.5)
      #   # Label for y-axis
      #   mtext("Confidence Interval Width", side = 2, line=0, adj=0.5, padj = -17, cex = 1.5)
      legend(-30, -0.08, xpd = NA, legend = "Calculated CI width", fill = c("black", "red", "green"))
      # Label for x-axis
      mtext("Number of samples", side = 1, line=3, adj=5,  cex = 1.5)
      # Label for y-axis
      mtext("Confidence Interval Width", side = 2, line=0, adj=0.5, padj = -17, cex = 1.5)
    }
  }
}
dev.off()


pdf(file = paste0(imagesDirectory, "14CIWidthplotslow.pdf"), width = 8.51, height = 7.27)
par(mfrow=c(2,3), omi=c(0.8,0.3,0,0))
for (i in 1:14) {
  x <- calcs_for_quercus14[,,i]
  plot(xlab = "", ylab = "", main = species_name[i], x[,"ciwidth"], xlim = c(400,500), ylim = c(0,0.2), pch=16)
  meanCI[i] <- mean(x[,"ciwidth"])
  if(i==6|i==12){
    legend(-250, -0.08, xpd = NA, legend = "test", fill = c("black", "red", "green"))
    # Label for x-axis
    mtext("Number of samples", side = 1, line=3, adj=11.5,  cex = 1.5)
    # Label for y-axis
    mtext("Confidence Interval Width", side = 2, line=12.0, adj=-1.0, padj = -24, cex = 1.5)
  }else{
    if(i==14){
      legend(-200, -0.08, xpd = NA, legend = "test", fill = c("black", "red", "green"))
      # Label for x-axis
      mtext("Number of samples", side = 1, line=3, adj=11.5,  cex = 1.5)
      # Label for y-axis
      mtext("Confidence Interval Width", side = 2, line=12.0, adj=-1.0, padj = -24, cex = 1.5)
    }
  }
}
dev.off()


# 04/22/2022
# Three species confidence interval plots
#
# imagesDirectory is an object that is a designated file path you use to paste images of plots
imagesDirectory <- "C:/Users/gsalas/Documents/resampling_CIs/Code/Images/"
# create pdf using pdf(), specify the file path by pasting images directory with the .pdf title of your plot and specify dimensions
pdf(file=paste0(imagesDirectory, "CPCposter3CIplots.pdf"), width = 11, height = 4)
# png(file=paste0(imagesDirectory, "test.png"), units = "in", width = 11, height = 4, res = 300)
# create pdf using pdf(), specify the file path by pasting images directory with the .pdf title
# x <- calcs_for_quercus14
# par(mfrow=c(1,2), omi=c(0,0,0,0))
# plot(xlab = "", ylab = "", x[,,1][,"meanRepValues"], xlim = c(0,525), ylim = c(0,1), pch=16, col='blue')
# abline(h = 0.95, lty = "dotted", col = "pink", lwd = 2)
# lines(x[,,1][,"upper95"], col = "blue",lwd = 2, lty = "dashed")
# lines(x[,,1][,"lower95"], col = "blue",lwd = 2, lty = "dashed")
# points(xlab = "", ylab = "", x[,,2][,"meanRepValues"], xlim = c(0,525), ylim = c(0,1), pch=16, col = 'red')
# lines(x[,,2][,"upper95"], col = "red",lwd = 2, lty = "dashed")
# lines(x[,,2][,"lower95"], col = "red",lwd = 2, lty = "dashed")
# points(xlab = "", ylab = "", x[,,7][,"meanRepValues"], xlim = c(0,525), ylim = c(0,1), pch=16, col = 'orange')
# lines(x[,,7][,"upper95"], col = "orange",lwd = 2, lty = "dashed")
# lines(x[,,7][,"lower95"], col = "orange",lwd = 2, lty = "dashed")
# plot(xlab = "", ylab = "", x[,,1][,"meanRepValues"], xlim = c(0,525), ylim = c(0.6,1), pch=16, col='blue')
# abline(h = 0.95, lty = "dotted", col = "pink", lwd = 2)
# lines(x[,,1][,"upper95"], col = "blue",lwd = 2, lty = "dashed")
# lines(x[,,1][,"lower95"], col = "blue",lwd = 2, lty = "dashed")
# points(xlab = "", ylab = "", x[,,2][,"meanRepValues"], pch=16, col = 'red')
# lines(x[,,2][,"upper95"], col = "red",lwd = 2, lty = "dashed")
# lines(x[,,2][,"lower95"], col = "red",lwd = 2, lty = "dashed")
# points(xlab = "", ylab = "", x[,,7][,"meanRepValues"], pch=16, col = 'orange')
# lines(x[,,7][,"upper95"], col = "orange",lwd = 2, lty = "dashed")
# lines(x[,,7][,"lower95"], col = "orange",lwd = 2, lty = "dashed")
# # Label for x-axis
# mtext("Proportion of allelic representation", side = 2, line = -1.5, outer=TRUE)
# # # Label for y-axis
# mtext("Number of randomly sampled individuals", side = 1, line = -2, outer = TRUE)
# legend("bottomright", legend = c("Q. acerifolia", "Q. boyntonii", "Q. engelmannii", "Confidence intervals"),
#        xpd = NA,
#        cex = 3,
#        col = c("blue", 'red', 'orange', "black"),
#        lty = c(NA,NA,NA,"dashed"),
#        pch = c(16,16,16,NA))
# dev.off()

# plot.new()
imagesDirectory <- "C:/Users/gsalas/Documents/resampling_CIs/Code/Images/"
png(file=paste0(imagesDirectory, "3CIplotCPC2024.png"), units = "in", width = 13.05, height = 5.57, res = 300)
dev.new()
x <- calcs_for_quercus14
min_95totavg<-(min(which(x[,,1] > 0.95)))
upperbound <- (min(which(x[,,1][,"upper95"] > 0.95)))
lowerbound <- (min(which(x[,,1][,"lower95"] > 0.95)))
par(mfrow=c(1,2), omi=c(0,0,0,0))
# par("mar" = c(0,0,0,0))
plot(xlab = "", ylab = "", x[,,1][,"meanRepValues"], xlim = c(0,525), ylim = c(0,1), pch=16, col='blue')
lines(x[,,1][,"upper95"], col = "blue",lwd = 2, lty = "dashed")
lines(x[,,1][,"lower95"], col = "blue",lwd = 2, lty = "dashed")
points(xlab = "", ylab = "", x[,,2][,"meanRepValues"], xlim = c(0,525), ylim = c(0,1), pch=16, col = 'red')
lines(x[,,2][,"upper95"], col = "red",lwd = 2, lty = "dashed")
lines(x[,,2][,"lower95"], col = "red",lwd = 2, lty = "dashed")
points(xlab = "", ylab = "", x[,,7][,"meanRepValues"], xlim = c(0,525), ylim = c(0,1), pch=16, col = 'orange')
lines(x[,,7][,"upper95"], col = "orange",lwd = 2, lty = "dashed")
lines(x[,,7][,"lower95"], col = "orange",lwd = 2, lty = "dashed")
abline(h = 0.95, lty = "dotted", col = "pink", lwd = 2)
legend(x = 300, y = 0.25, legend = c("Q. acerifolia", "Q. boyntonii", "Q. engelmannii", "Confidence intervals"),
       xpd = NA,
       cex = 1,
       col = c("blue", 'red', 'orange', "black"),
       lty = c(NA,NA,NA,"dashed"),
       pch = c(16,16,16,NA))
plot(xlab = "", ylab = "", x[,,1][,"meanRepValues"], xlim = c(0,525), ylim = c(0.6,1), pch=16, col='blue')
lines(x[,,1][,"upper95"], col = "blue",lwd = 2, lty = "dashed")
lines(x[,,1][,"lower95"], col = "blue",lwd = 2, lty = "dashed")
points(xlab = "", ylab = "", x[,,2][,"meanRepValues"], pch=16, col = 'red')
lines(x[,,2][,"upper95"], col = "red",lwd = 2, lty = "dashed")
lines(x[,,2][,"lower95"], col = "red",lwd = 2, lty = "dashed")
points(xlab = "", ylab = "", x[,,7][,"meanRepValues"], pch=16, col = 'orange')
lines(x[,,7][,"upper95"], col = "orange",lwd = 2, lty = "dashed")
lines(x[,,7][,"lower95"], col = "orange",lwd = 2, lty = "dashed")
abline(h = 0.95, lty = "dotted", col = "pink", lwd = 2)
text(x = 300, y = 0.7, paste("Q. acerifolia MSSE = ", min_95totavg, "individuals ", "CI [",upperbound,', ',lowerbound,']'), col="blue", cex = 1)
# Label for x-axis
mtext("Proportion of allelic representation", cex = 1,side = 2, line = -2, outer=TRUE)
# # Label for y-axis
mtext("Number of randomly sampled individuals", cex = 1, side = 1, line = -2, outer = TRUE)
dev.off()