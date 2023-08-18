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

2# This represents the amount of genetic diversity for a sample size of 
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
leg.txt1 = c("Total Mean Genetic Diversity")
legend(250, 0.7, legend = leg.txt1,
       fill = c("black"))

# graph of genetic diversity, replicates 1 - 3 are plotted along with a line showing the tmean gd
plot(xlab = "Sample size", ylab = "Genetic Diversity", final_quercus_results[,1,1], pch = 16)
points(final_quercus_results[,2,1], col="blue", pch = 16)
points(final_quercus_results[,3,1], col="green", pch = 16)
# recall this is the average and the above points are the individual replicates plotted on one graph
lines(meanRepValues, col="red", lwd=2)
leg.txt2 <- c("Replicate 1", "Replicate 2", "Replicate 3", "Total Mean of Genetic Diversity")
legend(200, 0.7, legend = leg.txt2,
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
leg.txt3 <- "Sample size"
legend(650, 250, legend = leg.txt3,
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
plot(xlab = "Sample Size", ylab = "Genetic Diversity", meanRepValues, ylim = c(0,1))
leg.txt4 = c("Total Mean Genetic Diversity", "95% Upper Limit", "95% Lower Limit")
lines(upper95, col = "red",lwd = 2, lty = "dashed")
lines(lower95, col = "green",lwd = 2, lty = "dashed")
abline(h = 0.95, lty = "dotted", col = "orange")
abline(v = min(which(meanRepValues > 0.95)), lty = "dotted", col = "orange")
legend(250, 0.7, legend = leg.txt4,
       fill = c("black", "red", "green"))


# testing for species 3
# meanRepValuess3 <- vector()
# upper953 <- vector()
# lower953 <- vector()
# for (n in 1:nrow(final_quercus_results)) {
#   meanRepValuess3[n] <- mean(final_quercus_results[n,,3])
#   upper953[n] <- quantile(final_quercus_results[n,,3],0.95)
#   lower953[n] <- quantile(final_quercus_results[n,,3],0.05)
# }
# 
# # add the lines to the legend, add an asymptote of 0.95 horizontal, and add the line at which sample size reaches the 0.95 benchmark vertically
# plot(xlab = "Sample Size", ylab = "Genetic Diversity", meanRepValuess3, ylim = c(0,1))
# leg.txt4 = c("Total Mean Genetic Diversity", "95% Upper Limit", "95% Lower Limit")
# lines(upper953, col = "red",lwd = 2, lty = "dashed")
# lines(lower953, col = "green",lwd = 2, lty = "dashed")
# abline(h = 0.95, lty = "dotted", col = "orange")
# abline(v = min(which(meanRepValuess3 > 0.95)), lty = "dotted", col = "orange")
# legend(250, 0.7, legend = leg.txt4,
#        fill = c("black", "red", "green"))

# # array1<-array(c(meanRepValues, upper95, lower95), dim =c(500,3,1))
# matrix1<-cbind(meanRepValues, upper95, lower95),500,3
# matrix3<-cbind(meanRepValuess3,upper953, lower953),500,3
# 
# resultsArray1 <- array(data=c(matrix1, matrix3), dim = c(500, 3, 2),
#                        dimnames=c("Samples","Stats","Species"))
# str(resultsArray1)


# now, we want to declare a higher dimension object for the 14 slices (spp.) of the array for the 3 vectors
resultsArray <- array(dim = c(500,4,14))
meansppvalue <- vector()
upper95spp <- vector()
lower95spp <- vector()

# for (n in 1:nrow(final_quercus_results)) {
#   browser()
#   meanRepValuess3[n] <- mean(final_quercus_results[n,,3])
#   upper953[n] <- quantile(final_quercus_results[n,,3],0.95)
#   lower953[n] <- quantile(final_quercus_results[n,,3],0.05)
# }

for (q in 1:14) {
  for (i in 1:nrow(final_quercus_results)) {
    # 
    meansppvalue[i] <- mean(final_quercus_results[i,,q])
    upper95spp[i] <- quantile(final_quercus_results[i,,q],0.95)
    lower95spp[i] <- quantile(final_quercus_results[i,,q],0.05)
    CIwidth <- upper95spp - lower95spp
    # TO DO: put the calculation of the CI width here!!!
    # You'll have to change the dimensions of the resultsArray
    # (4 columns, instead of 3)
  }
  # Bind vectors together into a matrix
  speciesMat <- cbind(meansppvalue, upper95spp, lower95spp, CIwidth)
  # Pass the matrix into a slot of the array
  resultsArray[,,q] <- speciesMat
  
}



str(resultsArray)

par(mfrow=c(3,2))
for (i in 1:14) {
  x <- resultsArray[,,i]
  plot(x[,1], ylim = c(0,1))
  lines(x[,2], col = "red",lwd = 2, lty = "dashed")
  lines(x[,3], col = "green",lwd = 2, lty = "dashed")
  abline(h = 0.95, lty = "dotted", col = "orange")
  abline(v = min(which(x[,1] > 0.95)), lty = "dotted", col = "orange")
  legend(250, 0.7, legend = leg.txt4,
         fill = c("black", "red", "green"))
}

# width of CI
upper95s <-resultsArray[1:500,2,]    
str(upper95s)
lower95s <- resultsArray[1:500,3,]
str(lower95s)
width <- upper95s - lower95s
str(width)

par(mfrow=c(2,3))
meanCI <- vector()
for (i in 1:14) {
  plot(x=1:nrow(resultsArray), y = width[,i], pch=16)
  meanCI[i] <- mean(width[,i])
}

pdf(file="14plots.pdf", width = 9, height = 7.5)
for (i in 1:14) {
  plot(x=1:nrow(resultsArray), y = width[,i], pch=16)
  meanCI[i] <- mean(width[,i])
}
dev.off()

png(file="14plots.png", width = 480, height = 480)
for (i in 1:14) {
  plot(x=1:nrow(resultsArray), y = width[,i], pch=16)
  meanCI[i] <- mean(width[,i])
}
dev.off()

meanCI
matrix(meanCI, nrow=14)

matrix(width[25,1:14])
matrix(width[50,1:14])
matrix(width[100,1:14])
matrix(width[200,1:14])

# plot(x=1:nrow(resultsArray), y=width[,1], pch=16)

# for (i in 1:14) {
#   plot(resultsArray[,1,i])
#   lines(resultsArray[,2,i])
#   lines(resultsArray[,3,i])
# }


# for (q in 1:14) {
#   for (n in 1:nrow(final_quercus_results)) {
#     resultsArray[n,,q] <- c(matrix(cbind(mean(final_quercus_results[n,,q]), 
#     quantile(final_quercus_results[n,,q],0.95),
#     quantile(final_quercus_results[n,,q],0.05), dim = c(n,1000,q))))
#   }
# }

# for (q in 1:14) {
#   for (n in 1:nrow(final_quercus_results)) {
#     matrix[n] <- matrix(cbind(mean(final_quercus_results[n,,q]), 
#                  quantile(final_quercus_results[n,,q],0.95), 
#                  quantile(final_quercus_results[n,,q],0.05),500,3))
#     resultsArray[,,q] <- array(c(matrix[n]), dim = c(500,3,q))
#   }
# }


# meansppvalue
# 
# 
# resultsArray
# 
# mean95 <- list()
# for(a in 1:nrow(final_quercus_results)){
#   mean95[[a]] <- mean(final_quercus_results[a,,1])
# }
# # creates a vector from the list so we can find min sample size
# unlist(mean95)
# # function to more efficiently find min 0.95 value?
# 
# abline(v = 191, lty = "dotted", col = "orange")
# 
# 
# abline(v = mean(min_samp95), lty ="dotted", col = "orange")