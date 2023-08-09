##########################################
# 2023/08/01 Exploring Resampling Dataset#
##########################################

setwd("C:/Users/gsalas/Documents/resampling_CIs/Code/")

load("Datasets/quercus_final_results_orig.Rdata")

ls()

dim(final_quercus_results)

final_quercus_results[,,1]

# This represents samples in a matrix of 1-100 for one replicate of one species.
final_quercus_results[1:100,1,1]

# This represents one randomly selected sample replicated one hundred times for one species
final_quercus_results[1,1:100,1]

# This represents one hundred randomly selected samples
final_quercus_results[1:100,1,2]
 
# this plot displays
plot(final_quercus_results[,1,1],col="blue"); points(final_quercus_results[,1,2],col="red")


plot(final_quercus_results[1,,2])

# this gives us the mean genetic diversity values for a range of sample sizes, 
# range of replicates, and species' along the rows(lines 30 - 35 show how Kaylee did her work)
sp <- 1
# we made a for loop alternative for the apply() function
apply(final_quercus_results[,,1],1,mean)
dim(final_quercus_results)
plot(final_quercus_results[,1,sp])
lines(apply(final_quercus_results[,,sp],1,mean),col="red",lwd=2)

# (loop alt) assign a variable to the vector to eventually capture all the means across the replicates
e <- vector()

# nrow is soft code that detects how many rows there are in an object
# we use it to get the amount of genetic diversity across replicates 
# we then set an index for e in order to get a mean.
for(i in 1:nrow(final_quercus_results)){
  e[i] <- mean(final_quercus_results[i,,1])
  }
# mean genetic diversity across all replicates for species 1
e
# graph of e 
plot(e)

plot(final_quercus_results[,1,1])
points(final_quercus_results[,2,1], col="blue")
points(final_quercus_results[,3,1], col="green")
# recall this is the average and the above points are the individual replicates plotted on one graph
lines(e, col="red", lwd=2)

# get the 95% CI of plot, but first go thru these ideas
# IDEA 1
# this gives us the position in the vector of the min 
# genetic diversity value greater than 0.95
min(which(final_quercus_results[,1,sp]>0.95))
sp<-1; min_samp95<-vector(length = 1000)
for (r in 1:1000) {
  min_samp95[r]<-min(which(final_quercus_results[,r,sp]>0.95))
  
}
#this gives mean across reps of the first individual to cross 95%
boxplot(min_samp95)

# distribution of min sample sizes across replicates
# not what we quite want however...
plot(min_samp95)

# IDEA 2
p<-1
# this shows the genetic diversity value in 95 percentile of the values across the replicates for a sample size of one for species one
quantile(final_quercus_results[p,,1],0.95)
