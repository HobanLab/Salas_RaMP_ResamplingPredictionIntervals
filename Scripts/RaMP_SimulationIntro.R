rm(list = ls())
library(strataG)
library(adegenet)
# %%% STRATAG APPROACH ----
# Specify the filepath to the Arlequin file in question
arpFilepath <- "C:/Users/gsalas/Desktop/fsc28_win64/MSAT_04pop_migLow/"
arpFiles <- list.files(path = arpFilepath, pattern = ".arp", full.names =  TRUE)
arpFileslist <- list(arpFiles)
# Read in the Arlequin file using read.arlequin
arpObject <- lapply(arpFileslist[[1]], arlequinRead)
# Convert the Arlequin object to a gtypes object
gTypeObject <- list()
genindObject <- list()
locivec <- vector()
meanLocivec <- vector()
Allelevec <- vector()
meanAllelevec <- vector()
for (i in 1:length(arpObject)) {
  gTypeObject[[i]] <- arp2gtypes(arpObject[[i]])
  genindObject[[i]] <- gtypes2genind(gTypeObject[[i]])
  locivec[i] <- c(sum(nLoc(genindObject[[i]]))) 
  meanLocivec <- mean(locivec)
  Allelevec[i] <- c(sum(nAll(genindObject[[i]])))
  meanAllelevec <- mean(Allelevec)
}
meanLocivec
meanAllelevec