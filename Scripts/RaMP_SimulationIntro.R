rm(list = ls())
library(strataG)
library(adegenet)
# %%%STRATAG APPROACH ----
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

# WRITING THE UPDATED PAR FILE ----
setwd("C:/Users/gsalas/Desktop/fsc28_win64/RaMP_SimulationIntro_DemoParFiles/")

deme0 <- fscDeme(deme.size = 600, sample.size = 4)

demes <- fscSettingsDemes(deme0, ploidy = 2)

fscEvent(event.time = 5e-4, source = 0, sink = 0, prop.migrants = 1, new.size = 1, new.growth = 0, migr.mat = 1)

msats <- fscBlock_microsat(num.loci = 1, recomb.rate = 0, mut.rate = 1e-3)


genetics <- fscSettingsGenetics(msats, num.chrom = 25)

updatedMSAT_04pop_migLow.params <- fscWrite(demes = demes, genetics = genetics, label ="updatedMSAT_04pop_migLow", use.wd = TRUE)
