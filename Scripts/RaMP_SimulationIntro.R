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
# Double the population size
rm(list = ls())
setwd("C:/Users/gsalas/Desktop/fsc28_win64/RaMP_SimulationIntro_DemoParFiles/")

deme0 <- fscDeme(deme.size = 600, sample.size = 600)

demes <- fscSettingsDemes(deme0,deme0,deme0,deme0, ploidy = 2)

events0 <- fscEvent(event.time = 5e4, source = 0, sink = 0, prop.migrants = 0, new.size = 1, new.growth = 0, migr.mat = 1)
events1 <- fscEvent(event.time = 5e4, source = 1, sink = 0, prop.migrants = 1, new.size = 1, new.growth = 0, migr.mat = 1)
events2 <- fscEvent(event.time = 5e4, source = 2, sink = 0, prop.migrants = 1, new.size = 1, new.growth = 0, migr.mat = 1)
events3 <- fscEvent(event.time = 5e4, source = 3, sink = 0, prop.migrants = 1, new.size = 1, new.growth = 0, migr.mat = 1)

events <- fscSettingsEvents(events0, events1, events2, events3)

migmat1 <- matrix(data = c(0, 0.001, 0.001, 0.001,
                0.001, 0, 0.001, 0.001,
                0.001, 0.001, 0, 0.001,
                0.001, 0.001, 0.001, 0), nrow = 4, ncol = 4)
migmat2 <- matrix(data = c(0, 0, 0, 0,
                           0, 0, 0, 0,
                           0, 0, 0, 0,
                           0, 0, 0, 0), nrow = 4, ncol = 4)

migration <- fscSettingsMigration(migmat1,migmat2)

msats <- fscBlock_microsat(num.loci = 1, recomb.rate = 0, mut.rate = 1e-3)

genetics <- fscSettingsGenetics(msats, num.chrom = 20)

updatedMSAT_04pop_migLowdoubledPop.params <- fscWrite(demes = demes, genetics = genetics, events = events, migration = migration, label ="MSAT_04pop_migLow_doubledPop", use.wd = TRUE)

# Add five more loci
rm(list = ls())
setwd("C:/Users/gsalas/Desktop/fsc28_win64/RaMP_SimulationIntro_DemoParFiles/")

deme0 <- fscDeme(deme.size = 300, sample.size = 300)

demes <- fscSettingsDemes(deme0,deme0,deme0,deme0, ploidy = 2)

events0 <- fscEvent(event.time = 5e4, source = 0, sink = 0, prop.migrants = 0, new.size = 1, new.growth = 0, migr.mat = 1)
events1 <- fscEvent(event.time = 5e4, source = 1, sink = 0, prop.migrants = 1, new.size = 1, new.growth = 0, migr.mat = 1)
events2 <- fscEvent(event.time = 5e4, source = 2, sink = 0, prop.migrants = 1, new.size = 1, new.growth = 0, migr.mat = 1)
events3 <- fscEvent(event.time = 5e4, source = 3, sink = 0, prop.migrants = 1, new.size = 1, new.growth = 0, migr.mat = 1)

events <- fscSettingsEvents(events0, events1, events2, events3)

migmat1 <- matrix(data = c(0, 0.001, 0.001, 0.001,
                           0.001, 0, 0.001, 0.001,
                           0.001, 0.001, 0, 0.001,
                           0.001, 0.001, 0.001, 0), nrow = 4, ncol = 4)
migmat2 <- matrix(data = c(0, 0, 0, 0,
                           0, 0, 0, 0,
                           0, 0, 0, 0,
                           0, 0, 0, 0), nrow = 4, ncol = 4)

migration <- fscSettingsMigration(migmat1,migmat2)

msats <- fscBlock_microsat(num.loci = 1, recomb.rate = 0, mut.rate = 1e-3)

genetics <- fscSettingsGenetics(msats, num.chrom = 25)

MSAT_04pop_migLow5Loci.params <- fscWrite(demes = demes, genetics = genetics, events = events, migration = migration, label ="MSAT_04pop_migLow_5Loci", use.wd = TRUE)