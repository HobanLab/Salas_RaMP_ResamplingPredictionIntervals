rm(list = ls())
library(strataG)
library(adegenet)
# %%% STRATAG APPROACH ----
# Specify the filepath to the Arlequin file in question
arpFilepath <- "C:/Users/gsalas/Desktop/fsc28_win64/MSAT_04pop_migLow/MSAT_04pop_migLow_1_5.arp"
arpFiles <- list.files(path = arpFilepath, pattern = ".arp", full.names =  TRUE)


# Read in the Arlequin file using read.arlequin
arpObject <- arlequinRead(arpFilepath)
# Convert the Arlequin object to a gtypes object
gTypeObject <- arp2gtypes(arpObject)
# Convert the gtypes to a genind object
genindObject <- gtypes2genind(gTypeObject)

genindObject@tab

arpObjectList <- list()
for (i in 1:length(arpFiles)) {
  arpObjectList[i] <- arlequinRead(arpFiles[i])
}

arpObjectList[1]
