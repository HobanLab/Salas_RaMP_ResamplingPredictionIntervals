################################################
# GDS 2024/01/10 QUAC QUBO Loci Bootstrapping #
################################################
setwd("C:/Users/gsalas/Documents/resampling_CIs/Code/Datasets")

datasets <- list.files(path = "LociBootstrapping_Datasets", pattern = "Wild", full.names = TRUE)

gm_list <- list(
  readRDS(datasets[1]),
  readRDS(datasets[2]),
  readRDS(datasets[3]),
  readRDS(datasets[4])
)
QUAC.MSAT.Wild.genind <- gm_list[[1]]
QUAC.SNP.Wild.genind <- gm_list[[2]]
QUBO.MSAT.Wild.genind <- gm_list[[3]]
QUBO.SNP.Wild.genind <- gm_list[[4]]
library(adegenet)
## MSAT  ##
# 5 loci
# Loop 25 times for 25 sets (5 loci in one set) of randomly selected loci
MSAT_5loc_function <- function(test_genind){
  resamp_category5loc <- array(dim = c(nrow(test_genind@tab),4,25))
  for (i in 1:25) {
    samp_5loc <- sample(locNames(test_genind), size = (length(locNames(test_genind))*(1/3)), replace = FALSE)
    MSAT.5loc.WILD.genind <- test_genind[, loc = samp_5loc]
    wildSamp_5loc <- MSAT.5loc.WILD.genind@tab
    wildSamp_5loc <- wildSamp_5loc[, which(colSums(wildSamp_5loc, na.rm = TRUE) != 0)]
    wildComplete <- colSums(wildSamp_5loc, na.rm = TRUE) / (nrow(wildSamp_5loc) * 2)
    wildSubset <- wildComplete[wildComplete > 0]
    total <- vector(length = nrow(wildSamp_5loc))
    common <- vector(length = nrow(wildSamp_5loc))
    lowfreq <- vector(length = nrow(wildSamp_5loc))
    rare <- vector(length = nrow(wildSamp_5loc))
    # Loop for each sample size
    for (j in 1:nrow(wildSamp_5loc)) {
      samp <- sample(nrow(wildSamp_5loc), size = j, replace = FALSE)
      samp <- wildSamp_5loc[samp,]
      # Measure the proportion of allelic representation in that sample
      if (j == 1) {
        total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(samp > 0)))]) / length(names(wildSubset))
        common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.1))
        lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
        rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(samp > 0)))) / length(which(wildSubset < 0.01))
      } else {
        total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))]) / length(names(wildSubset))
        common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.1))
        lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
        rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset < 0.01))
      }
      categorymat_5loc <- cbind(total, common, lowfreq, rare)
      category <- colnames(categorymat_5loc)
      dimnames(resamp_category5loc) <- list(paste0("sample ", 1:nrow(categorymat_5loc)), category, paste0("replicate ", 1:25))
      # Store the results for the current replicate in resamp_category5loc
      resamp_category5loc[, , i] <- categorymat_5loc
    }
  }
  return(resamp_category5loc)
}
MSAT_5loc_function(QUAC.MSAT.Wild.genind)
MSAT_5loc_function(QUBO.MSAT.Wild.genind)

# 10 loci
MSAT_10loc_function <- function(test_genind){
  resamp_category5loc <- array(dim = c(nrow(test_genind@tab),4,25))
  for (i in 1:25) {
    samp_5loc <- sample(locNames(test_genind), size = (length(locNames(test_genind))*(2/3)), replace = FALSE)
    MSAT.5loc.WILD.genind <- test_genind[, loc = samp_5loc]
    wildSamp_5loc <- MSAT.5loc.WILD.genind@tab
    wildSamp_5loc <- wildSamp_5loc[, which(colSums(wildSamp_5loc, na.rm = TRUE) != 0)]
    wildComplete <- colSums(wildSamp_5loc, na.rm = TRUE) / (nrow(wildSamp_5loc) * 2)
    wildSubset <- wildComplete[wildComplete > 0]
    total <- vector(length = nrow(wildSamp_5loc))
    common <- vector(length = nrow(wildSamp_5loc))
    lowfreq <- vector(length = nrow(wildSamp_5loc))
    rare <- vector(length = nrow(wildSamp_5loc))
    # Loop for each sample size
    for (j in 1:nrow(wildSamp_5loc)) {
      samp <- sample(nrow(wildSamp_5loc), size = j, replace = FALSE)
      samp <- wildSamp_5loc[samp,]
      # Measure the proportion of allelic representation in that sample
      if (j == 1) {
        total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(samp > 0)))]) / length(names(wildSubset))
        common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.1))
        lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
        rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(samp > 0)))) / length(which(wildSubset < 0.01))
      } else {
        total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))]) / length(names(wildSubset))
        common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.1))
        lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
        rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset < 0.01))
      }
      categorymat_5loc <- cbind(total, common, lowfreq, rare)
      category <- colnames(categorymat_5loc)
      dimnames(resamp_category5loc) <- list(paste0("sample ", 1:nrow(categorymat_5loc)), category, paste0("replicate ", 1:25))
      # Store the results for the current replicate in resamp_category5loc
      resamp_category5loc[, , i] <- categorymat_5loc
    }
  }
  return(resamp_category5loc)
}
MSAT_10loc_function(QUAC.MSAT.Wild.genind)
MSAT_10loc_function(QUBO.MSAT.Wild.genind)

# total loci
MSAT_totalloc_function <- function(test_genind){
  resamp_category5loc <- array(dim = c(nrow(test_genind@tab),4,25))
  for (i in 1:25) {
    samp_5loc <- sample(locNames(test_genind), size = (length(locNames(test_genind))), replace = FALSE)
    MSAT.5loc.WILD.genind <- test_genind[, loc = samp_5loc]
    wildSamp_5loc <- MSAT.5loc.WILD.genind@tab
    wildSamp_5loc <- wildSamp_5loc[, which(colSums(wildSamp_5loc, na.rm = TRUE) != 0)]
    wildComplete <- colSums(wildSamp_5loc, na.rm = TRUE) / (nrow(wildSamp_5loc) * 2)
    wildSubset <- wildComplete[wildComplete > 0]
    total <- vector(length = nrow(wildSamp_5loc))
    common <- vector(length = nrow(wildSamp_5loc))
    lowfreq <- vector(length = nrow(wildSamp_5loc))
    rare <- vector(length = nrow(wildSamp_5loc))
    # Loop for each sample size
    for (j in 1:nrow(wildSamp_5loc)) {
      samp <- sample(nrow(wildSamp_5loc), size = j, replace = FALSE)
      samp <- wildSamp_5loc[samp,]
      # Measure the proportion of allelic representation in that sample
      if (j == 1) {
        total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(samp > 0)))]) / length(names(wildSubset))
        common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.1))
        lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
        rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(samp > 0)))) / length(which(wildSubset < 0.01))
      } else {
        total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))]) / length(names(wildSubset))
        common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.1))
        lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
        rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset < 0.01))
      }
      categorymat_5loc <- cbind(total, common, lowfreq, rare)
      category <- colnames(categorymat_5loc)
      dimnames(resamp_category5loc) <- list(paste0("sample ", 1:nrow(categorymat_5loc)), category, paste0("replicate ", 1:25))
      # Store the results for the current replicate in resamp_category5loc
      resamp_category5loc[, , i] <- categorymat_5loc
    }
  }
  return(resamp_category5loc)
}
MSAT_totalloc_function(QUAC.MSAT.Wild.genind)
MSAT_totalloc_function(QUBO.MSAT.Wild.genind)

## SNP QUAC ##
# 5k loci
SNP_5kloc_function <- function(test_genind){
  resamp_category5loc <- array(dim = c(nrow(test_genind@tab),4,25))
  for (i in 1:25) {
    samp_5loc <- sample(locNames(test_genind), size = 5000, replace = FALSE)
    MSAT.5loc.WILD.genind <- test_genind[, loc = samp_5loc]
    wildSamp_5loc <- MSAT.5loc.WILD.genind@tab
    wildSamp_5loc <- wildSamp_5loc[, which(colSums(wildSamp_5loc, na.rm = TRUE) != 0)]
    wildComplete <- colSums(wildSamp_5loc, na.rm = TRUE) / (nrow(wildSamp_5loc) * 2)
    wildSubset <- wildComplete[wildComplete > 0]
    total <- vector(length = nrow(wildSamp_5loc))
    common <- vector(length = nrow(wildSamp_5loc))
    lowfreq <- vector(length = nrow(wildSamp_5loc))
    rare <- vector(length = nrow(wildSamp_5loc))
    # Loop for each sample size
    for (j in 1:nrow(wildSamp_5loc)) {
      samp <- sample(nrow(wildSamp_5loc), size = j, replace = FALSE)
      samp <- wildSamp_5loc[samp,]
      # Measure the proportion of allelic representation in that sample
      if (j == 1) {
        total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(samp > 0)))]) / length(names(wildSubset))
        common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.1))
        lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
        rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(samp > 0)))) / length(which(wildSubset < 0.01))
      } else {
        total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))]) / length(names(wildSubset))
        common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.1))
        lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
        rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset < 0.01))
      }
      categorymat_5loc <- cbind(total, common, lowfreq, rare)
      category <- colnames(categorymat_5loc)
      dimnames(resamp_category5loc) <- list(paste0("sample ", 1:nrow(categorymat_5loc)), category, paste0("replicate ", 1:25))
      # Store the results for the current replicate in resamp_category5loc
      resamp_category5loc[, , i] <- categorymat_5loc
    }
  }
  return(resamp_category5loc)
}
SNP_5kloc_function(QUAC.SNP.Wild.genind)

# 10k loci
SNP_10kloc_function <- function(test_genind){
  resamp_category5loc <- array(dim = c(nrow(test_genind@tab),4,25))
  for (i in 1:25) {
    samp_5loc <- sample(locNames(test_genind), size = 10000, replace = FALSE)
    MSAT.5loc.WILD.genind <- test_genind[, loc = samp_5loc]
    wildSamp_5loc <- MSAT.5loc.WILD.genind@tab
    wildSamp_5loc <- wildSamp_5loc[, which(colSums(wildSamp_5loc, na.rm = TRUE) != 0)]
    wildComplete <- colSums(wildSamp_5loc, na.rm = TRUE) / (nrow(wildSamp_5loc) * 2)
    wildSubset <- wildComplete[wildComplete > 0]
    total <- vector(length = nrow(wildSamp_5loc))
    common <- vector(length = nrow(wildSamp_5loc))
    lowfreq <- vector(length = nrow(wildSamp_5loc))
    rare <- vector(length = nrow(wildSamp_5loc))
    # Loop for each sample size
    for (j in 1:nrow(wildSamp_5loc)) {
      samp <- sample(nrow(wildSamp_5loc), size = j, replace = FALSE)
      samp <- wildSamp_5loc[samp,]
      # Measure the proportion of allelic representation in that sample
      if (j == 1) {
        total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(samp > 0)))]) / length(names(wildSubset))
        common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.1))
        lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
        rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(samp > 0)))) / length(which(wildSubset < 0.01))
      } else {
        total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))]) / length(names(wildSubset))
        common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.1))
        lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
        rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset < 0.01))
      }
      categorymat_5loc <- cbind(total, common, lowfreq, rare)
      category <- colnames(categorymat_5loc)
      dimnames(resamp_category5loc) <- list(paste0("sample ", 1:nrow(categorymat_5loc)), category, paste0("replicate ", 1:25))
      # Store the results for the current replicate in resamp_category5loc
      resamp_category5loc[, , i] <- categorymat_5loc
    }
  }
  return(resamp_category5loc)
}
SNP_10kloc_function(QUAC.SNP.Wild.genind)

# total loci
SNP_15kloc_function <- function(test_genind){
  resamp_category5loc <- array(dim = c(nrow(test_genind@tab),4,25))
  for (i in 1:25) {
    samp_5loc <- sample(locNames(test_genind), size = length(locNames(test_genind)), replace = FALSE)
    MSAT.5loc.WILD.genind <- test_genind[, loc = samp_5loc]
    wildSamp_5loc <- MSAT.5loc.WILD.genind@tab
    wildSamp_5loc <- wildSamp_5loc[, which(colSums(wildSamp_5loc, na.rm = TRUE) != 0)]
    wildComplete <- colSums(wildSamp_5loc, na.rm = TRUE) / (nrow(wildSamp_5loc) * 2)
    wildSubset <- wildComplete[wildComplete > 0]
    total <- vector(length = nrow(wildSamp_5loc))
    common <- vector(length = nrow(wildSamp_5loc))
    lowfreq <- vector(length = nrow(wildSamp_5loc))
    rare <- vector(length = nrow(wildSamp_5loc))
    # Loop for each sample size
    for (j in 1:nrow(wildSamp_5loc)) {
      samp <- sample(nrow(wildSamp_5loc), size = j, replace = FALSE)
      samp <- wildSamp_5loc[samp,]
      # Measure the proportion of allelic representation in that sample
      if (j == 1) {
        total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(samp > 0)))]) / length(names(wildSubset))
        common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.1))
        lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
        rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(samp > 0)))) / length(which(wildSubset < 0.01))
      } else {
        total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))]) / length(names(wildSubset))
        common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.1))
        lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
        rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset < 0.01))
      }
      categorymat_5loc <- cbind(total, common, lowfreq, rare)
      category <- colnames(categorymat_5loc)
      dimnames(resamp_category5loc) <- list(paste0("sample ", 1:nrow(categorymat_5loc)), category, paste0("replicate ", 1:25))
      # Store the results for the current replicate in resamp_category5loc
      resamp_category5loc[, , i] <- categorymat_5loc
    }
  }
  return(resamp_category5loc)
}
SNP_15kloc_function(QUAC.SNP.Wild.genind)

## SNP QUBO ##
SNP_2kloc_function <- function(test_genind){
  resamp_category5loc <- array(dim = c(nrow(test_genind@tab),4,25))
  for (i in 1:25) {
    samp_5loc <- sample(locNames(test_genind), size = 5000, replace = FALSE)
    MSAT.5loc.WILD.genind <- test_genind[, loc = samp_5loc]
    wildSamp_5loc <- MSAT.5loc.WILD.genind@tab
    wildSamp_5loc <- wildSamp_5loc[, which(colSums(wildSamp_5loc, na.rm = TRUE) != 0)]
    wildComplete <- colSums(wildSamp_5loc, na.rm = TRUE) / (nrow(wildSamp_5loc) * 2)
    wildSubset <- wildComplete[wildComplete > 0]
    total <- vector(length = nrow(wildSamp_5loc))
    common <- vector(length = nrow(wildSamp_5loc))
    lowfreq <- vector(length = nrow(wildSamp_5loc))
    rare <- vector(length = nrow(wildSamp_5loc))
    # Loop for each sample size
    for (j in 1:nrow(wildSamp_5loc)) {
      samp <- sample(nrow(wildSamp_5loc), size = j, replace = FALSE)
      samp <- wildSamp_5loc[samp,]
      # Measure the proportion of allelic representation in that sample
      if (j == 1) {
        total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(samp > 0)))]) / length(names(wildSubset))
        common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.1))
        lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
        rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(samp > 0)))) / length(which(wildSubset < 0.01))
      } else {
        total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))]) / length(names(wildSubset))
        common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.1))
        lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
        rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset < 0.01))
      }
      categorymat_5loc <- cbind(total, common, lowfreq, rare)
      category <- colnames(categorymat_5loc)
      dimnames(resamp_category5loc) <- list(paste0("sample ", 1:nrow(categorymat_5loc)), category, paste0("replicate ", 1:25))
      # Store the results for the current replicate in resamp_category5loc
      resamp_category5loc[, , i] <- categorymat_5loc
    }
  }
  return(resamp_category5loc)
}
SNP_2kloc_function(QUBO.SNP.Wild.genind)

SNP_4kloc_function <- function(test_genind){
  resamp_category5loc <- array(dim = c(nrow(test_genind@tab),4,25))
  for (i in 1:25) {
    samp_5loc <- sample(locNames(test_genind), size = 4000, replace = FALSE)
    MSAT.5loc.WILD.genind <- test_genind[, loc = samp_5loc]
    wildSamp_5loc <- MSAT.5loc.WILD.genind@tab
    wildSamp_5loc <- wildSamp_5loc[, which(colSums(wildSamp_5loc, na.rm = TRUE) != 0)]
    wildComplete <- colSums(wildSamp_5loc, na.rm = TRUE) / (nrow(wildSamp_5loc) * 2)
    wildSubset <- wildComplete[wildComplete > 0]
    total <- vector(length = nrow(wildSamp_5loc))
    common <- vector(length = nrow(wildSamp_5loc))
    lowfreq <- vector(length = nrow(wildSamp_5loc))
    rare <- vector(length = nrow(wildSamp_5loc))
    # Loop for each sample size
    for (j in 1:nrow(wildSamp_5loc)) {
      samp <- sample(nrow(wildSamp_5loc), size = j, replace = FALSE)
      samp <- wildSamp_5loc[samp,]
      # Measure the proportion of allelic representation in that sample
      if (j == 1) {
        total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(samp > 0)))]) / length(names(wildSubset))
        common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.1))
        lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
        rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(samp > 0)))) / length(which(wildSubset < 0.01))
      } else {
        total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))]) / length(names(wildSubset))
        common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.1))
        lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
        rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset < 0.01))
      }
      categorymat_5loc <- cbind(total, common, lowfreq, rare)
      category <- colnames(categorymat_5loc)
      dimnames(resamp_category5loc) <- list(paste0("sample ", 1:nrow(categorymat_5loc)), category, paste0("replicate ", 1:25))
      # Store the results for the current replicate in resamp_category5loc
      resamp_category5loc[, , i] <- categorymat_5loc
    }
  }
  return(resamp_category5loc)
}
SNP_4kloc_function(QUBO.SNP.Wild.genind)

SNP_QUBO_totalloc_function <- function(test_genind){
  resamp_category5loc <- array(dim = c(nrow(test_genind@tab),4,25))
  for (i in 1:25) {
    samp_5loc <- sample(locNames(test_genind), size = length(locNames(test_genind)), replace = FALSE)
    MSAT.5loc.WILD.genind <- test_genind[, loc = samp_5loc]
    wildSamp_5loc <- MSAT.5loc.WILD.genind@tab
    wildSamp_5loc <- wildSamp_5loc[, which(colSums(wildSamp_5loc, na.rm = TRUE) != 0)]
    wildComplete <- colSums(wildSamp_5loc, na.rm = TRUE) / (nrow(wildSamp_5loc) * 2)
    wildSubset <- wildComplete[wildComplete > 0]
    total <- vector(length = nrow(wildSamp_5loc))
    common <- vector(length = nrow(wildSamp_5loc))
    lowfreq <- vector(length = nrow(wildSamp_5loc))
    rare <- vector(length = nrow(wildSamp_5loc))
    # Loop for each sample size
    for (j in 1:nrow(wildSamp_5loc)) {
      samp <- sample(nrow(wildSamp_5loc), size = j, replace = FALSE)
      samp <- wildSamp_5loc[samp,]
      # Measure the proportion of allelic representation in that sample
      if (j == 1) {
        total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(samp > 0)))]) / length(names(wildSubset))
        common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.1))
        lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(samp > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
        rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(samp > 0)))) / length(which(wildSubset < 0.01))
      } else {
        total[j] <- length(names(wildSubset)[which(names(wildSubset) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))]) / length(names(wildSubset))
        common[j] <- length(which(names(which(wildSubset > 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.1))
        lowfreq[j] <- length(which(names(which(wildSubset > 0.01 & wildSubset < 0.1)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset > 0.01 & wildSubset < 0.1))
        rare[j] <- length(which(names(which(wildSubset < 0.01)) %in% names(which(colSums(samp, na.rm = TRUE) > 0)))) / length(which(wildSubset < 0.01))
      }
      categorymat_5loc <- cbind(total, common, lowfreq, rare)
      category <- colnames(categorymat_5loc)
      dimnames(resamp_category5loc) <- list(paste0("sample ", 1:nrow(categorymat_5loc)), category, paste0("replicate ", 1:25))
      # Store the results for the current replicate in resamp_category5loc
      resamp_category5loc[, , i] <- categorymat_5loc
    }
  }
  return(resamp_category5loc)
}
SNP_QUBO_totalloc_function(QUBO.SNP.Wild.genind)


QUAC_MSAT_resamp5loc <- MSAT_5loc_function(QUAC.MSAT.Wild.genind)
QUBO_MSAT_resamp3loc <- MSAT_5loc_function(QUBO.MSAT.Wild.genind)
QUAC_MSAT_resamp10loc <- MSAT_10loc_function(QUAC.MSAT.Wild.genind)
QUBO_MSAT_resamp6loc <- MSAT_10loc_function(QUBO.MSAT.Wild.genind)
QUAC_MSAT_resamptotloc <- MSAT_totalloc_function(QUAC.MSAT.Wild.genind)
QUBO_MSAT_resamptotloc <- MSAT_totalloc_function(QUBO.MSAT.Wild.genind)
QUAC_SNP_resamp5kloc <- SNP_5kloc_function(QUAC.SNP.Wild.genind)
QUAC_SNP_resamp10kloc <- SNP_10kloc_function(QUAC.SNP.Wild.genind)
QUAC_SNP_resamptotloc <- SNP_15kloc_function(QUAC.SNP.Wild.genind)
QUBO_SNP_resamp2kloc <- SNP_2kloc_function(QUBO.SNP.Wild.genind)
QUBO_SNP_resamp4kloc <- SNP_4kloc_function(QUBO.SNP.Wild.genind)
QUBO_SNP_resamptotloc <- SNP_totalloc_function(QUBO.SNP.Wild.genind)


## Linear model ##
# pass all arrays to a dataframe using the function. define the function that takes an input 'data_array'
analyze_resampling_array <- function(data_array) {
  # linear model of resampling array. Extract the column named 'total' from the array.
  # concatenate the extracted column values into the vector.
  totalsVector <- c(data_array[,"total",]) 
  
  # Specify sample numbers column.
  # Create a vector of sample numbers from 1 to the number of rows in the 'total' column 
  gm_sampleNumbers <- 1:(nrow(data_array[,"total",]))
  # Repeat the sample numbers vector for the number of replicates
  gm_sampleNumbers <- rep(gm_sampleNumbers, dim(data_array)[[3]])
  
  # Create data frame from resampling array values
  gm_DF <- data.frame(sampleNumbers=gm_sampleNumbers, totalValues=totalsVector)
  
  # Build and analyze linear models
  gm_Model <- lm(sampleNumbers ~ I((totalValues)^3), data = gm_DF)
  # Create a new data fram 'gm_newData with a single column 'totalValues' containing the value 0.95
  gm_newData <- data.frame(totalValues=0.95)
  # Use the linear model to predict the response for the new data frame. Specify 'interval = prediction to obtain a prediction interval
  gm_95MSSEprediction <- predict(gm_Model, gm_newData, interval = "prediction")
  
  # Pass the gm_95MSSEprediction to the object storing our results 
  # Store the predicted values and predictino interval in the object named 'result' 
  result <- gm_95MSSEprediction
  # Calculate the width of the prediction interval by substracting the lower limit from the upper limit
  piWidth <- gm_95MSSEprediction[3] - gm_95MSSEprediction[2]
  # Return a list containing the predicted values and the width of the prediction interval
  return(list(result = result, piWidth = piWidth))
}

# this is a list of the resampling arrays which will be iterated in a loop to 
# execute the analyze_resampling_array function which will output the prediction interval
# values and prediction interval widths. 
array_list <- list(
  QUAC_MSAT_resamp5loc,
  QUAC_MSAT_resamp10loc,
  QUAC_MSAT_resamptotloc,
  QUBO_MSAT_resamp3loc,
  QUBO_MSAT_resamp6loc,
  QUBO_MSAT_resamptotloc,
  QUAC_SNP_resamp5kloc,
  QUAC_SNP_resamp10kloc,
  QUAC_SNP_resamptotloc,
  QUBO_SNP_resamp2kloc,
  QUBO_SNP_resamp4kloc,
  QUBO_SNP_resamptotloc
)
# this matrix will store the pi values and pi widths
# Create an empty matrix to store the results
results_matrix <- matrix(nrow = length(array_list), ncol = 4)
# Set column names for 'results_Matrix'
colnames(results_matrix) <- c("fit", "lower", "upper", "piWidth")
# Set row names for 'results_matrix'
rownames(results_matrix) <- c(
  "QUAC_MSAT_resamp5loc",
  "QUAC_MSAT_resamp10loc",
  "QUAC_MSAT_resamptotloc",
  "QUBO_MSAT_resamp3loc",
  "QUBO_MSAT_resamp6loc",
  "QUBO_MSAT_resamptotloc",
  "QUAC_SNP_resamp5kloc",
  "QUAC_SNP_resamp10kloc",
  "QUAC_SNP_resamptotloc",
  "QUBO_SNP_resamp2kloc",
  "QUBO_SNP_resamp4kloc",
  "QUBO_SNP_resamptotloc"
)

# Iterate through the arrays and store results in the matrix
# Initiate loop that iterates over the idicies of 'array_list'
for (i in 1:length(array_list)) {
  # Call the 'analyze_resampling_array' function on the ith element of 'array_list'. This function returns a list with 
  # resul and piWidth values.
  analyze_resampling_array(array_list[[i]])
  # Store results and piWidth values in the ith row of the matrix
  results_matrix[i, ] <- c(analyze_resampling_array(array_list[[i]])$result, 
                           analyze_resampling_array(array_list[[i]])$piWidth)
  
}
# Print final 'results_matrix' showing the fit, lower, upper and piWidth values for each resampling array
print(results_matrix)

write.csv(results_matrix, 
          file = "C:/Users/gsalas/Documents/resampling_CIs/Code/Outputs/QUAC_QUBO_resamp_loci.csv",
          row.names = TRUE)