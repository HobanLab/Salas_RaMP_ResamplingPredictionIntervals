# Overview
This repository contains code that builds upon Quercus_IUCN_samp_sims, a previous simulation project by Kaylee Rosenberger. The goal of this subproject is to assess the variation on guidelines for minimum sample sizes required to maintain genetic diversity of various ex situ oak collections. In this subproject, we calculated 95% confidence intervals around minimum sample sizes required for 95% allelic representation across 1000 replicates for 14 IUCN red list endangered oak species using genetic marker resampling arrays. By calculating the width of the confidence intervals, the confidence of the guidelines can be quantitatively determined. Additionally, we used linear regression, a statistically rigorous technique to predict the 95% allelic representation as a function of minimum sampling size for Q. acerifolia and Q. boyntonii. We also explore how the 95% MSSE changes based on sampling of few to many loci for alleles of different frequency categories using MSAT and SNP genetic marker datasets.
## Directory
### Datasets
QUAC.MSAT.Complete_resampArr.Rdata (a microsatellite resampling array containing allele frequency category data that does not filter loci shared between garden and wild populations of Q. acerifolia)  
QUAC.SNP.DN.R0.Complete (a singlenucleotide polymorphism resampling array containing allele frequency category data that used the de novo processing approach, does not filter the number of samples that share 0% of loci, and does not filter loci shared between garden and wild populations of Q. acerifolia.)
QUAC.SNP.DN.R80.Complete_resampArr.Rdata (a single nucleotide polymorphism resampling array containing allele frequency category data that used the de novo processing approach, filters the number of samples that share 80% of loci, and does not filter loci shared between garden and wild populations of Q. acerifolia.)
quercus_final_results_orig.Rdata  
**LociBootstrapping_datasets**: Contains MSAT and SNP genind objects for wild populations of Q. acerifolia and Q. boyntonii saved as R objects. The slot of the genind object for each genetic marker of each species containing allele counts were filtered to only include samples with common alleles.  
**Subset**: Contains MSAT and SNP resampling arrays that filter loci shared between garden and wild populations of Q. acerifolia. Both SNP resampling arrays use the reference alignment processing approach, while only one SNP resampling array filters the number of samples shared by 80% of loci.
### Outputs
14CIWidths.csv  
datasetCIwidths.csv  
gm_predict.csv  
QUAC_resamp_loci.csv  
QUBO_resamp_loci.csv  


**QUAC_MSSE_Quantiles.R**: script calculates means, quantiles, and generates plots for the total and other allele categories in order to create confidence intervals around 95% minimum sample size estimates.   
**QUAC_QUBO_loci_bootstrapping.R**: script builds resampling arrays that randomly sample loci at different ranges, calculates prediction intervals aroud the 95% minimum sample size estimate using MSAT and SNP genind objects for QUAC and QUBO, and builds a matrix that stores the PI values     
**MSSE_PredictionIntervals.R**: scripts calculates prediction interval values aroud the 95% minimum sample size estimate using MSAT and SNP resampling arrays for 14 oak species and builds a matrix that stores the PI values
