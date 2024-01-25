# Overview
This repository contains code that builds upon Quercus_IUCN_samp_sims, a previous simulation project by Kaylee Rosenberger. The goal of this subproject is to assess the variation on guidelines for minimum sample sizes required to maintain genetic diversity of various ex situ oak collections. In this subproject, we calculated 95% confidence intervals around minimum sample sizes required for 95% allelic representation across 1000 replicates for 14 IUCN red list endangered oak species using genetic marker resampling arrays. By calculating the width of the confidence intervals, the confidence of the guidelines can be quantitatively determined. Additionally, we used linear regression, a statistically rigorous technique to predict the 95% allelic representation as a function of minimum sampling size for Q. acerifolia and Q. boyntonii. We also explore how the 95% MSSE changes based on sampling of few to many loci for alleles of different frequency categories using MSAT and SNP genetic marker datasets.

### Datasets
**LociBootstrapping_datasets**: Contains MSAT and SNP genind objects for Q. acerifolia and Q. boyntonii saved as R objects. The slot of the genind object for each genetic marker of each species containing allele counts were filtered to only include samples in common with each other.  
**Subset**:   
The file in this folder contains the dataset used to calculate the confidence intervals and confidence interval widths for each of the 14 species.   
QUAC.MSAT.Complete_resampArr.Rdata  
QUAC.SNP.DN.R0.Complete_resampArr.Rdata  
quercus_final_results_orig.Rdata
**Outputs**

