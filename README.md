# Overview
This repository contains code that builds upon [Quercus_IUCN_samp_sims](https://github.com/HobanLab/Quercus_IUCN_samp_sims), a previous simulation project by Kaylee Rosenberger. The goal of this subproject is to assess the variation in minimum sample size estimates (MSSEs) required to maintain genetic diversity of various ex situ oak collections. In this subproject, we calculate prediction intervals around MSSEs required for 95% allelic representation. We also explore how the 95% MSSE changes based on sampling of few to many loci for alleles of different frequency categories using MSAT and SNP genetic marker datasets.

## Repository structure
### Scripts
1. `QUAC_MSSE_Quantiles.R`
  - Script calculates MSSE means and quantiles, and generates plots for the total allelic representation (and other categories of allelic frequency) in order to create confidence intervals around 95% minimum sample size estimates. The approach used in this script for calculating allelic representation confidence intervals is improved upon by using the `predict` function (see `MSSE_PredictionIntervals.R`).
    - **Inputs** 
      - `QUAC_Subset_resampArrs` folder--resampling arrays built from _Quercus acerifolia_ (QUAC) microsatellite (MSAT) and single nucleotide polymorphism (SNP) genetic data (for SNPs, R0 and R80). These datasets are all subset to the same number of samples, to allow for greater comparability between marker types and missing data levels.
      
2. `MSSE_PredictionIntervals.R`
  - Script calculates the prediction interval (PI) values around the 95% MSSE using two different datasets (QUAC and IUCN 14 oaks), and builds a matrix that stores the PI values
    - **Inputs**
      - _QUAC_: `QUAC_Subset_resampArrs` folder--resampling arrays built from _Quercus acerifolia_ (QUAC) microsatellite (MSAT) and single nucleotide polymorphism (SNP) genetic data (for SNPs, R0 and R80). These datasets are all subset to the same number of samples, to allow for greater comparability between marker types and missing data levels.
      - _IUCN 14 oaks_: `quercus_final_results_orig.Rdata`--a resampling array containing the total allelic representation values for oaks simulated by Kaylee Rosenberger; source code found in the [`Quercus_IUCN_samp_sims` repo](https://github.com/HobanLab/Quercus_IUCN_samp_sims)
    - **Outputs**
      - _QUAC_: `QUAC_PI_values.csv`
      - _IUCN 14 oaks_: `Quercus14_PI_values.csv`
      
3. `QUAC_QUBO_loci_bootstrapping.R` 
  - Script builds resampling arrays based on different ranges of randomly sampled loci, calculates the prediction intervals around the 95% MSSEs, and builds a matrix that stores the PI values
    - **Inputs**
      - `LociBootstrapping_Datasets` folder--`genpop` objects for wild populations of _Q. acerifolia_ (QUAC) and _Q. boyntonii_ (QUBO), saved as R objects. 
    - **Outputs**
      - `QUAC_MSSE_Quantiles.csv`

### Datasets
This folder contains the input files read in by the analyses in the `Scripts` folder (see outline of **Inputs** above). These files are typically either resampling arrays (sets of allelic representation values, for a given number of randomly drawn samples) or `genpop` objects (read in using the `adegenet` library) from which resampling arrays are built.

### Outputs
This folder contains the CSV outputs generated by the analyses in the `Scripts` folder (see outline of **Outputs** above). Generally, the contents of these CSVs are minimum sample size estimates and upper/lower the confidence intervals (CI) or prediction intervals (CI) bounding them.
