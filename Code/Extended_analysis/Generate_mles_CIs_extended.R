###############################################################################
#' This script is adapted from Generate_mles_CIs.R. Plan: run using frequencies
#' from 1) Taylor et al 2020 and check mles match; 2) run using frequencies
#' computed using the extended data. As of Dec 10th, waiting for complete data
#' from Angela (61 extra samplesand remove low quality) and updated metadata
###############################################################################
rm(list = ls())
set.seed(1)
library(ggplot2)
library(dplyr)
library(Rcpp)
library(doParallel)
library(doRNG)
source("~/Dropbox/IBD_IBS/PlasmodiumRelatedness/Code/simulate_data.R") # Download this script from https://github.com/artaylor85/PlasmodiumRelatedness
sourceCpp("~/Dropbox/IBD_IBS/PlasmodiumRelatedness/Code/hmmloglikelihood.cpp") # Download this script from https://github.com/artaylor85/PlasmodiumRelatedness
registerDoParallel(cores = detectCores()-1)
epsilon <- 0.001 # Fix epsilon throughout
nboot <- 10 # For CIs (with nboot = 5 using 3 cores = 3.436389 hours; 41.23667 mins per bootstrap; 69 mins for 100)
set.seed(1) # For reproducibility
Ps = c(0.025, 0.975) # CI quantiles

# Load the extended data set
load(file = "../../RData/snpdata_extended.RData")

freqs_to_use <- "Taylor2020" # "AllAvailable"

## Mechanism to compute MLE given fs, distances, Ys, epsilon
compute_rhat_hmm <- function(frequencies, distances, Ys, epsilon){
  ndata <- nrow(frequencies)
  ll <- function(k, r) loglikelihood_cpp(k, r, Ys, frequencies, distances, epsilon, rho = 7.4 * 10^(-7))
  optimization <- optim(par = c(50, 0.5), fn = function(x) - ll(x[1], x[2]))
  rhat <- optimization$par
  return(rhat)
}

## Mechanism to generate Ys given fs, distances, k, r, epsilon
simulate_Ys_hmm <- function(frequencies, distances, k, r, epsilon){
  Ys <- simulate_data(frequencies, distances, k = k, r = r, epsilon, rho = 7.4 * 10^(-7))
  return(Ys)
}

#=====================================
# Create indices for pairwise comparisons
#=====================================
individual_names <- names(snpdata)[-(1:2)]
nindividuals <- length(individual_names)
name_combinations <- matrix(nrow = nindividuals*(nindividuals-1)/2, ncol = 2)
count <- 0
for (i in 1:(nindividuals-1)){
  for (j in (i+1):nindividuals){
    count <- count + 1
    name_combinations[count,1] <- individual_names[i]
    name_combinations[count,2] <- individual_names[j]
  }
}

#=====================================
# Sort data by chromosome and position and extract frequencies
#=====================================
snpdata <- snpdata %>% arrange(chrom, pos) 
snpdata$dt <- c(diff(snpdata$pos), Inf)
pos_change_chrom <- 1 + which(diff(snpdata$chrom) != 0) # find places where chromosome changes
snpdata$dt[pos_change_chrom-1] <- Inf


#========================================================
if(freqs_to_use == "Taylor2020"){
  
  # Cannot use frequencies from TxtData/hmmInput.txt since it is encoded differently
  # s.t. some zeros in new encoding will be ones in old and vice versa without
  # knowing which
  load(file = "../../RData/SNPData.RData")
  sids <- as.character(SNPData$SAMPLE.CODE)
  if (!all(sids %in% colnames(snpdata))) stop("Some samples from Taylor2020 are missing")
  snpdata$fs = rowMeans(snpdata[,-(1:2)][,sids], na.rm = TRUE) # Calculate frequencies
  frequencies = cbind(1-snpdata$fs, snpdata$fs)

} else if (freqs_to_use == "AllAvailable") {
  
  # Calculate frequencies using all available data
  snpdata$fs = rowMeans(snpdata[,-(1:2)], na.rm = TRUE) 
  frequencies = cbind(1-snpdata$fs, snpdata$fs)
  
}
#========================================================



# Check ordering of markers
plot(snpdata$chrom, type = 'l')
plot(snpdata$pos, type = 'l')

# Check all frequencies in (0,1): 
all(snpdata$fs > 0 & snpdata$fs < 1)
nrow(snpdata) # One marker failed throughout


#=====================================
# Calculate mles 
#=====================================
system.time(
  
  # For each pair...  
  mle_CIs <- foreach(icombination = 1:nrow(name_combinations),.combine = rbind) %dorng% {
    
    # Let's focus on one pair of individuals
    individual1 <- name_combinations[icombination,1]
    individual2 <- name_combinations[icombination,2]
    
    # Indeces of pair
    i1 <- which(individual1 == names(snpdata))
    i2 <- which(individual2 == names(snpdata))
    
    # Extract data 
    subdata <- cbind(snpdata[,c("fs","dt")],snpdata[,c(i1,i2)])
    names(subdata) <- c("fs","dt","Yi","Yj")
    
    # Generate mle
    krhat_hmm <- compute_rhat_hmm(frequencies, distances = subdata$dt, 
                                  Ys = cbind(subdata$Yi, subdata$Yj), epsilon)
    
    # Generate parametric bootstrap mles 
    krhats_hmm_boot = foreach(iboot = 1:nboot, .combine = rbind) %dorng% {
      Ys_boot <- simulate_Ys_hmm(frequencies, distances = snpdata$dt, k = krhat_hmm[1], r = krhat_hmm[2], epsilon)
      compute_rhat_hmm(frequencies, subdata$dt, Ys_boot, epsilon)
    }
    
    CIs = apply(krhats_hmm_boot, 2, function(x)quantile(x, probs=Ps)) # Few seconds
    X = data.frame('individual1' = individual1, 'individual2' = individual2, 
                   rhat = krhat_hmm[2], 'r2.5%' = CIs[1,2], 'r97.5%' = CIs[2,2],
                   khat = krhat_hmm[1], 'k2.5%' = CIs[1,1], 'k97.5%' = CIs[2,1])
    X
  })

save(mle_CIs, file = sprintf("../../RData/mles_CIs_extended_freqs%s.RData", freqs_to_use))


