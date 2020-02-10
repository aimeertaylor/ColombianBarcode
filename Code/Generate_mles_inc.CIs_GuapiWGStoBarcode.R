###################################################################
# This script is adapted from Generate_mles_CIs.R
# At present, using frequencies calucted from the combined data sets
# As such, results for Diego's data set will differ slightly
###################################################################
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
nboot <- 100 # For CIs 
set.seed(1) # For reproducibility
Ps = c(0.025, 0.975) # CI quantiles

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

# Load and process data =============
data_set0 = read.delim("../TxtData/hmmInput.txt") # Original data set from Echeverry et al. 
data_set1 = read.csv("../OriginalData/Guapi_WGStoBarcode.csv") # Data set provided by Vladimir courtesy of Manuela 

# Match colnames of data_set1 to data_set0 
colnames(data_set1) = gsub("chr", "chrom", colnames(data_set1)) # Reformat colnames

# Match rownames of data_set1 to data_set0 
rownames(data_set1) = paste(data_set1$chrom, data_set1$pos, sep = "_")
rownames(data_set0) = paste(data_set0$chrom, data_set0$pos, sep = "_")
missing_SNPs = which(!rownames(data_set0) %in% rownames(data_set1))

# Add rows of missing SNPs to data_set1
X = cbind(data_set0$chrom[missing_SNPs], data_set1$pos[missing_SNPs],  
      matrix(NA, nrow = length(missing_SNPs), ncol = ncol(data_set1)-2))
rownames(X) = rownames(data_set0)[missing_SNPs]
colnames(X) = colnames(data_set1)
data_set1 = rbind(data_set1, X)

# Concatenate data and remove rownames
data_set = cbind(data_set0, data_set1[rownames(data_set0),-(1:2)])
rownames(data_set) = NULL
#=====================================

# Create indices for pairwise comparisons
individual_names <- names(data_set)[-(1:2)]
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

# Sort data by chromosome and position and add frequencies
data_set <- data_set %>% arrange(chrom, pos) 
data_set$fs = rowMeans(data_set[,-(1:2)], na.rm = TRUE) # Calculate frequencies
data_set$dt <- c(diff(data_set$pos), Inf)
pos_change_chrom <- 1 + which(diff(data_set$chrom) != 0) # find places where chromosome changes
data_set$dt[pos_change_chrom-1] <- Inf


frequencies = cbind(1-data_set$fs, data_set$fs)
system.time(
  
  # For each pair...  
  mle_CIs <- foreach(icombination = 1:nrow(name_combinations),.combine = rbind) %dorng% {
    
    # Let's focus on one pair of individuals
    individual1 <- name_combinations[icombination,1]
    individual2 <- name_combinations[icombination,2]
    
    # Indeces of pair
    i1 <- which(individual1 == names(data_set))
    i2 <- which(individual2 == names(data_set))
    
    # Extract data 
    subdata <- cbind(data_set[,c("fs","dt")],data_set[,c(i1,i2)])
    names(subdata) <- c("fs","dt","Yi","Yj")
    
    # Generate mle
    krhat_hmm <- compute_rhat_hmm(frequencies, subdata$dt, cbind(subdata$Yi, subdata$Yj), epsilon)
    
    # Generate parametric bootstrap mles 
    krhats_hmm_boot = foreach(iboot = 1:nboot, .combine = rbind) %dorng% {
      Ys_boot <- simulate_Ys_hmm(frequencies, distances = data_set$dt, k = krhat_hmm[1], r = krhat_hmm[2], epsilon)
      compute_rhat_hmm(frequencies, subdata$dt, Ys_boot, epsilon)
    }
    
    CIs = apply(krhats_hmm_boot, 2, function(x)quantile(x, probs=Ps)) # Few seconds
    X = data.frame('individual1' = individual1, 'individual2' = individual2, 
                   rhat = krhat_hmm[2], 'r2.5%' = CIs[1,2], 'r97.5%' = CIs[2,2],
                   khat = krhat_hmm[1], 'k2.5%' = CIs[1,1], 'k97.5%' = CIs[2,1])
    X
  })

save(mle_CIs, file = "../RData/mles_inc.GuapiWGStoBarcode.RData")


