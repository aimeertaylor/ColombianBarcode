###################################################################
# This script is adapted from PlasmodiumRelatedness/Generate_mles.R
# s.t. all mles and parametric bootstrap mles for Colombia only
# are calculated in one script.
# 134690 seconds on my Pro with nrepeat = 500. 
###################################################################
rm(list = ls())
library(ggplot2)
library(dplyr)
library(Rcpp)
library(doParallel)
library(doRNG)
source("~/Dropbox/IBD_IBS/PlasmodiumRelatedness/Code/simulate_data.R") # Download this script from https://github.com/artaylor85/PlasmodiumRelatedness
sourceCpp("~/Dropbox/IBD_IBS/PlasmodiumRelatedness/Code/hmmloglikelihood.cpp") # Download this script from https://github.com/artaylor85/PlasmodiumRelatedness
registerDoParallel(cores = detectCores()-1)
epsilon <- 0.001 # Fix epsilon throughout
nboot <- 100 # For CIs 100 with 3 cores 12hr
Ps = c(0.025, 0.975) # CI quantiles

## Mechanism to compute MLE given fs, distances, Ys, epsilon
compute_rhat_hmm <- function(frequencies, distances, Ys, epsilon){
  if(any(is.na(Ys))) stop ("NAs in Ys not allowed")
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

## Function to compute distances (stand-apart function because need to 
# recompute distances everytime there are NAs in Ys)
compute_distances <- function(x){
  # compute distances
  x$dt <- c(diff(x$pos), Inf)
  # find places where chromosome changes
  pos_change_chrom <- 1 + which(diff(x$chrom) != 0) 
  # replace chromosome limits with Inf
  x$dt[pos_change_chrom-1] <- Inf
  # end of function
  return(x$dt)
}

# Load data
data_set = read.delim("../TxtData/hmmInput.txt") 

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
frequencies <- cbind(1-data_set$fs, data_set$fs)
pos_chrom <- data_set[,c("pos", "chrom")]

system.time(
  
  # For each pair...  
  mle_CIs <- foreach(icombination = 1:nrow(name_combinations),.combine = rbind) %dorng% {
    
    # Let's focus on one pair of individuals
    individual1 <- name_combinations[icombination,1]
    individual2 <- name_combinations[icombination,2]
    
    # Indeces of pair
    i1 <- which(individual1 == names(data_set))
    i2 <- which(individual2 == names(data_set))
    
    # Extract data and check if informative
    subdata <- cbind(Yi = data_set[,i1], Yj = data_set[,i2])
    snps_to_keep_ind <- apply(subdata, 1, function(x) !any(is.na(x))) 
    snp_count <- sum(snps_to_keep_ind)
    
    if (snp_count == 0) next()
  
    # Extract non-NA observed genotypes
    Ys_ <- subdata[snps_to_keep_ind,,drop = F]
    frequencies_ <- frequencies[snps_to_keep_ind,,drop=F]
    distances_ <- compute_distances(pos_chrom[snps_to_keep_ind,])
    
    # Generate mle
    krhat_hmm <- compute_rhat_hmm(frequencies = frequencies_, 
                                  distances = distances_, 
                                  Ys = Ys_, 
                                  epsilon)
  
    # Generate parametric bootstrap mles 
    set.seed(1) # For reproducibility 
    krhats_hmm_boot = foreach(iboot = 1:nboot, .combine = rbind) %dorng% {
      Ys_boot <- simulate_Ys_hmm(frequencies = frequencies_, distances = distances_, 
                                 k = krhat_hmm[1], r = krhat_hmm[2], epsilon)
      compute_rhat_hmm(frequencies_, distances_, Ys_boot, epsilon)
    }
    
  
    CIs = apply(krhats_hmm_boot, 2, function(x)quantile(x, probs=Ps)) # Few seconds
    X = data.frame('individual1' = individual1, 'individual2' = individual2, 
                   rhat = krhat_hmm[2], 'r2.5%' = CIs[1,2], 'r97.5%' = CIs[2,2],
                   khat = krhat_hmm[1], 'k2.5%' = CIs[1,1], 'k97.5%' = CIs[2,1])
    X
  })

save(mle_CIs, file = "../RData/mles_CIs.RData")
     
     
     