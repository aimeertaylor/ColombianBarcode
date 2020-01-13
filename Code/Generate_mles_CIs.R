###################################################################
# This script is adapted from PlasmodiumRelatedness/Generate_mles.R
# s.t. all mles and parametric bootstrap mles for Colombia only
# are calculated in one script.
# 134690 seconds on my Pro with nrepeat = 500. 
# consider changing nboot to 100 and dropping those for k
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

funif = runif(250)
frequencies_unif <- cbind(funif, 1-funif) 
frequencies_true <- cbind(1-data_set$fs, data_set$fs)
f_list = list(frequencies_true = frequencies_true, frequencies_unif = frequencies_unif)

for(f in 2){
  
  frequencies = f_list[[f]]
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
  save(mle_CIs, file = sprintf("../RData/mles_%s.RData", names(f_list)[f]))
}

