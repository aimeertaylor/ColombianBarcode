##############################################################################################
# Script to generate fst estimates using Colombian data
# Adapted from script used for TM data
##############################################################################################

rm(list = ls())
require('foreach') # For running in parallel 
require('doMC') # For running in parallel 
require('rngtools') # For running in parallel 
require('abind') 
require('plyr') # For aplyr

# Source data
source('/Users/aimeet/Documents/BroadLaptop/TM_border/QuantLocalPfConnIBD/FunctionFiles/calculate_pairwise_Fst.R')
load('/Users/aimeet/Documents/BroadLaptop/ColombianBarcode/RData/SNPData.RData')
load('/Users/aimeet/Documents/BroadLaptop/ColombianBarcode/RData/geo_dist_info.RData')

attach(geo_dist_info, warn.conflicts = FALSE)
Pair_wise_site_comparisons <- pairwise_site_distance[,c(1,2)]
numComparisons <- nrow(Pair_wise_site_comparisons)

# Magic no.s/variables
registerDoMC(cores = 3) # For foreach
n_repeat <- 1000

# ============================================================================
# Function to calculate pairwise Fst estimates 
# =============================================================================
fst_calculations <- function(SNPDataBinary){
  
  numSNPs <- ncol(SNPDataBinary)
  numSamples <- nrow(SNPDataBinary)
  
  #----------------------------------------------------------------------------
  # Between citites averaging over years 
  #----------------------------------------------------------------------------
  Pair_wise_site_comparisons_Fst <- array(dim = numComparisons,dimnames = list(geo_order))
  Pair_wise_site_comparisons_Fst_perSNP <- array(dim=c(numSNPs,numComparisons),dimnames = list(NULL,geo_order))

  for(i in 1:numComparisons){
    
    # Indices for populations to compare
    P1 <- (SNPData$City == as.character(Pair_wise_site_comparisons[i,1]))
    P2 <- (SNPData$City == as.character(Pair_wise_site_comparisons[i,2]))
    P1_name <- as.character(Pair_wise_site_comparisons[i,1])
    P2_name <- as.character(Pair_wise_site_comparisons[i,2])
    
    # Calculate Fst
    X_Reich <- calculate_pairwise_reich(P1, P2, SNPDataBinary)
    
    # Save to matrix
    Pair_wise_site_comparisons_Fst[i] <- X_Reich$F_st
    Pair_wise_site_comparisons_Fst_perSNP[,i] <- X_Reich$F_st_snp
  }
  
  Fst_barcode <- list(Pair_wise_site_comparisons_Fst = Pair_wise_site_comparisons_Fst, 
                      Pair_wise_site_comparisons_Fst_perSNP = Pair_wise_site_comparisons_Fst_perSNP) 
  return(Fst_barcode)
}


# =============================================================================
# Function to calculate peturbed estimates for significance of estimates
# =============================================================================
fst_Permuted <- function(SNPDataBinary){
  
  numSNPs <- ncol(SNPDataBinary)
  numSamples <- nrow(SNPDataBinary)
  permute <- function(x){x[sample(length(x))]}
  
  #----------------------------------------------------------------------------
  # Between sites averaging over years 
  #----------------------------------------------------------------------------
  foreach_return <- foreach(n = 1:n_repeat) %dopar% {
    
    # Allocate memory
    Pair_wise_site_comparisons_Fst <- array(dim = numComparisons, dimnames = list(geo_order))
    Location_true <- SNPData$City  # Preserve copy of unpermuted
    
    for(i in 1:numComparisons){
      
      # We only want to permute the labels of clinics A and B, not clincs A, B, C and D
      Location_temp <- Location_true
      ClinicA <- as.character(Pair_wise_site_comparisons[i,1])   
      ClinicB <- as.character(Pair_wise_site_comparisons[i,2])
      ClinicABInd <- Location_true == ClinicA | Location_true == ClinicB 
      Location_temp[ClinicABInd] <- permute(Location_true[ClinicABInd]) 
      
      # Indices for populations to compare
      P1 <- (Location_temp == as.character(Pair_wise_site_comparisons[i,1]))
      P2 <- (Location_temp == as.character(Pair_wise_site_comparisons[i,2]))
      P1_name <- as.character(Pair_wise_site_comparisons[i,1])
      P2_name <- as.character(Pair_wise_site_comparisons[i,2])
      
      # Calculate Fst
      X_Reich <- calculate_pairwise_reich(P1, P2, SNPDataBinary)
 
      # Save to matrix
      Pair_wise_site_comparisons_Fst[i] <- X_Reich$F_st
      rm(Location_temp)
    }
    
    Fst_Permuted <- list(Pair_wise_site_comparisons_Fst = Pair_wise_site_comparisons_Fst)
    return(Fst_Permuted)
  }
  
  # Re-structure result in an array, rather than in a list
  X <- do.call(abind, args = list(sapply(foreach_return, FUN = function(x){x['Pair_wise_site_comparisons_Fst']}), along = 2))
  Fst_Permuted <- list(Pair_wise_site_comparisons_Fst = X)
  return(Fst_Permuted)
}


# =============================================================================
# Function to calculate bootstrap estimates for confidence intervals
# =============================================================================
fst_bootstrapped <- function(SNPDataBinary){
  
  numSNPs <- ncol(SNPDataBinary)
  numSamples <- nrow(SNPDataBinary)
 
  # Between sites averaging over years -----------------------------------------
  foreach_return <- foreach(n = 1:n_repeat) %dopar% {
    
    # Allocate memory
    Pair_wise_site_comparisons_Fst <- array(dim = numComparisons, dimnames = list(geo_order))
      
    # Bootstrap SNPs
    bootstrap_snps <- sample(numSNPs, replace = TRUE)
    SNPDataBootStrap <- SNPDataBinary[, bootstrap_snps]
    
    for(i in 1:numComparisons){
      
      # Indices for populations to compare
      P1 <- (SNPData$City == as.character(Pair_wise_site_comparisons[i,1]))
      P2 <- (SNPData$City == as.character(Pair_wise_site_comparisons[i,2]))
      P1_name <- as.character(Pair_wise_site_comparisons[i,1])
      P2_name <- as.character(Pair_wise_site_comparisons[i,2])
      
      # Calculate pairwise thetas Weir and Cockerham
      X_Reich <- calculate_pairwise_reich(P1, P2, SNPDataBootStrap)
      
      # Save to matrix
      Pair_wise_site_comparisons_Fst[i] <- X_Reich$F_st
    }
    
    Fst_bootstrapped <- list(Pair_wise_site_comparisons_Fst = Pair_wise_site_comparisons_Fst)
    
    return(Fst_bootstrapped)
  }
  
  # Pull out CIs and re-structure result in an array
  X <- do.call(abind, args = list(sapply(foreach_return, FUN = function(x){x['Pair_wise_site_comparisons_Fst']}), along = 2))
  Fst_bootstrapped <- list(Pair_wise_site_comparisons_Fst = X)
  return(Fst_bootstrapped)
}


# =============================================================================
# Generate Fst estimates using functions above
# =============================================================================
# 0.3 sec
system.time(Fst_barcode <- fst_calculations(SNPDataBinary = SNPData[,c(6:255)]))

# Add Permuted (2.568 sec for 10)
system.time(Fst_barcode_Permuted <- fst_Permuted(SNPDataBinary = SNPData[,c(6:255)]))
Fst_barcode$Pair_wise_site_comparisons_Fst_Permuted <- Fst_barcode_Permuted$Pair_wise_site_comparisons_Fst

# Calculate CIs (2.569 sec for 10)
system.time(Fst_barcode_bootstrapped <- fst_bootstrapped(SNPDataBinary = SNPData[,c(6:255)]))
A <- Fst_barcode_bootstrapped$Pair_wise_site_comparisons_Fst
A_deltas <- apply(A, 2, FUN = function(x){x[geo_order] - Fst_barcode$Pair_wise_site_comparisons_Fst[geo_order]}) # Calculate differences
A_percentiles <- apply(A_deltas, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
Fst_barcode$Pair_wise_site_comparisons_Fst_CIs <- apply(A_percentiles, 1, FUN = function(x){Fst_barcode$Pair_wise_site_comparisons_Fst[geo_order] - x[geo_order]})

# Save in one big list
save(Fst_barcode, file = '~/Documents/BroadLaptop/ColombianBarcode/RData/Fst_barcode.RData')

