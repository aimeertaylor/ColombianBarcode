#################################################################################
# Script to investigate the impact of down sampling the barcode per site samples
# Takes ~ 21 seconds to run for loop
#################################################################################



# ============================================================================
# Function to calculate pairwise Fst estimates 
# ============================================================================
fst_calculations <- function(SNPDataBinary, Data,){
  
  SNPDataBinary <- SNPData[,c(6:255)]
  numSNPs <- ncol(SNPDataBinary)
  numSamples <- nrow(SNPDataBinary)
  
  # Between sites averaging over years ---------------------------------------
  Pair_wise_site_comparisons_Fst <- array(dim = c(numComparisons),dimnames = list(geo_order))
  Sample_size_min <- array(dim = c(numComparisons),dimnames = list(geo_order))
  
  for(i in 1:numComparisons){
    
    # Indices for populations to compare
    P1_ind <- which(SNPData$City == as.character(Pair_wise_site_comparisons[i,1]))
    P2_ind <- which(SNPData$City == as.character(Pair_wise_site_comparisons[i,2]))
    P1 <- 1:numSamples %in% sample(P1_ind, size = sample_down)
    P2 <- 1:numSamples %in% sample(P2_ind, size = sample_down)
    P1_name <- as.character(Pair_wise_site_comparisons[i,1])
    P2_name <- as.character(Pair_wise_site_comparisons[i,2])
    
    # Store sample size
    Sample_size_min[i] <- min(sum(P1), sum(P2))
    
    # Calculate pairwise thetas Weir and Cockerham
    X_Reich <- calculate_pairwise_reich(P1, P2, SNPDataBinary)
    
    # Save to matrix
    Pair_wise_site_comparisons_Fst[i, 'Reich'] <- X_Reich$F_st
  }
  
  return(Pair_wise_site_comparisons_Fst)
}

# ============================================================================
# Generate Fst estimates
# ============================================================================
sample_downs <- c(25,50,116)
Fst_barcode_dwn <- array(dim = c(100, 6, length(sample_downs)), dimnames = list(NULL, NULL, sample_downs))

system.time(
  for(j in sample_downs){
    for(i in 1:100){
      Fst_barcode_dwn[i,,as.character(j)] <- fst_calculations(SNPDataBinary = Data_store$SNPData_no_multiclonal,
                                                              Data = Data_store$Data_no_multiclonal,
                                                              sample_down = j)
    }
  }
)

save(Fst_barcode_dwn, file = '../../RData/Fst_barcode_dwn.RData')

