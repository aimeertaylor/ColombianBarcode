##############################################################
# Script to re-generate Taylor et al CC membership for checking
# purposes only 
##############################################################
rm(list = ls())
library(igraph) # To make graph, construct components etc.
source('../igraph_functions.R')
eps <- 0.01 # threshold below which LCI is considered zero in Taylor et al. 2019

# Load relatedness results
load('../../RData/mles_CIs_extended_freqsTaylor2020_meta.RData') 
mle_CIs <- mle_CIs[!is.na(mle_CIs$rhat), ] # remove NAs 

# Load metadata
load(sprintf('../../RData/metadata_extended.RData'))

# Load DE meta data file
DE_sids <- metadata$SAMPLE.CODE[metadata$PloSGen2020]
DE_rhat_ind <- (mle_CIs$individual1 %in% DE_sids & mle_CIs$individual2 %in% DE_sids)
any(is.na(mle_CIs[DE_rhat_ind,"rhat"])) # Check for NAs before using graph_from_adjacency_matrix())

#===========================================================
# Get clonal component (CC) membership for all samples using 
# 1) the definition of a clonal used in Taylor et al. 2020
# 2) a more stringent definition of a clone 
#===========================================================
A_low <- construct_adj_matrix(mle_CIs[DE_rhat_ind,], Entry = 'r2.5.')  
A_high <- construct_adj_matrix(mle_CIs[DE_rhat_ind,], Entry = 'r97.5.') 
A_est <- construct_adj_matrix(mle_CIs[DE_rhat_ind,], Entry = 'rhat') 

for(definition in 1:2){
  
  if(definition == 1) {
    A_clonal <- A_high
    A_clonal[A_high < 1-eps] <- 0
    print(min(A_low[A_clonal > 0], na.rm = T)) 
  }
  else {
    A_clonal <- array(NA, dim = dim(A_high), dimnames = dimnames(A_high))
    A_clonal[A_high >= (1-eps) & A_low > 0.75] <- 1
    A_clonal[A_high < (1-eps) | A_low < 0.75] <- 0
  }
  
  # Check range
  print(range(A_clonal[A_clonal > 0], na.rm = T)) # Remove NA as A_clonal upper tri only 
  
  G_clonal <- graph_from_adjacency_matrix(A_clonal, mode='upper', diag=F, weighted=T) # Construct graph 
  C_clonal <- components(G_clonal) # Extract CCs from graph using cliques instead of components
  M_clonal <- C_clonal$membership # Extract CC membership per sample (numeric name of each CC)
  
  # Create a vector of all CCs with 2+ samples 
  CCs <- as.character(which(C_clonal$csize > 1))
  
  # Create a list with CC members for each cc
  Clonal_components <- lapply(CCs, function(cc){
    names(which(M_clonal == as.numeric(cc)))
  })
  
  # Extract the sid of the earliest detected sample per cc
  sids_1stperCC <- t(sapply(Clonal_components, function(cc){
    earliest_date <- min(metadata[cc,"COLLECTION.DATE"])
    sids <- cc[which(metadata[cc,"COLLECTION.DATE"] == earliest_date)]
    sid <- sids[1] # Take the first if there is more than one
  }))
  
  # Extract metadata, ordered by date, for the each of the earliest detected sample per cc
  sids_1stperCC_metadata <- dplyr::arrange(metadata[sids_1stperCC,], COLLECTION.DATE)
  
  # Create "CC#" character name order by the earliest detected sample per CC
  CC_chr_names <- paste0('CC',1:length(CCs))
  names(CC_chr_names) = as.character(M_clonal[sids_1stperCC_metadata$SAMPLE.CODE]) # Ensure CC names are ordered as memberships
  
  # Rename Clonal_components
  names(Clonal_components) <- CC_chr_names[CCs]
  
  # Reorder Clonal_components and save
  Clonal_components <- Clonal_components[CC_chr_names]
  save(Clonal_components, 
       file = sprintf("../../RData/Clonal_components_regenerated_definition%s.RData", definition))
}




#==============================================================
# Part two: check newly generated clonal components against old
#==============================================================
rm(list = ls())
load(file = "../../RData/Clonal_components.RData")
CC_original <- Clonal_components

#------------- Definition two ---------------
for(definition in 1:2){
  
  load(file = sprintf("../../RData/Clonal_components_regenerated_definition%s.RData", definition))
  CC_extended <- Clonal_components
  # Different CC counts: cannot attempt to align and compare
  
  # let's attempt to align and compare
  if(length(CC_extended) == length(CC_original)){
    
    # Definition 1) not exactly the same, but only two differences
    identical(CC_original, CC_extended)
    non_identical <- names(CC_original)[!sapply(names(CC_original), 
                                                function(CC){
                                                  identical(CC_extended[[CC]], 
                                                            CC_original[[CC]])})]
    
    # Definition 1) difference of give and take two samples
    lapply(non_identical, function(CC){ 
      c(extended = CC_extended[CC], original = CC_original[CC])})
    
  } else {
    
    # Do some overlap? Yes
    component_intersects <- list()
    for(i in names(CC_original)){
      for(j in names(CC_extended)){
        x <- intersect(CC_original[[i]],CC_extended[[j]]) 
        if(length(x) > 0) {
          component_intersects[[i]][[j]] <- x
        } else {
          next()
        }
      }
    }
    
    # Mapping of regenerated onto original
    mapping <- lapply(component_intersects, names)
    sapply(mapping, length) # mapping is one-to-one
    unlist(mapping) # CC10, CC17, CC23, CC32 original missing
  }
}






