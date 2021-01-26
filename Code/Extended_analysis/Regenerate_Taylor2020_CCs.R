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
mle_CIs <- mle_CIs[!(mle_CIs$r2.5. < eps & mle_CIs$r97.5. > (1-eps)), ] # Remove uninformative

# Load metadata
load(sprintf('../../RData/metadata_extended.RData'))

# Load DE meta data file
DE_sids <- metadata$SAMPLE.CODE[metadata$PloSGen2020]
DE_rhat_ind <- (mle_CIs$individual1 %in% DE_sids & mle_CIs$individual2 %in% DE_sids)

#===========================================================
# Get clonal component (CC) membership for all samples
# Remove low quality samples to explore effect on CC
#===========================================================
# Create graph of all samples to extract CCs
A_high <- construct_adj_matrix(mle_CIs[DE_rhat_ind,], Entry = 'r97.5.') # Adj. matrix unfiltered using 'r97.5.' 
A_high[A_high < 1-eps] <- 0 # Edit s.t. only clonal have weight
G_high <- graph_from_adjacency_matrix(A_high, mode='upper', diag=F, weighted=T) # Construct graph 
C_high <- components(G_high) # Extract CCs from graph using cliques instead of components
M_high <- C_high$membership # Extract CC membership per sample (numeric name of each CC)

# Create a vector of all CCs with 2+ samples 
CCs <- as.character(which(C_high$csize > 1)) 

# Create a list with CC members for each cc
Clonal_components <- lapply(CCs, function(cc){
  names(which(M_high == as.numeric(cc)))
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
names(CC_chr_names) = as.character(M_high[sids_1stperCC_metadata$SAMPLE.CODE]) # Ensure CC names are ordered as memberships

# Rename Clonal_components
names(Clonal_components) <- CC_chr_names[CCs]

# Reorder Clonal_components and save
Clonal_components <- Clonal_components[CC_chr_names]
save(Clonal_components, file = "../../RData/Clonal_components_extended_DE_sid.RData")

rm(list = ls())
load(file = "../../RData/Clonal_components.RData")
CC_original <- Clonal_components
load(file = "../../RData/Clonal_components_extended_DE_sid.RData")
CC_extended <- Clonal_components

# Not exactly the same
identical(CC_original, CC_extended)
non_identical <- names(CC_original)[!sapply(names(CC_original), 
                                            function(CC){
                                              identical(CC_extended[[CC]], 
                                                        CC_original[[CC]])})]

# Difference of give and take two samples
lapply(non_identical, function(CC){ 
    c(extended = CC_extended[CC], original = CC_original[CC])})







