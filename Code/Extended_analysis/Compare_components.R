##############################################################################
#' Script to compare the clonal components generated Generate_components.R. To
#' see how FSVC sids cluster with DE sids, we compare the clonal components
#' reported in Taylor et al. 2020 to those generated using all sample ids. To
#' see clonal components made using 250 SNP data compare with clusters based on
#' WGS data, we compare clusters geneterated from FSVC sids only against cluster
#' defined using WGS data (TO-DO).
#' 
#' The LCI threshold makes no difference to extended_components. It just breaks 
#' Up the original 46. 
##############################################################################


#=============================================================================
#' To compare the clonal components reported in Taylor et al. 2020 to those
#' generated using all sample ids.
#=============================================================================
rm(list = ls())
load(file = "../../RData/Clonal_components.RData")
CC_original <- Clonal_components

LCI_threshold <- 0.95 # Toggle between 0, 0.75 and 0.95 (make no difference)
load(file = sprintf("../../RData/Clonal_components_extended_all_LCIthrehold_%s.RData", LCI_threshold))
CC_extended <- Clonal_components

length(CC_extended)
length(CC_original)

# Define stores to categorise clonal components reported in Taylor et al. 2020
identical_components <- list()
extended_components <- list()
broken_components <- list()
  
# Cartegorise each of the clonal components reported in Taylor et al. 2020
for(i in names(CC_original)){
  for(j in names(CC_extended)){
    if(setequal(CC_original[[i]],CC_extended[[j]])) {
      identical_components[[i]][[j]] <- "Equal"
    } else if (all(CC_original[[i]] %in% CC_extended[[j]])) {
      extended_components[[i]][[j]][["original"]] <- CC_original[[i]]
      extended_components[[i]][[j]][["additional"]] <- setdiff(CC_extended[[j]],CC_original[[i]])
    } else {
      intersecting_sids <- intersect(CC_original[[i]],CC_extended[[j]])
      if(length(intersecting_sids) > 0) {
        broken_components[[i]][[j]][["intersect"]] <- intersecting_sids
        broken_components[[i]][[j]][["additional"]] <- setdiff(CC_extended[[j]],CC_original[[i]])
      } else {
        next()
      }
    }
  }
}

length(broken_components)
broken_components 

length(identical_components) 
identical_components_key <- sapply(identical_components, names)

length(extended_components) 
extended_components 

load(sprintf('../../RData/metadata_extended.RData'))
extended_components_metadata <- lapply(extended_components, function(x){
  metadata[unlist(x), ]
})





