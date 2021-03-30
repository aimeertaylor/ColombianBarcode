##############################################################################
#' Script to compare the clonal components reported in Taylor et al. 2020 with
#' those generated after reordering markers that were previously missordered
#' alphabetically. In summary, 45 instead of 46 CCs. Among 46, zero broken, 43
#' identical, old CC6 and CC19 grouped together in new CC6 and one sample added
#' to old CC25 in new CC24. Does this affect the paper, which hones in on the
#' five CCs with samples from both Buenaventura and Tumaco: old 1 (new 1), old
#' 12 (new 12), old 14 (new 14), old 20 (new 19) & old 40 (new 39). 
##############################################################################
rm(list = ls())

load(file = sprintf("../RData/PlosGen2020/Clonal_components.RData"))
CC_original <- Clonal_components

load(file = "../RData/Clonal_components.RData")
CC_reordered <- Clonal_components
range(sapply(CC_reordered, length)) # Size range 

load(file = "../RData/metadata_extended.RData") # To check number of CCs detected in more than one city
sum(sapply(CC_reordered, function(ccs) length(unique(metadata[ccs,"City"]))) > 1)


length(CC_reordered)
length(CC_original)
all(unlist(CC_original) %in% unlist(CC_reordered)) # Check all are within

# Define stores to categorise clonal components reported in Taylor et al. 2020
identical_components <- list()
extended_components <- list()
broken_components <- list()
  
# Cartegorise each of the clonal components reported in Taylor et al. 2020
for(i in names(CC_original)){
  for(j in names(CC_reordered)){
    if(setequal(CC_original[[i]],CC_reordered[[j]])) {
      identical_components[[i]][[j]] <- "Equal"
    } else if (all(CC_original[[i]] %in% CC_reordered[[j]])) {
      extended_components[[i]][[j]][["original"]] <- CC_original[[i]]
      extended_components[[i]][[j]][["additional"]] <- setdiff(CC_reordered[[j]],CC_original[[i]])
    } else {
      intersecting_sids <- intersect(CC_original[[i]],CC_reordered[[j]])
      if(length(intersecting_sids) > 0) {
        broken_components[[i]][[j]][["intersect"]] <- intersecting_sids
        broken_components[[i]][[j]][["additional"]] <- setdiff(CC_reordered[[j]],CC_original[[i]])
      } else {
        next()
      }
    }
  }
}

# Category breakdown
length(broken_components) # Zero broken
length(identical_components) #  44 identical; CC23 broken into singletons
length(extended_components) # CC25 extended

# Inspect
identical_components_key <- sapply(identical_components, names)
identical_components_key
broken_components 
extended_components 

# Identical but newly labeled
identical_components_key[identical_components_key != names(identical_components_key)]

