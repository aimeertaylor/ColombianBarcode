#################################################################
# This script 
# 1) loads the output from Generate_mles_CIs_inc.GuapiWGStoBarcode.R and adds metadata 
# 2) loads the mles with add data and filter
#################################################################

#====================================================
# 1) Add metadata 
#====================================================
rm(list = ls())
load("../../RData/mles_CIs_extended.RData") 

# Factors to characters
if(class(mle_CIs$individual1) == 'factor'){mle_CIs$individual1 = as.character(mle_CIs$individual1)}
if(class(mle_CIs$individual2) == 'factor'){mle_CIs$individual2 = as.character(mle_CIs$individual2)}

# Load Diego's meta data 
load('../../RData/SNPDataRecode.RData') # These are Diego's data
load('../../RData/geo_dist_info_cities.RData')
geo_dist_info_cities = geo_dist_info
load('../../RData/geo_dist_info_states.RData')
geo_dist_info_states = geo_dist_info

# Extract Guapi WGS to Barcode sample names
all_sids = unique(c(mle_CIs$individual1, mle_CIs$individual2))
VC_sids = all_sids[!all_sids %in% rownames(SNPDataRecode)]

# Create a template SNPData for VC_sids
SNPData_VC <- data.frame(array(NA, dim = c(length(VC_sids), ncol(SNPDataRecode)), 
                               dimnames = list(VC_sids, colnames(SNPDataRecode))), 
                         check.names = F) # Prevents conversion of "-" to "."
SNPData_VC$City = "Guapi"
SNPData_VC$STATE = "Cauca"
SNPData_VC$COLLECTION.DATE = format(Sys.Date(), format = "%Y-%m-%d")
SNPData_VC$Year = format(Sys.Date(), format = "%Y")

# Add SNPData_VC to SNPData
SNPData = rbind(SNPDataRecode, SNPData_VC)

# Add city comps
mle_CIs$City1 = as.character(SNPData[mle_CIs$individual1, 'City'])
mle_CIs$City2 = as.character(SNPData[mle_CIs$individual2, 'City'])
mle_CIs$City12 = apply(mle_CIs[,c('City1','City2')], 1, function(x)paste(sort(x),collapse="_"))

# Add state comps
mle_CIs$State1 = as.character(SNPData[mle_CIs$individual1, 'STATE'])
mle_CIs$State2 = as.character(SNPData[mle_CIs$individual2, 'STATE'])
mle_CIs$State12 = apply(mle_CIs[,c('State1','State2')], 1, function(x)paste(sort(x),collapse="_"))

# Add distance between cities
mle_CIs$geo_dist_cities = geo_dist_info_cities$pairwise_site_distance_all[mle_CIs$City12]
mle_CIs$geo_dist_states = geo_dist_info_states$pairwise_site_distance_all[mle_CIs$State12]

# Add dates and distance in time 
mle_CIs$date1 = as.character(SNPData[mle_CIs$individual1, 'COLLECTION.DATE'])
mle_CIs$date2 = as.character(SNPData[mle_CIs$individual2, 'COLLECTION.DATE'])
mle_CIs$time_dist = abs(difftime(mle_CIs$date1, mle_CIs$date2, units = 'weeks'))

# Add distance in time bins
interval <- 50 # The resolution of time bins in weeks
max_time_dist = max(as.numeric(mle_CIs$time_dist))
time_breaks <- seq(0, max_time_dist, interval) 
max_time_breaks = max(time_breaks)
if(max_time_breaks < max_time_dist){time_breaks = c(time_breaks, max_time_breaks + interval)}
time_bins = time_breaks[-length(time_breaks)] # The time bins value are the min for each bin
no_time_bins = length(time_bins)

# Add time bins 
mle_CIs$time_bins = NA
for(i in 1:no_time_bins){
  ind <- (as.numeric(mle_CIs$time_dist) >= time_breaks[i]) & (as.numeric(mle_CIs$time_dist) < time_breaks[i+1])
  mle_CIs$time_bins[ind] = time_bins[i]
}

# Add sample comp
mle_CIs$sample_comp = apply(mle_CIs[, c("individual1", "individual2")], 1, 
                            function(x)paste(sort(x), collapse = "_"))

# Save data frame
save(mle_CIs, file = "../../RData/mles_CIs_extended_meta.RData")
save(SNPData, file = "../../RData/SNPData_extended.RData")


#====================================================
# 2) Filter data set with added metadata, save and 
# delete intermediate file
#====================================================
rm(list = ls())
library(igraph)
source('../igraph_functions.R') # For rm_highly_related_within and construct_adj_matrix
load('../../RData/SNPData_extended.RData') # Load SNP data for cities
eps = 0.01 # Below which LCI considered close to zero (needed by rm_highly_related_within)
Cities = SNPData$City; names(Cities) = row.names(SNPData) # n x 1 vector of cities named by sample ID
load('../../RData/mles_CIs_extended_meta.RData') # Few mins

#===========================================================
# Filter and save results
#   - remove edges (almost all samples remain)
#   - remove vertices (removes all samples per CC per city except one)
#===========================================================
All_results = lapply(c(F,T), rm_highly_related_within, Result = mle_CIs, Cities = Cities)
All_results[[3]] = mle_CIs # Add unfiltered
names(All_results) = c('Filter by vertex', 'Filter by edge', 'Unfiltered')

# Extract summaries
EV_summary = sapply(All_results, function(x){c('Edge count' = nrow(x), 
                                               'Vertex count' = length(unique(c(x$individual1,x$individual2))))})

writeLines('Data set summary of edge and vertex count:')
print(EV_summary)

system('rm ../../RData/mles_CIs_extended_meta.RData') # Delete intermediate file
save(All_results, file = '../../RData/All_results_extended.RData')



