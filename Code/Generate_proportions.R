#######################################################
# This script generates proportions 
# given a single example threshold (e.g. 0.25), 
# using cities and time as the partion
#######################################################
rm(list = ls())

# Load and summarise raw data 
load('../RData/All_results.RData') 
load('../RData/geo_dist_info_cities.RData')
mle_CIs <- All_results$Unfiltered

# Extract sites in ascending order distance-wise
inter = geo_dist_info$geo_order  
intra = c('Guapi_Guapi', 
          'Tado_Tado', 
          'Tumaco_Tumaco', 
          'Buenaventura_Buenaventura', 
          'Quibdo_Quibdo')
intra_inter = c(intra, inter)

#=================================================================
# Calculate proportions by time and site comp
#=================================================================
nrep <- 100 # For bootstrap confidence intervals
r_threshold = 0.25
set.seed(1) # For reproducibility 

site_comps = intra_inter 
time_bins = sort(unique(mle_CIs$time_bins))
no_site_comps = length(site_comps)
no_time_bins = length(time_bins)

# Stores (inc. those where intervals are partioned into site_comps and time bins)
proportions_time = array(0, dim = c(3, no_time_bins), dimnames = list(c('mean', '2.5%', '97.5%'), time_bins))
proportions_geo = array(0, dim = c(3, no_site_comps), dimnames = list(c('mean', '2.5%', '97.5%'), site_comps))
proportions_time_grouped = array(0, dim = c(no_site_comps, no_time_bins, 2), dimnames = list(site_comps, time_bins, c('all', 'r_threshold')))
proportions_geo_grouped = array(0, dim = c(no_time_bins, no_site_comps, 2), dimnames = list(time_bins, site_comps, c('all', 'r_threshold')))

# Calculation proportions with time
for(i in time_bins){
  
  # Extract indices, data, etc. 
  index_i = as.character(i)
  ind <- mle_CIs$time_bins == i
  n_time_bin <- sum(ind)
  rhats <- mle_CIs$`r2.5.`[ind]
  
  # Proportions of different site comps within
  site_comp_breakdown_all = table(mle_CIs$City12[ind])
  site_comp_breakdown_r_threshold = table(mle_CIs$City12[ind][rhats > r_threshold])
  proportions_time_grouped[names(site_comp_breakdown_all), index_i, 'all'] = site_comp_breakdown_all/n_time_bin
  proportions_time_grouped[names(site_comp_breakdown_r_threshold), index_i, 'r_threshold'] = site_comp_breakdown_r_threshold/n_time_bin
  
  # Boostrap proportions wrt time
  prob_b <- sapply(1:nrep, function(b)mean(sample(rhats, size = n_time_bin, replace = TRUE) > r_threshold))
  
  # Observed proportion overall wrt time and quantiles
  proportions_time[,index_i] = c(mean(rhats > r_threshold), quantile(prob_b, probs = c(0.025, 0.975)))
}

# # Check proptions wrt time the same if grouped or not: Yes 
# cbind(colSums(proportions_time_grouped[,,'r_threshold'], na.rm = TRUE),proportions_time['mean',])

# Calculate proportions with site
for(i in site_comps){
  
  # Extract data
  ind <- mle_CIs$City12 == i 
  n_city12 = sum(ind)
  rhats <- mle_CIs$`r2.5.`[ind]
  
  # Proportions of different times within
  time_bin_breakdown_all = table(mle_CIs$time_bins[ind])
  time_bin_breakdown_r_threshold = table(mle_CIs$time_bins[ind][rhats > r_threshold])
  proportions_geo_grouped[names(time_bin_breakdown_all), i, 'all'] = time_bin_breakdown_all/n_city12
  proportions_geo_grouped[names(time_bin_breakdown_r_threshold), i, 'r_threshold'] = time_bin_breakdown_r_threshold/n_city12
  
  # Boostrap proportions 
  prob_b <- sapply(1:nrep, function(b)mean(sample(rhats, size = n_city12, replace = TRUE) > r_threshold))
  
  # Observed proportion 
  proportions_geo[,i] <- c('mean' =  mean(rhats > r_threshold), quantile(prob_b, probs = c(0.025, 0.975)))
}

# # Check proptions wrt time the same if grouped or not: Yes 
# cbind(colSums(proportions_geo_grouped[,,'r_threshold'], na.rm = TRUE), proportions_geo['mean',])

# Save 
save(proportions_time, proportions_time_grouped, file = '../RData/proportions_time.RData')
save(proportions_geo, proportions_geo_grouped, file = '../RData/proportions_geo.RData')












