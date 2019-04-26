#################################
# To-do: could change to state
# Get colombian vivax
#################################
rm(list = ls())

# Load and summarise raw data 
load('../RData/mle_CIs.RData') 
load('../RData/edge_cols.RData') # for histograms partioned by clone
load('../RData/geo_dist_info_cities.RData')


# Extract sites in ascending order distance-wise
inter = geo_dist_info$geo_order  
intra = c('Guapi_Guapi', 'Tado_Tado', 'Tumaco_Tumaco', 'Buenaventura_Buenaventura', 
          'Quibdo_Quibdo')
intra_inter = c(intra, inter)

#=================================================================
# Calculate proportions by time and site comp
#=================================================================
nrep <- 100 # For bootstrap confidence intervals
r_threshold = 0.25
set.seed(1) # For reproducibility

# Add membership by colour
mle_CIs$sample_comp = apply(mle_CIs[, c("individual1", "individual2")], 1, function(x)paste(sort(x), collapse = "_"))
mle_CIs$comp_col = "#D3D3D3FF" # gray
if(all(names(edge_cols) %in% mle_CIs$sample_comp)){
  ind = mle_CIs$sample_comp %in% names(edge_cols) # reduce df to relevant comp of 
  for(comp in names(edge_cols)){
    ind_comp = which(mle_CIs$sample_comp[ind] == comp)
    mle_CIs$comp_col[ind][ind_comp] = edge_cols[comp]
  }} else {stop('Not all edges present')}

site_comps = intra_inter 
time_bins = sort(unique(mle_CIs$time_bins))
comp_cols = unique(mle_CIs$comp_col)
no_site_comps = length(site_comps)
no_time_bins = length(time_bins)
no_comp_cols = length(comp_cols)

# Stores (inc. those where intervals are partioned into site_comps and time bins)
proportions_time = array(0, dim = c(3, no_time_bins), dimnames = list(c('mean', '2.5%', '97.5%'), time_bins))
proportions_geo = array(0, dim = c(3, no_site_comps), dimnames = list(c('mean', '2.5%', '97.5%'), site_comps))
proportions_time_grouped = array(0, dim = c(no_site_comps, no_time_bins, 2), dimnames = list(site_comps, time_bins, c('all', 'r_threshold')))
proportions_geo_grouped = array(0, dim = c(no_time_bins, no_site_comps, 2), dimnames = list(time_bins, site_comps, c('all', 'r_threshold')))
proportions_time_cloned = array(0, dim = c(no_comp_cols, no_time_bins, 2), dimnames = list(comp_cols, time_bins, c('all', 'r_threshold')))
proportions_geo_cloned = array(0, dim = c(no_comp_cols, no_site_comps, 2), dimnames = list(comp_cols, site_comps, c('all', 'r_threshold')))

# Calculation proportions with time
for(i in time_bins){
  
  # Extract indices, data, etc. 
  index_i = as.character(i)
  ind <- mle_CIs$time_bins == i
  n_time_bin <- sum(ind)
  rhats <- mle_CIs$`2.5%`[ind]
  
  # Proportions of different site comps within
  site_comp_breakdown_all = table(mle_CIs$City12[ind])
  site_comp_breakdown_r_threshold = table(mle_CIs$City12[ind][rhats > r_threshold])
  proportions_time_grouped[names(site_comp_breakdown_all), index_i, 'all'] = site_comp_breakdown_all/n_time_bin
  proportions_time_grouped[names(site_comp_breakdown_r_threshold), index_i, 'r_threshold'] = site_comp_breakdown_r_threshold/n_time_bin
  
  # Proportions of different clonal components within
  clone_comp_breakdown_all = table(mle_CIs$comp_col[ind])
  clone_comp_breakdown_r_threshold = table(mle_CIs$comp_col[ind][rhats > r_threshold])
  proportions_time_cloned[names(clone_comp_breakdown_all), index_i, 'all'] = clone_comp_breakdown_all/n_time_bin
  proportions_time_cloned[names(clone_comp_breakdown_r_threshold), index_i, 'r_threshold'] = clone_comp_breakdown_r_threshold/n_time_bin
  
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
  rhats <- mle_CIs$`2.5%`[ind]
  
  # Proportions of different times within
  time_bin_breakdown_all = table(mle_CIs$time_bins[ind])
  time_bin_breakdown_r_threshold = table(mle_CIs$time_bins[ind][rhats > r_threshold])
  proportions_geo_grouped[names(time_bin_breakdown_all), i, 'all'] = time_bin_breakdown_all/n_city12
  proportions_geo_grouped[names(time_bin_breakdown_r_threshold), i, 'r_threshold'] = time_bin_breakdown_r_threshold/n_city12
  
  # Proportions of different clonal components within
  clone_comp_breakdown_all = table(mle_CIs$comp_col[ind])
  clone_comp_breakdown_r_threshold = table(mle_CIs$comp_col[ind][rhats > r_threshold])
  proportions_geo_cloned[names(clone_comp_breakdown_all), i, 'all'] = clone_comp_breakdown_all/n_city12
  proportions_geo_cloned[names(clone_comp_breakdown_r_threshold), i, 'r_threshold'] = clone_comp_breakdown_r_threshold/n_city12
  
  # Boostrap proportions 
  prob_b <- sapply(1:nrep, function(b)mean(sample(rhats, size = n_city12, replace = TRUE) > r_threshold))
  
  # Observed proportion 
  proportions_geo[,i] <- c('mean' =  mean(rhats > r_threshold), quantile(prob_b, probs = c(0.025, 0.975)))
}

# # Check proptions wrt time the same if grouped or not: Yes 
# cbind(colSums(proportions_geo_grouped[,,'r_threshold'], na.rm = TRUE), proportions_geo['mean',])

# CIs for intra clonal proportion
intra = c('Guapi_Guapi', 'Tado_Tado', 'Buenaventura_Buenaventura','Tumaco_Tumaco','Quibdo_Quibdo')
CIs_clonal_intra = sapply(intra_inter, function(i){
  
  # Extract data
  ind <- mle_CIs$City12 == i 
  n_city12 = sum(ind)
  clones <- mle_CIs$comp_col[ind]
  
  # Boostrap proportions not "#D3D3D3FF" (i.e. not singlton)
  prob_b <- sapply(1:nrep, function(b)mean(sample(clones, size = n_city12, replace = TRUE) != "#D3D3D3FF"))
  quantile(prob_b, probs = c(0.025, 0.975))
})
colnames(CIs_clonal_intra) = intra_inter


# Save 
save(proportions_time, proportions_time_grouped, proportions_time_cloned, file = '../RData/proportions_time.RData')
save(proportions_geo, proportions_geo_grouped, proportions_geo_cloned, file = '../RData/proportions_geo.RData')
save(CIs_clonal_intra, file = '../RData/CIs_clonal_intra.RData' )












