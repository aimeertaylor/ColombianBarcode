##############################################
# This script simply loads a very large file 
# and saves more manageable versions
##############################################
rm(list = ls())
load('../RData/mles.RData') # Few mins
mle_core = mle_df[,1:4] # Extract core

# Save a more manageable version withs only 
# 100 rhat bootstrapped quantities 
rhats_boot_100 = lapply(mle_df$rhats_boot, function(x) x[1:100])
mle_core$rhats_boot = rhats_boot_100
mle_boot = mle_core
save(mle_boot, file = '../RData/mle_boot.RData')

#---------------------------------------------
# Save a version withs rhat CIs 
#---------------------------------------------
P = c(0.025, 0.975) # CI quantiles
rhat_CIs = t(apply(mle_df, 1, function(x)quantile(x$rhats_boot, probs=P))) # Few seconds
mle_CIs = cbind(mle_core, rhat_CIs)

# Do some formating 
if(class(mle_CIs$individual1) == 'factor'){mle_CIs$individual1 = as.character(mle_CIs$individual1)}
if(class(mle_CIs$individual2) == 'factor'){mle_CIs$individual2 = as.character(mle_CIs$individual2)}

# Load data to add meta data 
load('../RData/SNPData.RData') 
load('../RData/geo_dist_info.RData')

# Add site comps
mle_CIs$City1 = as.character(SNPData[mle_CIs$individual1, 'City'])
mle_CIs$City2 = as.character(SNPData[mle_CIs$individual2, 'City'])
mle_CIs$City12 = apply(mle_CIs[,c('City1','City2')], 1, function(x)paste(sort(x),collapse="_"))
mle_CIs$geo_dist = geo_dist_info$pairwise_site_distance_all[mle_CIs$City12]

# Add site comps
mle_CIs$date1 = as.character(SNPData[mle_CIs$individual1, 'COLLECTION.DATE'])
mle_CIs$date2 = as.character(SNPData[mle_CIs$individual2, 'COLLECTION.DATE'])
mle_CIs$time_dist = abs(difftime(mle_CIs$date1, mle_CIs$date2, units = 'weeks'))

# Add sample comp
mle_CIs$sample_comp = apply(mle_CIs[, c("individual1", "individual2")], 1, 
                            function(x)paste(sort(x), collapse = "_"))

# Save data frame
save(mle_CIs, file = '../RData/mle_CIs.RData')
