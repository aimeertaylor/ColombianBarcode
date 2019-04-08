#############################################
# This script generates proportions of 
# highly related parasite sample comparisons
# using various thresholds for high relatedness
# with and without filtering 
#############################################
rm(list = ls())
load('../RData/All_results.RData') 

#-----------------------------------------------
# Fractions highly related - to finish
#-----------------------------------------------
cities = unique(All_results[["Unfiltered"]]$City12)
states = unique(All_results[["Unfiltered"]]$State12)
times = unique(All_results[["Unfiltered"]]$time_bins)
Summaries = c('mean', '2.5%', '97.5%') # Quanties of interest
Thresholds = c(0.01, 0.25, 0.5)
nrep = 100 # For bootstrap confidence intervals
set.seed(1) # For reproducibility

# Stores 
proportions_cities = array(0, dim = c(length(Summaries),length(cities),length(Thresholds),length(All_results)), 
                                     dimnames = list(Summaries,cities,Thresholds,names(All_results)))
proportions_states = array(0, dim = c(length(Summaries),length(states),length(Thresholds),length(All_results)), 
                           dimnames = list(Summaries,states,Thresholds,names(All_results)))
proportions_times = array(0, dim = c(length(Summaries),length(times),length(Thresholds),length(All_results)), 
                                     dimnames = list(Summaries,times,Thresholds,names(All_results)))

# Function to retuen mean and CIs given inds
return_mean_CIs = function(inds){
  Mean = mean(X$`2.5%`[inds] >= threshold)
  prop_bootstrapped = sapply(1:nrep, function(b){
    booti_ind = sample(inds, length(inds), replace = T)
    mean(X$`2.5%`[booti_ind] >= threshold)})
  CIs = quantile(prop_bootstrapped, probs = c(0.025, 0.975))
  return(c('mean' = Mean, CIs))
}


for(l in names(All_results)){
  
  X = All_results[[l]]
  
  for(threshold in Thresholds){
    
    # Calculate proportions over city comparisons 
    prop_cities = sapply(cities, function(x){
      return_mean_CIs(which(X$City12 == x))
    })
    
    # Calculate proportions over state comparisons 
    prop_states = sapply(states, function(x){
      return_mean_CIs(which(X$State12 == x))
    })
    
    # Calculate proportions over time comparisons 
    prop_times = sapply(as.character(times), function(x){
      return_mean_CIs(which(X$time_bins == as.numeric(x)))
    })
    
    proportions_cities[,,as.character(threshold), l] = prop_cities
    proportions_states[,,as.character(threshold), l] = prop_states
    proportions_times[,,as.character(threshold), l] = prop_times
}}

save(proportions_cities, proportions_states, proportions_times, 
     file = '../RData/proportions_sensitivities')


