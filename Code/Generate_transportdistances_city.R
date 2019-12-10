##############################################################
# In this script we generate transport ("1-Wasserstein") distances 
# between parasite populations sampled in different cities.
#
# The transport approach is directly applicable to our between-city
# spatial analysis: it is trivial to construct a fully populated 
# adjacency matrix of comparisons across samples that fall into two 
# discrete populations (e.g. samples collected in Tumaco vesus 
# Buenaventura). 
#
# The transport approach is not directly applicable to our within-city
# spatial analysis: to construct an adjacency matrix of comparisons 
# across samples from a single city, we randomly split the set of 
# samples from that city into two. Since within city results are based 
# on a random split they are random and included as a sanity check only. 
##############################################################
rm(list = ls())
library(transport) # transport package
load('../RData/All_results.RData')
load('../RData/geo_dist_info.RData')
attach(geo_dist_info)
source('./transport_functions.R')
nrepeats = 100 # Number of bootstrap repeats (79 seconds)
set.seed(10) # For reproduciblity within sites

# Get inter and intra city names
inter_c = geo_dist_info$geo_order  
intra_c = c('Guapi_Guapi', 
            'Tado_Tado', 
            'Tumaco_Tumaco', 
            'Buenaventura_Buenaventura', 
            'Quibdo_Quibdo')
intra_inter_c = c(intra_c, inter_c)

All_W_results_c = lapply(All_results, function(mle_CIs){
  
  # Change characters to factors
  for(j in 1:ncol(mle_CIs)){
    if(class(mle_CIs[,j]) == 'factor'){mle_CIs[,j] = as.character(mle_CIs[,j])}
  }
  
  # Generate a full rhat adjacency matrix
  sample_names = unique(c(mle_CIs$individual1, mle_CIs$individual2))
  adj_matrix = array(data = 0, dim = rep(length(sample_names), 2), 
                     dimnames = list(sample_names, sample_names))
  
  for(i in 1:nrow(mle_CIs)){
    indi = mle_CIs$individual1[i]
    indj = mle_CIs$individual2[i]
    adj_matrix[indi, indj] = adj_matrix[indj, indi] = 1-mle_CIs[i, 'rhat']
  }
  
  #=============================================================================
  # Calculate the "1-Wasserstein" distance within and between cities
  #============================================================================
  W_results_c = sapply(intra_inter_c, function(city_comp){
    
    cities = strsplit(city_comp, split = '_')[[1]] # Extract city names
    inds = (mle_CIs$City1 %in% cities) & (mle_CIs$City2 %in% cities) # indices for cities
    mle_cities = mle_CIs[inds,] # Extract relevent mles
    
    if (cities[1] == cities[2]){
      
      # Split the sample into two, uniformly at random
      samples = unique(c(mle_cities$individual1, mle_cities$individual2))
      samples_c1 = sample(samples, size = floor(length(samples)/2), replace = F)
      samples_c2 = setdiff(samples, samples_c1)
      
    } else {
      
      # Samples per different city
      samples_c1 = unique(c(mle_cities$individual1[mle_cities$City1 %in% cities[1]], 
                            mle_cities$individual2[mle_cities$City2 %in% cities[1]]))
      samples_c2 = unique(c(mle_cities$individual1[mle_cities$City1 %in% cities[2]], 
                            mle_cities$individual2[mle_cities$City2 %in% cities[2]]))
    }
    
    # Extravt distance matrix
    dist_matrix = adj_matrix[samples_c1, samples_c2]
    
    # Caclulate 1-W and CIs
    cost <- generate_1_w(dist_matrix)
    cost_CI <- generate_1_w_CI(dist_matrix, nrepeats)
    return(c(W = cost, cost_CI))
  })
  colnames(W_results_c) = intra_inter_c
  return(W_results_c)
})


# Save
save(All_W_results_c, file = '../RData/All_W_results_c.RData')

# Quick plot
W_results_c = All_W_results_c$Unfiltered
X = barplot(W_results_c[1,], las = 2, ylim = c(0, max(W_results_c)))
segments(x0 = X, x1 = X, y0 = W_results_c[2,], y1 = W_results_c[3,])









