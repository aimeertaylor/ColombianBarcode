##############################################################
# In this script we generate and plot transport distances between 
# parasite populations sampled from different cities on the 
# Colombian Pacific coast
# Rename
##############################################################

rm(list = ls())
library(transport) # transport package
load('../RData/All_results.RData')
load('../RData/SNPData.RData')
load('../RData/geo_dist_info_cities.RData')
attach(geo_dist_info)
nrepeats = 100 # Number of bootstrap repeats (79 seconds)
Gen_dist = T # Set to true to regenerate distances, otherwise just plot

#==========================================================
# Function to create an adjacency matrix
#=========================================================
construct_adj_matrix = function(Result, Entry = 'rhat', dm){
  
  # Cities 
  cities = unique(c(Result$City1,Result$City2))
  
  # Samples per city
  samples_c1 = unique(c(Result$individual1[Result$City1 %in% cities[1]], 
                        Result$individual2[Result$City2 %in% cities[1]]))
  samples_c2 = unique(c(Result$individual1[Result$City1 %in% cities[2]], 
                        Result$individual2[Result$City2 %in% cities[2]]))
  
  sample_names = c(samples_c1, samples_c2) # City 1 then city 2
  sample_count = length(sample_names)
  
  adj_matrix = array(data = 0, dim = c(sample_count, sample_count), 
                     dimnames = list(sample_names, sample_names))
  
  for(i in 1:nrow(Result)){
    indi = Result$individual1[i]
    indj = Result$individual2[i]
    adj_matrix[indi, indj] = adj_matrix[indj, indi] = 1-Result[i, Entry]
  }
  
  if(dm == T){ # If distance matrix true
    to_return = adj_matrix[samples_c1, samples_c2]
  } else {
    to_return = list(adj_matrix = adj_matrix, 
                     samples_c1 = samples_c1,  
                     samples_c2 = samples_c2)
  }
  return(to_return)
}


system.time(if(Gen_dist){
  All_W_results = lapply(All_results, function(mle_CIs){
    
    if(class(mle_CIs$individual1) == 'factor'){mle_CIs$individual1 = as.character(mle_CIs$individual1)}
    if(class(mle_CIs$individual2) == 'factor'){mle_CIs$individual2 = as.character(mle_CIs$individual2)}
    
    #=============================================================================
    # Calculate the "1-Wasserstein" distance between population 1 and population 2
    #=============================================================================
    W_results = sapply(geo_order, function(city_comp){
      
      cities = strsplit(city_comp, split = '_')[[1]]
      inds = (mle_CIs$City1 %in% cities) & (mle_CIs$City2 %in% cities) # indices for cities
      
      # Actual distance matrix
      dist_matrix = construct_adj_matrix(Result = mle_CIs[inds, ], dm = T)
      
      # Bootstrapping rows then columns (samples from citites 1 and 2) 
      dist_matrix_boot = lapply(1:nrepeats, function(b){
        row_boot = sample(nrow(dist_matrix), nrow(dist_matrix), replace = T)
        col_boot = sample(ncol(dist_matrix), ncol(dist_matrix), replace = T)
        dm_boot = rbind(dist_matrix[row_boot,])
        dm_boot = cbind(dm_boot[row_boot,])
        return(dm_boot)
      })
      
      # Weights (could adapted to propogate uncertainity)
      w1 <- rep(1/nrow(dist_matrix), nrow(dist_matrix))
      w2 <- rep(1/ncol(dist_matrix), ncol(dist_matrix))
      
      # Caclulate 1-W
      a <- transport(w1, w2, costm = dist_matrix, method = "shortsimplex")
      cost <- 0 
      for (i in 1:nrow(a)){
        cost <- cost + dist_matrix[a$from[i], a$to[i]] * a$mass[i]}
      
      # Calculate boostrapped 1-W
      costs_boot = sapply(dist_matrix_boot, function(dist_matrix){
        a <- transport(w1, w2, costm = dist_matrix, method = "shortsimplex")
        cost <- 0
        for (i in 1:nrow(a)){
          cost <- cost + dist_matrix[a$from[i], a$to[i]] * a$mass[i]}
        return(cost)})
      
      # Return value and CIs
      c('cost' = cost, quantile(costs_boot, probs = c(0.025,0.975)))
    })
    colnames(W_results) = as.character(geo_order)
    return(W_results)
  })
  save(All_W_results, file = '../RData/All_W_results.RData')
})



# ### now Sinkhorn distance implemented in the winference package
# library(devtools)
# devtools::install_github("pierrejacob/winference") # Error :-(
# library(winference)
# 
# C <- dist_matrix
# eps <- 0.01
# p <- 1
# epsilon <- eps * median(C^p)
# wass <- winference::wasserstein(w1, w2, C^p, epsilon, 1e3)
# wass$distances
# ## wass$distances above is the "1-Sinkhorn" divergence between population 1 and population 2





