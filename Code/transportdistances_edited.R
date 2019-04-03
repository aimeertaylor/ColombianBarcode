##############################################################
# To-do
# Try with Pierre's different measures
##############################################################

#rm(list = ls())
library(transport) # transport package
library(snowboot)
load('../RData/mle_CIs.RData')
load('../RData/SNPData.RData')
load('../RData/geo_dist_info.RData')

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



par(mfrow = c(2,2), family = 'serif', mar = c(6,4,3,1))
for(j in 1:length(All_results)){
  
  mle_CIs = All_results[[j]]
  if(class(mle_CIs$individual1) == 'factor'){mle_CIs$individual1 = as.character(mle_CIs$individual1)}
  if(class(mle_CIs$individual2) == 'factor'){mle_CIs$individual2 = as.character(mle_CIs$individual2)}
  
  # Add site comps
  mle_CIs$City1 = as.character(SNPData[mle_CIs$individual1, 'City'])
  mle_CIs$City2 = as.character(SNPData[mle_CIs$individual2, 'City'])
  mle_CIs$City12 = apply(mle_CIs[,c('City1','City2')], 1, function(x)paste(sort(x),collapse="_"))
  
  # Create adjacency matrix for each combination of cities
  costs = array(NA, length(geo_dist_info$geo_order), dimnames = list(geo_dist_info$geo_order))
  
  ## Calculate the "1-Wasserstein" distance between population 1 and population 2
  for(city_comp in geo_dist_info$geo_order){
    cities = strsplit(city_comp, split = '_')[[1]]
    inds = (mle_CIs$City1 %in% cities) & (mle_CIs$City2 %in% cities) # indices for cities
    #inds = mle_CIs$City12 == city_comp
    dist_matrix = construct_adj_matrix(Result = mle_CIs[inds, ], dm = T)
    w1 <- rep(1/nrow(dist_matrix), nrow(dist_matrix))
    w2 <- rep(1/ncol(dist_matrix), ncol(dist_matrix))
    a <- transport(w1, w2, costm = dist_matrix, method = "shortsimplex")
    cost <- 0 
    for (i in 1:nrow(a)){
      cost <- cost + dist_matrix[a$from[i], a$to[i]] * a$mass[i]
    }
    costs[city_comp] = cost
  }
  
  # Plot 
  X = barplot(costs, las = 2, xaxt = 'n', ylab = "1-Wasserstein distance", 
              main = names(All_results)[j])
  text(x = X, y = -0.02, srt = 40, adj= 1, xpd = TRUE,
       labels =  gsub('_', ' & ', names(costs)), cex = 0.7)
}


#================================================================
# Calculate the "1-Wasserstein" distance between population 1 and population 2
#================================================================

nrepeats = 100 # Number of bootstrap repeats

CIs = sapply(geo_dist_info$geo_order, function(city_comp){
  
  mle_CIs = All_results[[j]]
  if(class(mle_CIs$individual1) == 'factor'){mle_CIs$individual1 = as.character(mle_CIs$individual1)}
  if(class(mle_CIs$individual2) == 'factor'){mle_CIs$individual2 = as.character(mle_CIs$individual2)}
  
  # Add site comps
  mle_CIs$City1 = as.character(SNPData[mle_CIs$individual1, 'City'])
  mle_CIs$City2 = as.character(SNPData[mle_CIs$individual2, 'City'])
  mle_CIs$City12 = apply(mle_CIs[,c('City1','City2')], 1, function(x)paste(sort(x),collapse="_"))
  
  cities = strsplit(city_comp, split = '_')[[1]]
  inds = (mle_CIs$City1 %in% cities) & (mle_CIs$City2 %in% cities) # indices for cities

  #--------------------------------------------------------
  # Using vertboot - I don't understand exactly what vertboot is doing seems to only work
  # with integers: e.g. vertboot(m1 = matrix(rnorm(10,10),10,10), boot_rep = 10) 
  # Scale to map adj_matrix onto integer values for vertboot (see example above),
  # Scale = 10^8
  #
  # # Return adj matrix
  # adj_matrix_list = construct_adj_matrix(Result = mle_CIs[inds, ], dm = F)
  # no_c1 = length(adj_matrix_list$samples_c1) # sumarise sample count
  # no_c2 = length(adj_matrix_list$samples_c2) # sumarise sample count
  # 
  # # Bootstrap adj matrix
  # adj_matrix_boot = vertboot(m1 = adj_matrix_list[[1]]*Scale, boot_rep = nrepeats)
  # 
  # # Extract distance matrices
  # dis_matrix_boot = lapply(adj_matrix_boot, function(m)m[1:no_c1,no_c1+1:no_c2]/Scale)
  #--------------------------------------------------------
  
  #--------------------------------------------------------
  # Simply bootstrapping contents  
  dist_matrix = construct_adj_matrix(Result = mle_CIs[inds, ], dm = T)
  dist_matrix_boot = lapply(1:nrepeats, function(b){
    row_boot = sample(nrow(dist_matrix), nrow(dist_matrix), replace = T)
    col_boot = sample(ncol(dist_matrix), ncol(dist_matrix), replace = T)
    dm_boot = rbind(dist_matrix[row_boot,])
    dm_boot = cbind(dm_boot[row_boot,])
    return(dm_boot)
  })
  #--------------------------------------------------------
  
  # Need to boostrap the distance matrix somehow
  costs_boot = sapply(dist_matrix_boot, function(dist_matrix){
    w1 <- rep(1/nrow(dist_matrix), nrow(dist_matrix))
    w2 <- rep(1/ncol(dist_matrix), ncol(dist_matrix))
    a <- transport(w1, w2, costm = dist_matrix, method = "shortsimplex")
    cost <- 0
    for (i in 1:nrow(a)){
      cost <- cost + dist_matrix[a$from[i], a$to[i]] * a$mass[i]
    }
    return(cost)
  })
  
  # Return quantiles
  quantile(costs_boot, probs = c(0.025,0.975))
})

segments(x0 = X[,1], x1 = X[,1], y0 = CIs['2.5%',], y1 = CIs['97.5%',])





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





       