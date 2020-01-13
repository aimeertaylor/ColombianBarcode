library(transport) # transport package

# Function to calculate 1-W distance 
generate_1_w = function(dist_matrix){
  
  # Weights (could adapted to propogate uncertainity)
  w1 <- rep(1/nrow(dist_matrix), nrow(dist_matrix))
  w2 <- rep(1/ncol(dist_matrix), ncol(dist_matrix))
  
  # Caclulate 1-W
  a <- transport(w1, w2, costm = dist_matrix, method = "shortsimplex")
  cost <- 0 
  
  for (i in 1:nrow(a)){
    cost <- cost + dist_matrix[a$from[i], a$to[i]] * a$mass[i]
  } 
  return(cost)
}

# Function to calculate 1-W distance and bootstrapped 
generate_1_w_CI = function(dist_matrix, nrepeats){
  
  # Bootstrapping rows then columns (samples from citites 1 and 2) 
  dist_matrix_boot = lapply(1:nrepeats, function(b, dist_matrix){
    row_boot = sample(nrow(dist_matrix), nrow(dist_matrix), replace = T)
    col_boot = sample(ncol(dist_matrix), ncol(dist_matrix), replace = T)
    dm_boot = rbind(dist_matrix[row_boot,])
    dm_boot = cbind(dm_boot[,col_boot])
    return(dm_boot)
  }, dist_matrix)
  
  # Calculate 1-W  boostrapped
  costs_boot <- sapply(dist_matrix_boot, generate_1_w)
  
  # Return value and CIs
  to_return <- quantile(costs_boot, probs = c(0.025,0.975))
  
  return(to_return)
}


