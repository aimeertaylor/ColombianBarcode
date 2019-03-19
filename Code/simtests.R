# Similar to as.rtest from ade4 package, but 
# 1) allows for vectors as well as scalers
# 2) returns both one and two-tailed, where latter is the summation over left and right
mcarlo_rtest <- function(sim, obs){ 
  if(!is.matrix(sim)){sim <- matrix(sim, nrow = 1)}
  res <- list(sim = sim, obs = obs)
  res$nrep <- ncol(sim)
  res$pvalue_left <- (rowSums(sim <= pmin(obs, -obs)) + 1)/(res$nrep + 1) # one-sided
  res$pvalue_right <- (rowSums(sim >= pmax(obs, -obs)) + 1)/(res$nrep + 1) # one-sided
  res$pvalue <- res$pvalue_left + res$pvalue_right # two-sided
  class(res) <- "rtest"
  return(res)}

# Similar to as.rtest from ade4 package, but 
# 1) exact (assumes all perturbations were enumerated)
# 2) allows for vectors as well as scalers
# 3) returns both one and two-tailed, where latter is the summation over left and right
exact_rtest <- function(sim, obs){ 
  if(!is.matrix(sim)){sim <- matrix(sim, nrow = 1)}
  res <- list(sim = sim, obs = obs)
  res$nrep <- ncol(sim)
  res$pvalue_left <- rowSums(sim <= pmin(obs, -obs))/res$nrep # one-sided
  res$pvalue_right <- rowSums(sim >= pmax(obs, -obs))/res$nrep # one-sided
  res$pvalue <- res$pvalue_left + res$pvalue_right # two-sided
  class(res) <- "rtest"
  return(res)}


# Perturbation test for univariate linear regression
permute.univlm <- function(x, y, nrep = 1000){ 
  
  library(combinat) # for permn
  n <- factorial(length(y))
  obs <- lm(y ~ x)$coefficient['x']
  
  if(n < nrep){ # Exact p-value
    allperm <- permn(y)
    sim <- sapply(allperm, function(y){lm(y ~ x)$coefficient['x']})
    res <- exact_rtest(sim, obs)
  } else { # Monte Carlo p-value
    z <- vector('list', length = 1000)
    sim <- sapply(z, function(z){
      yrep <- y[sample(length(y))] # Permute
      lm(yrep ~ x)$coefficient['x']}) # Calculate test statistic
    res <- mcarlo_rtest(sim, obs)
  }
  
  result <- list(pvalue = res$pvalue, exact = n < nrep, n = n, obs = obs) 
  return(result)
}

