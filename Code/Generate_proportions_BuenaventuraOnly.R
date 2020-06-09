#######################################################
# This script generates proportions 
# given a single example threshold (e.g. 0.25)
# per year for Buenaventura 
#######################################################
rm(list = ls())

# Load and summarise raw data 
load('../RData/All_results.RData') 
load('../RData/geo_dist_info_cities.RData')
mle_CIs <- All_results$Unfiltered
nrep <- 100 # For bootstrap confidence intervals
r_threshold = 0.25
set.seed(1) # For reproducibility 

## Overall
# Restrict the Buenaventura only
Buenaventura_inds <- mle_CIs$City12 == 'Buenaventura_Buenaventura'
mle_CIs <- mle_CIs[Buenaventura_inds, ]
rhats <- mle_CIs$`r2.5.`
n_Buenaventura <- sum(Buenaventura_inds)

# Boostrap proportions 
prob_b <- sapply(1:nrep, function(b){ 
  mean(sample(rhats, size = n_Buenaventura, replace = TRUE) > r_threshold)})

# Observed proportion 
proportions_overall <- c('mean' =  mean(rhats > r_threshold), 
                          quantile(prob_b, probs = c(0.025, 0.975)))

## By year
# Add years 
mle_CIs$year1 = do.call(rbind, strsplit(mle_CIs$date1, split = "-"))[,1]
mle_CIs$year2 = do.call(rbind, strsplit(mle_CIs$date2, split = "-"))[,1]

# Restrict to same year only
SameYr_inds <- mle_CIs$year1 == mle_CIs$year2
mle_CIs <- mle_CIs[SameYr_inds, ]

# Unique years
years <- sort(unique(mle_CIs$year1))

# Stores (inc. those where intervals are partioned into site_comps and time bins)
proportions_year = array(0, dim = c(3, length(years)), dimnames = list(c('mean', '2.5%', '97.5%'), years))

# Calculate proportions with site
for(i in years){
  
  # Extract data
  ind <- mle_CIs$year1 == i
  if (unique(mle_CIs$year2[ind]) != i) stop ("Something wrong")
  n_yr <- sum(ind)
  rhats <- mle_CIs$`r2.5.`[ind]
  print(n_yr)
  
  # Boostrap proportions 
  prob_b <- sapply(1:nrep, function(b) mean(sample(rhats, size = n_yr, replace = TRUE) > r_threshold))
  
  # Observed proportion 
  proportions_year[,i] <- c('mean' =  mean(rhats > r_threshold), quantile(prob_b, probs = c(0.025, 0.975)))
}

proportions <- cbind(proportions_year, overall = proportions_overall)


# Plot (n.b. the proportion overall is not an average of the proportions per year 
# as overall inc. comparisons across years whereas per year do not)
x <- barplot(proportions["mean",], ylim = range(proportions), 
             ylab = expression('Fraction of highly-related'~italic('P. falciparum')~'sample pairs from Buenaventura'))
segments(x0 = x, x1 = x, y0 = proportions["2.5%",], y1 = proportions["97.5%",])







