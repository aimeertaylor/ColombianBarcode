##############################################################
# In this script we generate transport ("1-Wasserstein") distances 
# between parasite populations sampled in five different time 
# intervals (five is arbitrary - chosen to match the city count). 
# 
# The transport approach is directly applicable to our between-city
# spatial analysis: it is trivial to construct a fully populated 
# adjacency matrix of comparisons across samples that fall into two 
# discrete populations (e.g. samples collected in Tumaco versus 
# Buenaventura). 
# 
# The transport approach is not directly applicable to our temporal 
# analysis: it is not possible to construct a fully populated and  
# representative adjacency matrix with distinct rows and columns from
# sample pairs that are labelled relatively (e.g. sample pairs collected one 
# week apart) ---- For example, if we were to randomly split all samples 
# collected one week apart into two distinct sets and thus construct a 
# adjacency matrix with distinct rows and columns, only a subset of samples 
# comparisons within the adjacency matrix would correspond to one week's 
# separation. The rest could either be given a NA value, resulting in a 
# sparse adjacency matrix that cannot be evaluated, or their relatedness value, 
# resulting in a non-representative matrix whose one-week signal is buried 
# in noise. Since almost all samples have at least one comparison of any 
# given time separation, noisy adjacency matrices are likely no more informative 
# that any random split of the entire data set. --- We thus create 
# pseudo-discrete time populations by splitting the time-line into intervals 
# and perform this analysis as a sanity check only (it isn't worthy of 
# publication). Note that within time-interval results are random since 
# they are based on a random split of the time-interval. 
##############################################################
rm(list = ls())
load('../RData/All_results_true.RData')
source('./transport_functions.R')
nrepeats = 100 # Number of bootstrap repeats 
set.seed(10)

# Unfiltered
mle_CIs = All_results$Unfiltered

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
# Calculate the "1-Wasserstein" distance over same versus different time intervals
# (not possible to otherwise construct a adjacency matrix)
#=============================================================================
min_date = min(mle_CIs$date1, mle_CIs$date2)
mle_CIs$week1 = difftime(time2 = min_date, time1 = mle_CIs$date1, units = "weeks")
mle_CIs$week2 = difftime(time2 = min_date, time1 = mle_CIs$date2, units = "weeks")
indivs = unique(c(mle_CIs$individual1, mle_CIs$individual2))

# Extract no. of weeks for n = 325 unique samples 
indiv_weeks = unique(mle_CIs[,c('individual1', 'week1')])
weeks = c(indiv_weeks$week1, sapply(setdiff(indivs, indiv_weeks$individual1), function(indiv){
  ind = which(mle_CIs$individual2 == indiv)[1]
  as.numeric(mle_CIs$week2[ind])
}))

# Using length.out = 6, create n = 5 intervals that are 
# 1) evenly populated using quatile() 
# 2) evenly spaced using seq() - does not work for n = 5 so use n = 4
week_bin_breaks = quantile(weeks, probs = seq(0,1,length.out = 6))
# week_bin_breaks = seq(min(weeks), max(weeks), length.out = 5)

# Allocate bins (.bincode is fast for assigning numeric bins)
mle_CIs$week_bin1 <- .bincode(as.numeric(mle_CIs$week1), 
                              breaks = week_bin_breaks, include.lowest = T)
mle_CIs$week_bin2 <- .bincode(as.numeric(mle_CIs$week2), 
                              breaks = week_bin_breaks, include.lowest = T)

# Extract bins, combinations thereof and order
week_bins = sort(unique(c(mle_CIs$week_bin1, mle_CIs$week_bin2)))
week_coms = gtools::combinations(n = length(week_bins), r = 2, v = week_bins, 
                                 repeats.allowed = TRUE)
week_coms_ordered = sort.int(abs(week_coms[,1] - week_coms[,2]), index.return = TRUE)

# Calculate 1-W distance for different combinations of time intervals
W_results_t = sapply(week_coms_ordered$ix, function(i){
  
  week_com = week_coms[i,]
  inds = (mle_CIs$week_bin1 %in% week_com) & (mle_CIs$week_bin2 %in% week_com) 
  if (length(inds))
  mle_weeks = mle_CIs[inds,] # Extract relevent mles
  
  if (week_com[1] == week_com[2]){
    
    # Split the sample into two, uniformly at random
    samples = unique(c(mle_weeks$individual1, mle_weeks$individual2))
    samples_c1 = sample(samples, size = floor(length(samples)/2), replace = F)
    samples_c2 = setdiff(samples, samples_c1)
    
  } else {
    
    # Samples per different city
    samples_c1 = unique(c(mle_weeks$individual1[mle_weeks$week_bin1 %in% week_com[1]], 
                          mle_weeks$individual2[mle_weeks$week_bin2 %in% week_com[1]]))
    samples_c2 = unique(c(mle_weeks$individual1[mle_weeks$week_bin1 %in% week_com[2]], 
                          mle_weeks$individual2[mle_weeks$week_bin2 %in% week_com[2]]))
  }
  
  # Extract distance matrix
  dist_matrix = adj_matrix[samples_c1, samples_c2]
  
  # Caclulate 1-W and CIs
  cost <- generate_1_w(dist_matrix)
  cost_CI <- generate_1_w_CI(dist_matrix, nrepeats)
  return(c(W = cost, cost_CI))
})

save(W_results_t, file = '../RData/W_results_t.RData')

# Quick plot 
par(mfrow = c(1,1))
X = barplot(W_results_t[1,], ylim = c(0, max(W_results_t)), 
            names.arg = sapply(week_coms_ordered$ix, function(i){
              paste(week_coms[i,], collapse = ' ')
            }), las = 2)
segments(x0 = X, x1 = X, y0 = W_results_t[2,], y1 = W_results_t[3,])

# Plot of time intervals
weeks_ordered = sort(weeks) 
plot(weeks_ordered, pch = 20, ylab = 'weeks', xlab = 'sample index', 
     col = .bincode(weeks_ordered, breaks = week_bin_breaks, include.lowest = T))

round(week_bin_breaks)

# Unsuprising that comparisons with one are more costly given
# the range of weeks interval one spans






