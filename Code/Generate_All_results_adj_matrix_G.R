rm(list = ls())
load('../RData/mles_frequencies_true.RData')
load('../RData/SNPData.RData') # Load SNP data for cities
source('./igraph_functions.R') # For rm_highly_related_within and construct_adj_matrix
eps = 0.01 # Below which LCI considered close to zero (needed by rm_highly_related_within)

# n x 1 vector of cities named by sample ID
Cities = SNPData$City; names(Cities) = row.names(SNPData) 

#===========================================================
# Filter and save results
#   - remove edges (almost all samples remain)
#   - remove vertices (removes all samples per CC per city except one)
#===========================================================
All_results = lapply(c(F,T), rm_highly_related_within, Result = mle_CIs, Cities = Cities)
All_results[[3]] = mle_CIs # Add unfiltered
names(All_results) = c('Filter by vertex', 'Filter by edge', 'Unfiltered')
save(All_results, file = '../RData/All_results.RData')

# Extract summaries
sapply(All_results, function(x){c('Edge count' = nrow(x), 
                                  'Vertex count' = length(unique(c(x$individual1,x$individual2))))})




