##############################################################
# Script to plot igraph results 
#' Remove low quality samples if CCs very sensitive or 
#' use cliques instead of components. Also, consider how to fix 
#' CCs as per Taylor et al. 2020 and add to
#' Since rm_highly_related_within uses SNPData, need to rename 
#' metadata "SNPData" (or change rm_highly_related_within) and 
#' make sure it is in the same format
##############################################################
rm(list = ls())
library(igraph) # To make graph, construct components etc.
library(RColorBrewer) # For colours
library(kableExtra)
source('../igraph_functions.R')

freqs_used <- "Taylor2020"
load('../../RData/SNPData.RData') # Needed for rm_highly_related_within
load(sprintf('../../RData/mles_CIs_meta_extended_freqs%s.RData', freqs_used)) # Load All_results

eps <- 0.01 # Below which LCI considered close to zero
a = 0.05; b = 1; c = 2 # Scaling parameters for edge width and transparancy 
cities <- unique(c(mle_CIs$City1, mle_CIs$City2))
cols_cities <- brewer.pal(length(cities), 'Spectral') 
names(cols_cities) <- as.character(cities)
PDF <- F # Set to TRUE to plot graphs to pdf

#===========================================================
# Get clonal component (CC) membership for all samples
# Remove low quality samples to explore effect on CC
#===========================================================

#-----------------------------------------------------------
# Create graph of all samples to extract CCs
A_high <- construct_adj_matrix(mle_CIs, Entry = 'r97.5.') # Adj. matrix unfiltered using 'r97.5.' 
A_high[A_high < 1-eps] <- 0 # Edit s.t. only clonal have weight
G_high <- graph_from_adjacency_matrix(A_high, mode='upper', diag=F, weighted=T) # Construct graph 

# Last attempt using components resulted in one component with 506 samples
# Extract CCs from graph using cliques instead of components
C_high <- components(G_high) # Extract CCs from graph using cliques instead of components
M_high <- C_high$membership # Extract CC membership per sample (numeric name of each CC)
#----------------------------------------------------------

#----------------------------------------------------------
# Create a vector of all CCs with 2+ samples 
CCs <- as.character(which(C_high$csize > 1)) 


#----------------------------------------------------------
# Create CC character names and order by earliest detected sample per CC
CC_chr_names <- paste0('CC',1:length(CCs)) 
# Map CC character names to dates of earliest detected sample per CC
# Remove all duplicate CCs regardless of city using a dummy city name
Rhat_filtCC = rm_highly_related_within(Result = mle_CIs, Edge = F, 
                                       Cities = sapply(unique(c(mle_CIs$individual1, mle_CIs$individual2)),
                                                       function(x)"DummyCity"))

# Extract remaining sample ids
sids_noCCdupl = unique(c(as.character(Rhat_filtCC$individual1), 
                         as.character(Rhat_filtCC$individual2))) 

# Keep only those in a CC with one or more sample
sids_1stperCC = sids_noCCdupl[M_high[sids_noCCdupl] %in% CCs]

# Extract dates for filtered sample ids 
sample_dates <- SNPData[sids_1stperCC, 'COLLECTION.DATE'] # extract dates of 
order_of_sample_dates <- sort.int(sample_dates, index.return = T)$ix # extract order 
sids_1stperCC_ordered <- sids_1stperCC[order_of_sample_dates] # reorder sample ID by date 
names(CC_chr_names) = as.character(M_high[sids_1stperCC_ordered]) # Ensure CC names are ordered as memberships
#--------------------------------------------------------

#--------------------------------------------------------
# Extract Table of site and date of earliest sample per CC 
First_sample_per_CC_table = cbind(CC = CC_chr_names, 
                                  Earliest_sample = sids_1stperCC_ordered, 
                                  Date_earliest_sample = as.character(SNPData[sids_1stperCC_ordered, 'COLLECTION.DATE']), 
                                  Site_earliest_sample = as.character(SNPData[sids_1stperCC_ordered, 'City']))
kableExtra::kable(First_sample_per_CC_table, format = "latex")
#--------------------------------------------------------

#--------------------------------------------------------
# Adjacency over all related with rhat statistically distinguishable from zero
# (no filter)
A_related = construct_adj_matrix(All_results$Unfiltered, Entry = 'rhat')
A_low = construct_adj_matrix(All_results$Unfiltered, Entry = 'r2.5.')
A_related[A_low < eps] = 0 # Edit s.t. only not stat diff from clonal have weight
#--------------------------------------------------------


#==================================================================================
# CC pie-graph whose edges are weighted according to average relatedness 
# averaged of those that are statistically distinguishable from zero
#==================================================================================.
if(PDF){pdf(sprintf('../../Plots/All_CCs_WGStoBarcode.pdf'), height = 8, width = 8)}
par(mfrow = c(1,1), family = 'serif')

for(Remove_singletons in c(TRUE,FALSE)){ # If true, remove CCs with only one parasite sample 
  
  # For each pair of CCs with 2 or more samples, extract samples and calculate average relatedness
  CCs_inc <- if(Remove_singletons){CCs}else{as.character(1:C_high$no)}
  CCpairs <- gtools::combinations(n = length(CCs_inc), r = 2, v = CCs_inc) # Get CC pairs
  colnames(CCpairs) = c('individual1','individual2')
  
  # Average over relatedness estimates that are statistically distinguishable from zero
  R_comp <- data.frame(CCpairs, 'rhat_av' = NA, stringsAsFactors = F) 
  for(i in 1:nrow(R_comp)){
    
    # Extract samples per CC
    samples_CC1 <- names(M_high)[as.character(M_high) == R_comp$individual1[i]]
    samples_CC2 <- names(M_high)[as.character(M_high) == R_comp$individual2[i]]
    
    # Construct sample pairs that bridge CCs
    sample_pairs_between_CCs <- expand.grid(samples_CC1, samples_CC2)
    
    # Extract rhats 
    rhats = apply(sample_pairs_between_CCs, 1, function(sids){
      rhat = A_related[sids[1], sids[2]]
      # Thus if rhat is na try swapping rows and columns:
      if(is.na(rhat)){rhat <- A_related[sids[2], sids[1]]}
      return(rhat)
    })
    
    # Average rhats
    R_comp[i, 'rhat_av'] = mean(rhats) 
  }
  
  # Create adjacency matrix over clonal components
  Comp_A = construct_adj_matrix(R_comp, Entry = 'rhat_av')
  # Create graph
  Comp_G = graph_from_adjacency_matrix(Comp_A, mode='upper', diag=F, weighted=T)
  
  # Aside: related components among CCs 
  if(Remove_singletons){
    writeLines(sprintf("Among the clonal components, there are %s related components size %s", 
                       components(Comp_G)$no, 
                       paste(components(Comp_G)$csize, collapse = " and ")))
    
    # Find and inspect parasite members of the anamolous CC
    anomaly_CC_sID1 = names(which(components(Comp_G)$membership != 1)) # The sample ID retained in Comp_G
    anomaly_CC_sIDs = names(which(M_high == anomaly_CC_sID1)) # All sample IDs
    
    
    # Extract and print to screen their data (check with Diego )
    anomaly_CC_SNPData = SNPData[anomaly_CC_sIDs,]
    writeLines(sprintf("The anomalous CC, %s, contains %s samples:", 
                       CC_chr_names[anomaly_CC_sID1], 
                       length(anomaly_CC_sIDs)))
    print(anomaly_CC_SNPData[,1:5])
  }
  
  # For each CC with 2 or more samples, extract number of samples per city (order consistently)
  sample_count_per_city_per_CC <- lapply(CCs_inc, function(CC){
    
    # Initialise empty vector of city counts
    city_count_per_CC_inc_zeros = array(0, dim = length(cities), dimnames = list(cities))
    
    # Extract samples per CC 
    samples_per_CC <- names(M_high)[as.character(M_high) == CC]
    
    # Extract number of samples per CC per city
    city_count_per_CC_exc_zeros <- table(SNPData[samples_per_CC,'City'])
    
    # Populate initial empty vector of city counts
    city_count_per_CC_inc_zeros[names(city_count_per_CC_exc_zeros)] <- city_count_per_CC_exc_zeros
    
    return(city_count_per_CC_inc_zeros)
  })
  
  # Name by CC numeric name
  names(sample_count_per_city_per_CC) = CCs_inc
  
  # Print non-zero edge weights before rescaling
  writeLines(sprintf("Relatedness values of edges range from %s to %s", 
                     round(min(E(Comp_G)$weight),3), round(max(E(Comp_G)$weight),3)))
  
  # Scale edge widths and transparancy to range a, b, c
  weights_rescaled = a + (b-a) * (E(Comp_G)$weight - min(E(Comp_G)$weight)) / (max(E(Comp_G)$weight) - min(E(Comp_G)$weight))
  vertex_size = sapply(sample_count_per_city_per_CC[V(Comp_G)$name], sum)
  vertex_size[vertex_size > 1] <- vertex_size[vertex_size > 1] + 5
  vertex_size[vertex_size == 1] <- vertex_size[vertex_size == 1] + 2
  
  # Plot as pie (Tim's suggestion)
  set.seed(1)
  plot(Comp_G, 
       layout = layout_with_fr, 
       vertex.shape = "pie", 
       vertex.pie = sample_count_per_city_per_CC[V(Comp_G)$name],
       vertex.pie.color = list(cols_cities[cities]),
       vertex.size = vertex_size, 
       vertex.label = CC_chr_names[V(Comp_G)$name], 
       vertex.label.cex = 0.35, 
       vertex.label.color = 'black',
       vertex.frame.color = NA,
       edge.width = weights_rescaled * c, 
       edge.color = sapply(weights_rescaled, function(x)adjustcolor('black', alpha.f = x)))
  
  # Legend of cities
  legend('left', pch = 16, bty = 'n', cex = 0.7, pt.cex = 1.5, col = cols_cities, 
         legend = names(cols_cities), inset = 0.1)
}
if(PDF){dev.off()}





