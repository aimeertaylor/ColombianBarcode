##############################################################
# Script to plot igraph results 
# Note width < 1 does not show on pdfs so map r to transpancy
#
# To-do: reduce redundancy (e.g. A_related)
##############################################################

rm(list = ls())
library(igraph) # To make graph, construct components etc.
library(RColorBrewer) # For colours
load('../RData/All_results.RData') # Load All_results
load('../RData/SNPData.RData') # Load SNP data for cities
load('../RData/geo_dist_info.RData')
source('./igraph_functions.R')
eps = 0.01 # Below which LCI considered close to zero
cols = colorRampPalette(brewer.pal(12, "Paired")) # Function to create colours
FILTER = F # To re-generate filtered results
PDF = F # Set to TRUE to plot graphs to pdf

mle_CIs <- All_results$Unfiltered # Could just use All_results$Unfiltered


# 5 x 1 vector of city colours named by city
cities = unique(SNPData$City)
cols_cities = brewer.pal(length(cities), 'Spectral') 
names(cols_cities) = rev(cities) 

#===========================================================
# Using unfiltered results, generate CC adj matrix and graph
#===========================================================
# Adj. matrix unfiltered using 'r97.5.' (not rhat as in All_adj_matrix$Unfiltered)
A_high = construct_adj_matrix(mle_CIs, Entry = 'r97.5.')
A_high[A_high < 1-eps] = 0 # Edit s.t. only not stat diff from clonal have weight
G_high = graph_from_adjacency_matrix(A_high, mode='upper', diag=F, weighted=T) # Construct graph 
C_high = components(G_high) # Extract CC from graph
M_high = C_high$membership # Extract CC membership per sample (numeric name of each CC)

# Create a vector of all CCs with two or more samples and assign each a unique colour
CCs = cols(sum(C_high$csize > 1)) # First get unique colours 
names(CCs) = (1:C_high$no)[C_high$csize > 1] # Name CCs using numeric CC name given igraph

# Create a vector of vertex colours by CC
M_cols = sapply(M_high, function(x){ # Return white for singleton 
  ifelse(C_high$csize[x]==1,'#FFFFFF',CCs[as.character(x)])})


############################################################
# 2) In this section we plot the various graphs and components
############################################################
if(PDF){pdf('../Plots/All_CCs.pdf', height = 8, width = 8)}
par(mfrow = c(1,1), family = 'serif')

#++++++++++++++++++++++++++++++++++++++++++++
# Put dates and location of first observed in table? 
# To access dates: '\n', SNPData[V(Comp_G)$name, 'COLLECTION.DATE'])
# Extract vertices per edge: v_names = do.call(rbind, strsplit(attributes(E(G_high))$vnames, split = "\\|")) 
# James said that edges hard to see when printed
#===========================================================
# Plot the relatedness between components
# with edge only if statistically distinguishable from zero
# and edge colour proportion to remaining rhats 
#===========================================================

# To construct a CC graph whose edges are weighted according to the 
# earliest sample per clonal component...
# Remove all but but the first sample per CC regardless of city
# using a dummy vector that says all samples come from one city
R_comp = rm_highly_related_within(Result = mle_CIs, Edge = F, 
                                  Cities = sapply(row.names(SNPData), function(x)"DummyCity"))
# Create adjacency matrix over clonal components
A_comp = construct_adj_matrix(R_comp, Entry = 'rhat') 
# Set edges that are statistically indistinguishable from zero to zero
A_comp[construct_adj_matrix(R_comp, Entry = 'r2.5.') < eps] = 0
# Create adjacency matrix over clonal components
A_comp = construct_adj_matrix(R_comp, Entry = 'rhat_av')

# Create graph
G_comp = graph_from_adjacency_matrix(A_comp, mode='upper', diag=F, weighted=T)

# Colour vertices by CC
V(G_comp)$color = M_cols[V(G_comp)$name]

# Based on colour, extract subgraph of those with two or more members
Comp_G = induced_subgraph(G_comp, vids = which(V(G_comp)$color!="#FFFFFF"))

# Annotate the subgraph
E(Comp_G)$width <- E(Comp_G)$weight
E(Comp_G)$colour <- 'black'
V(Comp_G)$date <- SNPData[V(Comp_G)$name, 'COLLECTION.DATE']
V(Comp_G)$site <- SNPData[V(Comp_G)$name, 'City']
V(Comp_G)$cc <- as.character(M_high[V(Comp_G)$name])

# Create CC names for CCs with two or more only in order of date of earliest sample per CC
C_names = paste0('CC',1:length(CCs)) # Create CC names based on colour
ordered_date_index <- sort.int(V(Comp_G)$date, index.return = T)$ix
sid_ordered_date = V(Comp_G)$name[ordered_date_index] # reorder sample ID by date 
names(C_names) = as.character(M_high[sid_ordered_date]) # Make sure the CC names are same order as memberships

# Check order and print mapping
cbind(C_names[as.character(M_high[V(Comp_G)$name[ordered_date_index]])], 
      as.character(SNPData[V(Comp_G)$name[ordered_date_index], 'COLLECTION.DATE']))

# CCs 12 and 13 have earliest parasite samples detected on the same date: "2002-04-03"
duplicate <- C_names[as.character(M_high[V(Comp_G)$name[duplicated(V(Comp_G)$date)]])]
writeLines(sprintf('Note that %s has a duplicate date', duplicate))

# -----------------------------------------
# Aside: related components among CCs 
writeLines(sprintf("Among the clonal components, there are %s related components size %s", 
                   components(Comp_G)$no, 
                   paste(components(Comp_G)$csize, collapse = " and ")))

# Find and inspect parasite members of the anamolous CC
anomaly_CC_sID1 = names(which(components(Comp_G)$membership != 1)) # The sample ID retained in Comp_G
anomaly_CC_sIDs = which(M_high == M_high[anomaly_CC_sID1]) # All sample IDs
C_names[as.character(M_high[anomaly_CC_sID1])] # Check name

# Extract and print to screen their data (check with Diego )
anomaly_CC_SNPData = SNPData[anomaly_CC_sIDs,]
writeLines(sprintf("The anomalous CC (unrelated to the other CCs) contains %s samples:", 
                   length(anomaly_CC_sIDs)))
print(anomaly_CC_SNPData[,1:5])
save(anomaly_CC_SNPData, file = '../RData/anomaly_CC_SNPData.RData')
# -----------------------------------------

# scale weights to zero one for maximal visualisation
weights_rescaled = (E(Comp_G)$weight - min(E(Comp_G)$weight)) / (max(E(Comp_G)$weight) - min(E(Comp_G)$weight))
range(E(Comp_G)$weight)

# Plot without pie
set.seed(1)
plot(Comp_G, 
     layout = layout_with_fr, 
     vertex.size = table(M_cols)[V(Comp_G)$color]+5, # Number of members per CC
     vertex.label = C_names[V(Comp_G)$cc], # Label by CC name 
     vertex.label.cex = 0.5, # Label size
     vertex.label.color = 'black', # Label colour
     vertex.color = cols_cities[V(Comp_G)$site], # Colour by city of first sample per CC
     vertex.frame.color = 'white', # Add white border for visual separtion
     edge.color = sapply(weights_rescaled, adjustcolor, col = 'black') # Use transparencey not weight
)


#==================================================================================
# Construct a CC pie-graph whose edges are weighted according to average relatedness 
#==================================================================================.
A_hat = construct_adj_matrix(mle_CIs, Entry = 'rhat')
A_low = construct_adj_matrix(mle_CIs, Entry = 'r2.5.')
A_hat[A_low < eps] = 0 # Edit s.t. only not stat diff from clonal have weight

# For each pair of CCs with 2 or more samples, extract samples and calculate average relatedness
CCpairs <- gtools::combinations(n = length(CCs), r = 2, v = names(CCs))
colnames(CCpairs) = c('individual1','individual2')
R_comp <- data.frame(CCpairs, 'rhat_av' = NA, stringsAsFactors = F)
for(i in 1:nrow(R_comp)){

  # Extract samples per CC
  samples_CC1 <- names(M_high)[as.character(M_high) == R_comp$individual1[i]]
  samples_CC2 <- names(M_high)[as.character(M_high) == R_comp$individual2[i]]
  
  # Construct sample pairs that bridge CCs
  sample_pairs_between_CCs <- expand.grid(samples_CC1, samples_CC2)
  
  # Extract rhats 
  rhats = apply(sample_pairs_between_CCs, 1, function(sids){
    rhat = A_hat[sids[1], sids[2]]
    # Thus if rhat is na try swapping rows and columns:
    if(is.na(rhat)){rhat <- A_hat[sids[2], sids[1]]}
    return(rhat)
  })
  
  R_comp[i, 'rhat_av'] = mean(rhats) # Average rhats
}
# Create adjacency matrix over clonal components
Comp_A = construct_adj_matrix(R_comp, Entry = 'rhat_av')
# Create graph
Comp_G = graph_from_adjacency_matrix(Comp_A, mode='upper', diag=F, weighted=T)

# For each CC with 2 or more samples, extract number of samples per city (order consistently)
sample_count_per_city_per_CC <- lapply(names(CCs), function(CC){
  
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
names(sample_count_per_city_per_CC) = names(CCs) 

# Scale edge weights to zero one 
weights_rescaled = (E(Comp_G)$weight - min(E(Comp_G)$weight)) / (max(E(Comp_G)$weight) - min(E(Comp_G)$weight))
E(Comp_G)$width = weights_rescaled *2

# Plot as pie (Tim's suggestion)
set.seed(1)
plot(Comp_G, #layout = layout_with_fr, 
     vertex.shape = "pie", 
     vertex.pie = sample_count_per_city_per_CC[V(Comp_G)$name],
     vertex.pie.color = list(cols_cities[cities]),
     vertex.size = sapply(sample_count_per_city_per_CC[V(Comp_G)$name], sum) + 5, 
     vertex.label = C_names[V(Comp_G)$name], 
     vertex.label.cex = 0.5, 
     vertex.label.color = 'black',
     vertex.frame.color = NA,
     edge.color = 'black')


# Legend of cities
legend('left', pch = 16, bty = 'n', cex = 0.7, pt.cex = 1.5, col = cols_cities, 
       legend = names(cols_cities), inset = -0.1)

if(PDF){dev.off()}





#===========================================================
# For each site comparison vizualise effect of filter with
# Dates and edges proportional to rhats 
#===========================================================
if(PDF){pdf('../Plots/Graphs_timespace.pdf', height = 5, width = 10)}
par(mfrow = c(1,1), mar = c(,4,3,0.1), family = 'serif', bg = 'white')
filter_name = c("Unfiltered" = "Unfiltered", "Filter by vertex" = "Filter CCs")
SHAPES = c('circle', 'square') # Different shapes distinguish different sites 
for(l in c("Unfiltered","Filter by vertex")){
  
  # Extract adjacency 
  X = All_results[[l]]
  A_related_lci = construct_adj_matrix(Result = X, Entry = 'r2.5.')
  A_related = construct_adj_matrix(Result = X, Entry = 'rhat')
  A_related[A_related_lci < eps] = 0 # Edit s.t. only not stat diff from related have weight
  G_related = graph_from_adjacency_matrix(A_related, mode='upper', diag=F, weighted=T)
  # Add meta data
  V(G_related)$site = SNPData[V(G_related)$name, 'City']
  V(G_related)$date = SNPData[V(G_related)$name, 'COLLECTION.DATE']
  
  
  for(i in 1:nrow(geo_dist_info$pairwise_site_distance)){
    
    # Extract cities
    s1 = as.character(geo_dist_info$pairwise_site_distance[i,'X1'])
    s2 = as.character(geo_dist_info$pairwise_site_distance[i,'X2'])
    
    # Filter graph by city
    G = extract_city_graph(x = G_related, city1 = s1, city2 = s2, 'black')
    
    # Calculate years
    colldates = SNPData[V(G)$name, "COLLECTION.DATE"]
    min_dt = min(colldates)
    max_dt = max(colldates)
    years1stjan = seq(as.Date(round.POSIXt(min_dt, units = "years")), max_dt, by = 'year') 
    yr01 = as.numeric(years1stjan-min_dt)/as.numeric(max_dt - min_dt)
    
    # Plot graph
    plot.igraph(G, layout = attributes(G)$layout, vertex.size = 0, edge.color = NA,
                vertex.label = NA, asp = 0) # For layout
    abline(v = -1 + yr01 * 2, lty = 'dotted', col = 'lightgray') # Add grid line
    plot.igraph(G, layout = attributes(G)$layout, vertex.shape = SHAPES[as.numeric(V(G)$site == s1)+1], 
                vertex.size = 3, vertex.label = NA, add = T) # For layout
    
    # Overlay coloured only
    E(G)$color[grepl("#000000", E(G)$color)] = NA
    plot.igraph(G, layout = attributes(G)$layout, vertex.shape = SHAPES[as.numeric(V(G)$site == s1)+1], 
                vertex.size = 3, vertex.label = NA, add = T)
    
    # Annotate
    mtext(side = 3, line = 2, adj = 0, sprintf('%s %s %s', filter_name[l], s1, s2), cex = 1)
    # legend('left', pch = 23, bty = 'n', inset = -0.07, 
    #        pt.bg = CCs[CCs %in% unique(V(G)$color[V(G)$color!="#FFFFFF"])], 
    #        cex = 0.7, pt.cex = 1,
    #        legend = C_names[as.character(names(CCs[CCs %in% unique(V(G)$color[V(G)$color!="#FFFFFF"])]))])
    # legend('topleft', pch = c(0,1), pt.bg = 'white', legend = c(s1,s2), 
    #        bty = 'n', inset = -0.07)
    axis(side = 1, labels = years1stjan, las = 2, cex.axis = 0.7, at = yr01 * 2 - 1,
         line = -1.2, tick = F) # Years
  }
}

if(PDF){dev.off()}



#===============================================
# Plot without dates and with edges if
# statistically distinguishable from zero with transparency
# proportional to rhat
#===============================================
if(PDF){pdf('../Plots/Graphs.pdf', height = 5, width = 10)}
par(mfrow = c(1,1), mar = c(0,4,0,0), family = 'serif', bg = 'white', pty = 'm')
# Create matrixA_related = construct_adj_matrix(Result = mle_CIs, Entry = 'r2.5.')
A_related = construct_adj_matrix(Result = mle_CIs, Entry = 'rhat')
A_related_lci = construct_adj_matrix(Result = mle_CIs, Entry = 'r2.5.')
A_related[A_related_lci < eps] = 0 # Edit s.t. only not stat diff from related have weight

# Create graphs 
G_related = graph_from_adjacency_matrix(A_related, mode='upper', diag=F, weighted=T)
V(G_related)$site = SNPData[V(G_related)$name, 'City']

# ---------------------------------------------------------------------
# Aside: note that there are two related components across the entire graph
# The second consists of the DD2 clonal component plus a loosely related 
# sample, "sid95". This is the outlier in the Buenaventura Tumaco plot 
# (check by doing components(G) when l = "Unfiltered" and i = 6 in the above for loop)
components(G_related)
sids_related_to_sid95 = which(A_related[rownames(A_related) == "sid95",] > 0)
A_related["sid95",sids_related_to_sid95]
SNPData[c(names(sids_related_to_sid95),"sid95"),1:5]
# ---------------------------------------------------------------------

set.seed(1)
SHAPES = c('circle', 'square') # Different shapes distinguish different sites 

for(i in 1:nrow(geo_dist_info$pairwise_site_distance)){
  
  # Extract cities
  s1 = as.character(geo_dist_info$pairwise_site_distance[i,'X1'])
  s2 = as.character(geo_dist_info$pairwise_site_distance[i,'X2'])
  
  # Filter graph by city
  G = extract_city_graph(x = G_related, city1 = s1, city2 = s2, 'black')
  
  # # Remove the sid95 outlier to make better use of spacef for more clear visualisation
  # if (s1 %in% c("Buenaventura", "Tumaco") & s2 %in% c("Buenaventura", "Tumaco")) {
  #   G = delete_vertices(G, "sid95")
  # }
  
  # Plot with legend
  set.seed(1)
  # To do the igraph equivalent of pty = "m" (looks ugly, prefer square)
  #l = layout_(G, with_fr()) # set layout using Fruchterman-Reingold
  #l = norm_coords(l, par()$usr[1], par()$usr[2], par()$usr[3], par()$usr[4]) # normalise to devise coordinates not square
  plot.igraph(G, vertex.size = 3, vertex.label = NA, layout = layout_with_fr, 
              vertex.shape = SHAPES[as.numeric(V(G)$site == s1)+1])
  legend('left', pch = 23, bty = 'n', inset = 0.1, y.intersp = 0.7,
         pt.bg = CCs[CCs %in% unique(V(G)$color[V(G)$color!="#FFFFFF"])], cex = 0.8, pt.cex = 1, 
         legend = C_names[as.character(names(CCs[CCs %in% unique(V(G)$color[V(G)$color!="#FFFFFF"])]))])
  legend('topleft', pch = c(0,1), pt.bg = 'white', legend = c(s1,s2), bty = 'n', inset = 0.7)
}

if(PDF){dev.off()}


#===========================================================
# Looking specifically at clonal components found across 
# Buenaventura and Tumaco (some of also found in Guapi)
# Edge transparency proportion to rhat if statistically
# distinguishable from zero (since based on Comp_G)
#===========================================================
par(bg = 'white')
G_B = induced_subgraph(G_related, vids = which(grepl("Buenaventura", V(G_related)$site)))
G_T = induced_subgraph(G_related, vids = which(grepl("Tumaco", V(G_related)$site)))

# Some NAs due to only components with two or more samples having names
BT_CCs = intersect(C_names[as.character(M_high[V(G_B)$name])], C_names[as.character(M_high[V(G_T)$name])])
BT_CCs = BT_CCs[!is.na(BT_CCs)]
unique_M_BT = as.numeric(names(C_names[C_names %in% BT_CCs]))

# Extract all those with said membership (includes some from guapi)

G_clonal_comp_BT = induced_subgraph(All_G$Unfiltered, vids = which(M_high %in% unique_M_BT)) 

# Plot all those with said membership
set.seed(2)
plot.igraph(G_clonal_comp_BT, vertex.size = 3, vertex.label = NA, 
            vertex.color = M_cols[V(G_clonal_comp_BT)$name], 
            edge.color = sapply(E(G_clonal_comp_BT)$weight, adjustcolor, col = 'black')) # Make transparency a function of rhat)

# Extract edge vertex names as matrix
EV_names = do.call(rbind, strsplit(attributes(E(G_clonal_comp_BT))$vname, split = "\\|"))

# Extract average relatedness between clonal components (includes some guapi)
Mean_sd_matrix = array(NA, dim = c(length(unique_M_BT), length(unique_M_BT), 6), 
                       dimnames = list(C_names[as.character(unique_M_BT)], C_names[as.character(unique_M_BT)], 
                                       c('r mean','r 2.5%','r 97.5%',
                                         'k mean','k 2.5%','k 97.5%')))
for(j in 1:(length(unique_M_BT)-1)){
  for(i in (j+1):length(unique_M_BT)){
    members_i = names(which(M_high == unique_M_BT[i])) 
    members_j = names(which(M_high == unique_M_BT[j])) 
    inds = (EV_names[,1] %in% members_i & EV_names[,2] %in% members_j) | EV_names[,1] %in% members_j & EV_names[,2] %in% members_i 
    Mean_sd_matrix[i,j,'r mean'] = mean(E(G_clonal_comp_BT)$weight[inds])
    
    # Extract the mean lower confidence interval 
    inds_mle = mle_CIs$sample_comp %in% apply(EV_names[inds,], 1, function(x)paste(x[1],x[2],sep='_')) | mle_CIs$sample_comp %in% apply(EV_names[inds,], 1, function(x)paste(x[2],x[1],sep='_'))
    if (Mean_sd_matrix[i,j,'r mean'] == mean(mle_CIs$rhat[inds_mle])){ # check
      Mean_sd_matrix[i,j,'r 2.5%'] = max(mle_CIs$`r2.5.`[inds_mle]) # To see highest lower bound
      Mean_sd_matrix[i,j,'r 97.5%'] = max(mle_CIs$`r97.5.`[inds_mle]) 
      Mean_sd_matrix[i,j,'k mean'] = max(mle_CIs$khat[inds_mle])
      Mean_sd_matrix[i,j,'k 2.5%'] = max(mle_CIs$`k2.5.`[inds_mle]) # To see highest lower bound
      Mean_sd_matrix[i,j,'k 97.5%'] = max(mle_CIs$`k97.5.`[inds_mle]) 
    } 
  }
}

# Order matrices by CC name before printing
order_by_CC_name = sort.int(as.numeric(gsub('CC', '',C_names[as.character(unique_M_BT)])), index.return = T)$ix

# Print to screen
round(Mean_sd_matrix[order_by_CC_name,order_by_CC_name,], 3)

# Print for overleaf
require(kableExtra)
kable(round(Mean_sd_matrix[order_by_CC_name,order_by_CC_name,'r mean'],3), format = 'latex') 
round(Mean_sd_matrix[order_by_CC_name,order_by_CC_name, 'r 2.5%'], 3)

# Aside: relatedness between CC20, CC15 and CC14: valid example of recombination? Likely via some intermediates
Ms_CC20_15_14 <- names(C_names[C_names %in% c("CC20", "CC15", "CC14")])
Ms_Comp_G <- as.character(M_high[V(Comp_G)$name])
earliest_sids_CC20_15_14 <- V(Comp_G)$name[Ms_Comp_G %in% Ms_CC20_15_14]
C_names[as.character(M_high[earliest_sids_CC20_15_14])] # check
G_earliest_sids_CC20_15_14 <- induced_subgraph(Comp_G, vids = earliest_sids_CC20_15_14)
plot(G_earliest_sids_CC20_15_14, 
     vertex.label = C_names[as.character(M_high[V(G_earliest_sids_CC20_15_14)$name])])
E(G_earliest_sids_CC20_15_14)$weight
