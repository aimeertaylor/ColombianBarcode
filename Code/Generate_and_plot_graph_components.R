##############################################################
# Script to 
#
# 1) filter results 
#   - remove edges and see if trend remains (almost all samples remain)
#   - remove vertices and see if trend remains (samples deplete)
# 2) plot igraph results (not width < 1 does not show on pdfs so map r to transpancy)
#
# To-do: reduce redundancy (e.g. A_related)
##############################################################

rm(list = ls())
library(igraph) # To make graph, construct components etc.
library(RColorBrewer)
load('../RData/mles_true.RData')
load('../RData/SNPData.RData') # Load SNP data for cities
load('../RData/geo_dist_info.RData')
source('./igraph_functions.R')
eps = 0.01 # Below which LCI considered close to zero
cols = colorRampPalette(brewer.pal(12, "Paired")) # Function to create colours

FILTER = F # To re-generate filtered results
PDF = T # Set to TRUE to plot graphs to pdf




############################################################
# 1) In this section we filter results, create graphs 
# and memberships etc. 
############################################################
Cities = SNPData$City; names(Cities) = row.names(SNPData) 
cols_cities = brewer.pal(length(unique(Cities)), 'Spectral')
names(cols_cities) = rev(unique(Cities))

#===========================================================
# Filter and save
#===========================================================
if (FILTER) {
  All_results = lapply(c(F,T), rm_highly_related_within, Result = mle_CIs, Cities = Cities)
  All_results[[3]] = mle_CIs
  names(All_results) = c('Filter by vertex', 'Filter by edge', 'Unfiltered')
  
  # Extract summaries
  sapply(All_results, function(x){c('Edge count' = nrow(x), 
                                    'Vertex count' = length(unique(c(x$individual1,x$individual2))))})
  
  # Save list of results filtered and not 
  save(All_results, file = '../RData/All_results.RData')
} else {
  load('../RData/All_results.RData')
}

#===========================================================
# Convert filtered results into graphs
#===========================================================
# Create adjacency matrices 
All_adj_matrix = lapply(All_results, construct_adj_matrix, Entry = 'rhat')

# Create graphs 
All_G = lapply(All_adj_matrix, graph_from_adjacency_matrix, mode='upper', diag=F, weighted=T)

# Add meta data
All_G = lapply(All_G, function(x){
  V(x)$site = SNPData[V(x)$name, 'City']
  V(x)$date = SNPData[V(x)$name, 'COLLECTION.DATE']
  return(x)
})

#===========================================================
# Based on unfiltered results, generate jitter 
# and clonal component membership and colour 
#===========================================================
# Adjaceny matrix unfiltered (same as All_adj_matrix_v[[1]] and All_adj_matrix_e[[1]]
A_high = construct_adj_matrix(mle_CIs, Entry = 'r97.5.')
A_high[A_high < 1-eps] = 0 # Edit s.t. only not stat diff from clonal have weight
G_high = graph_from_adjacency_matrix(A_high, mode='upper', diag=F, weighted=T) # Construct graph 
C_high = components(G_high) # Extract components from graph
M_high = C_high$membership # Extract membership of vertices==

# Layout within sites
Jitter = M_high*10^-3 # For graph layout (function on M)
Jitter[M_high %in% which(C_high$csize < 2)] = -0.15 
Jitter = Jitter + rnorm(length(Jitter), 0 , 0.05)

# Create clonal component colours 
Cols = cols(sum(C_high$csize > 1)) # Enumerate colours
names(Cols) = (1:C_high$no)[C_high$csize > 1] # Names cols agrees with memberships set by igraph

# Create a vector of vertex colours
M_cols = sapply(M_high, function(x){ # Return white for singleton 
  ifelse(C_high$csize[x]==1,'#FFFFFF',Cols[as.character(x)])})
names(M_cols) = names(M_high)

# Extract colours of comparisons for histogram
v_names = do.call(rbind, strsplit(attributes(E(G_high))$vnames, split = "\\|")) # Extract vertices per edge
edge_col_ind = M_cols[v_names[,1]] == M_cols[v_names[,2]] & M_cols[v_names[,1]] != "#000000"
edge_cols = array(M_cols[v_names[,1]][edge_col_ind], dim = sum(edge_col_ind), # Creat a vector grays
                  dimnames = list(apply(v_names[edge_col_ind, ], 1, function(x)paste(sort(x), collapse = '_'))))

save(edge_cols, file = '../RData/edge_cols.RData')




############################################################
# 2) In this section we plot the various graphs and components
############################################################
if(PDF){pdf('../Plots/All_CCs.pdf', height = 8, width = 8)}
par(mfrow = c(1,1), family = 'serif')

#===========================================================
# Plot the relatedness between components
# with edge only if statistically distinguishable from zero
# and edge colour proportion to remaining rhats 
#===========================================================
PYRIMID = F # For pyrimid layout where clones seen in multiple sites are on top
R_comp = rm_highly_related_within(Result = mle_CIs, Edge = F,
                                  Cities = sapply(row.names(SNPData), function(x)return("City")))
A_comp = construct_adj_matrix(R_comp, Entry = 'rhat')
A_comp_related = A_comp = construct_adj_matrix(R_comp, Entry = 'r2.5.')
A_comp[A_comp_related < eps] = 0

G_comp = graph_from_adjacency_matrix(A_comp, mode='upper', diag=F, weighted=T)
V(G_comp)$color = M_cols[V(G_comp)$name]
Comp_G = induced_subgraph(G_comp, vids = which(V(G_comp)$color!="#FFFFFF"))
E(Comp_G)$width <- E(Comp_G)$weight
E(Comp_G)$colour <- 'black'
V(Comp_G)$date <- SNPData[V(Comp_G)$name, 'COLLECTION.DATE']
V(Comp_G)$site <- SNPData[V(Comp_G)$name, 'City']
set.seed(150)

# Create clonal component names in order of date of earliest sample per CC
C_names = paste0('CC',1:length(Cols)) # Create CC names
ordered_date_index <- sort.int(V(Comp_G)$date, index.return = T)$ix
sid_ordered_date = V(Comp_G)$name[ordered_date_index] # reorder sample ID by date 
names(C_names) = as.character(M_high[sid_ordered_date]) # Make sure the CC names are same order or memberships

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


if(PYRIMID){
  
  # Work out where each clone is seen
  cites_per_clone = lapply(V(Comp_G)$color, function(COL){
    inds = mle_CIs$sample_comp %in% names(edge_cols)[edge_cols == COL]
    unique(c(mle_CIs$City1[inds], mle_CIs$City1[inds]))
  })
  
  # Count the number of sites observed per clone
  clone_site_no = sapply(cites_per_clone, length) 
  V(Comp_G)$label = sapply(cites_per_clone, function(x){
    paste(sapply(strsplit(x, split = ''), function(y) paste(y[1:2], collapse = '')), collapse = '+')})
  attr(Comp_G, 'layout') = cbind(X = V(Comp_G)$date, 
                                 Y = clone_site_no + rnorm(length(V(Comp_G)), 0, 0.2)) # Full layout
  plot.igraph(Comp_G, vertex.size = 10, 
       vertex.label.cex = 0.5, 
       layout = attributes(Comp_G)$layout, # uncomment for pyrimid layout 
       vertex.label.color = 'black',
       vertex.frame.color = V(Comp_G)$color,
       edge.color = sapply(E(Comp_G)$weight, adjustcolor, col = 'black'))
  
  # Add dates if pyrimid layout
  unique_yrs = as.numeric(unique(SNPData[V(Comp_G)$name, "Year"]))
  min_yr = min(unique_yrs); max_yr = max(unique_yrs)
  yr01 = (unique_yrs-min_yr)/(max_yr - min_yr)
  axis(side = 1, labels = unique_yrs, las = 1, cex.axis = 0.7, 
       at = -1 + yr01 * 2, line = -1, tick = F)
} else {
  
  set.seed(1)
  plot(Comp_G, layout = layout_with_fr, 
       vertex.size = table(M_cols)[V(Comp_G)$color]+5, 
       vertex.label.cex = 0.5, 
       vertex.label = paste0(C_names[as.character(M_high[V(Comp_G)$name])], 
                             '\n', SNPData[V(Comp_G)$name, 'COLLECTION.DATE']), 
       vertex.label.color = 'black',
       vertex.frame.color = 'white',
       vertex.color = cols_cities[V(Comp_G)$site], 
       edge.color = sapply(weights_rescaled, adjustcolor, col = 'black')
       )
}
range(E(Comp_G)$weight)


# # Legend of names of clonal components
# legend('left', pch = 16, bty = 'n', cex = 0.7, pt.cex = 1.5, col = Cols, 
#        legend = C_names[names(Cols)], inset = -0.1)

# Legend of cities
legend('left', pch = 16, bty = 'n', cex = 0.7, pt.cex = 1.5, col = cols_cities, 
       legend = names(cols_cities), inset = -0.1)


if(PDF){dev.off()}





#===========================================================
# For each site comparison vizualise effect of filter with
# Dates and edges proportional to rhats 
#===========================================================
if(PDF){pdf('../Plots/Graphs_timespace.pdf', height = 5, width = 10)}
par(mfrow = c(1,1), mar = c(3,4,3,0.1), family = 'serif', bg = 'white')
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
    #        pt.bg = Cols[Cols %in% unique(V(G)$color[V(G)$color!="#FFFFFF"])], 
    #        cex = 0.7, pt.cex = 1,
    #        legend = C_names[as.character(names(Cols[Cols %in% unique(V(G)$color[V(G)$color!="#FFFFFF"])]))])
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
         pt.bg = Cols[Cols %in% unique(V(G)$color[V(G)$color!="#FFFFFF"])], cex = 0.8, pt.cex = 1, 
         legend = C_names[as.character(names(Cols[Cols %in% unique(V(G)$color[V(G)$color!="#FFFFFF"])]))])
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
