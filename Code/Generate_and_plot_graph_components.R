##############################################################
# Script to 
#
# 1) filter results 
#   - remove edges and see if trend remains (almost all samples remain)
#   - remove vertices and see if trend remains (samples deplete)
# 2) plot igraph results (not width < 1 does not show on pdfs so map r to transpancy)
# 3) Generate and plot PCA and compare with gravity network (maybe remove)
##############################################################

rm(list = ls())
library(igraph) # To make graph, construct components etc.
library(RColorBrewer)
load('../RData/mles_true.RData')
load('../RData/SNPData.RData') # Load SNP data for cities
load('../RData/geo_dist_info.RData')
source('./igraph_functions.R')
PDF = T # Set to TRUE to plot graphs to pdf
eps = 0.01 # Below which LCI considered close to zero
cols = colorRampPalette(brewer.pal(12, "Paired")) # Functn to create colours


############################################################
# 1) In this section we filter results, create graphs 
# and memberships etc. 
############################################################
Cities = SNPData$City; names(Cities) = row.names(SNPData) 
cols_cities = brewer.pal(length(unique(Cities)), 'Spectral')
names(cols_cities) = rev(unique(Cities))

#===========================================================
# Filter 
#===========================================================
All_results = lapply(c(F,T), rm_highly_related_within, Result = mle_CIs, Cities = Cities)
All_results[[3]] = mle_CIs
names(All_results) = c('Filter by vertex', 'Filter by edge', 'Unfiltered')

# Extract summaries
sapply(All_results, function(x){c('Edge count' = nrow(x), 
                                  'Vertex count' = length(unique(c(x$individual1,x$individual2))))})

# Save list of results filtered and not 
save(All_results, file = '../RData/All_results.RData')


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
# Adjency matrix unfiltered (same as All_adj_matrix_v[[1]] and All_adj_matrix_e[[1]]
A_high = construct_adj_matrix(mle_CIs, Entry = 'r97.5.')
A_high[A_high < 1-eps] = 0 # Edit s.t. only not stat diff from clonal have weight
G_high = graph_from_adjacency_matrix(A_high, mode='upper', diag=F, weighted=T) # Construct graph 
C_high = components(G_high) # Extract components from graph
M_high = C_high$membership # Extract membership of vertices

# Layout within sites
Jitter = M_high*10^-3 # For graph layout (function on M)
Jitter[M_high %in% which(C_high$csize < 2)] = -0.15 
Jitter = Jitter + rnorm(length(Jitter), 0 , 0.05)

# Create clonal component colours
Cols = cols(sum(C_high$csize > 1)) # Enumerate colours
names(Cols) = unique((1:C_high$no)[C_high$csize > 1])
C_names = paste0('CC',1:length(Cols))
names(C_names) = names(Cols)

# Create a vector of vertex colours
M_cols = sapply(M_high, function(x){ # Return white for singleton 
  ifelse(C_high$csize[x]==1,'#FFFFFF',Cols[as.character(x)])})
names(M_cols) = names(M_high)

# Extract colours of comparisons for histogram
v_names = do.call(rbind, strsplit(attributes(E(G_high))$vnames, split = "\\|")) # Extract vertices per edge
edge_col_ind = M_cols[v_names[,1]] == M_cols[v_names[,2]] & M_cols[v_names[,1]] != "#FFFFFF"
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
  plot(Comp_G, vertex.size = 10, 
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
  plot(Comp_G, 
       vertex.size = table(M_cols)[V(Comp_G)$color]+5, 
       vertex.label.cex = 0.5, 
       vertex.label = paste0(C_names[as.character(M_high[V(Comp_G)$name])], 
                             '\n', SNPData[V(Comp_G)$name, 'City'], 
                             '\n', SNPData[V(Comp_G)$name, 'COLLECTION.DATE']), 
       vertex.label.color = 'black',
       vertex.frame.color = 'white',
       vertex.color = cols_cities[V(Comp_G)$site], 
       edge.color = sapply(E(Comp_G)$weight, adjustcolor, col = 'black'))
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
# Dates and edges proportional to rhats (regardless if statistically significant or not)
#===========================================================
if(PDF){pdf('../Plots/Graphs_timespace.pdf', height = 8, width = 8)}
par(mfrow = c(2,1), mar = c(3,2,3,2), family = 'serif', bg = 'white')
filter_name = c("Unfiltered" = "Unfiltered", "Filter by vertex" = "Filter CCs")
for(i in 1:nrow(geo_dist_info$pairwise_site_distance)){
  for(l in c("Unfiltered","Filter by vertex")){
    
    # Extract graph
    X = All_G[[l]] 
    
    # Extract cities
    s1 = as.character(geo_dist_info$pairwise_site_distance[i,'X1'])
    s2 = as.character(geo_dist_info$pairwise_site_distance[i,'X2'])
    
    # Filter graph by city
    G = extract_city_graph(x = X, city1 = s1, city2 = s2)
    
    # Plot graph
    plot.igraph(G, layout = attributes(G)$layout, vertex.size = 3, vertex.label = NA, asp=0) # For layout
    rect(xleft = par("usr")[1], ybottom = par("usr")[3], xright = par("usr")[2], ytop = par("usr")[4], col = "lightgray")
    plot.igraph(G, layout = attributes(G)$layout, vertex.size = 3, vertex.label = NA, add = T)
    
    # Overlay coloured only
    E(G)$color[grepl("#FFFFFF", E(G)$color)] = NA
    plot.igraph(G, layout = attributes(G)$layout, vertex.size = 3, vertex.label = NA, add = T)
    
    # Annotate
    legend('bottomleft', pch = 21, bty = 'n', cex = 0.5, 
           pt.bg = Cols[Cols %in% unique(V(G)$color[V(G)$color!="#FFFFFF"])], y.intersp = 0.7,
           col = Cols[Cols %in% unique(V(G)$color[V(G)$color!="#FFFFFF"])], 
           legend = C_names[as.character(names(Cols[Cols %in% unique(V(G)$color[V(G)$color!="#FFFFFF"])]))])
    mtext(side = 3, line = 0.2, adj = 0, sprintf('%s', filter_name[l]), cex = 1.5)
    mtext(side = 1, sprintf('%s', s2), line = 1.5, cex = 1.5)
    mtext(side = 3, sprintf('%s', s1), line = 0.5, cex = 1.5)
    
    # Add years
    unique_yrs = as.numeric(unique(SNPData[V(G)$name, "Year"]))
    min_yr = min(unique_yrs)
    max_yr = max(unique_yrs)
    yr01 = (unique_yrs-min_yr)/(max_yr - min_yr)
    axis(side = 1, labels = unique_yrs, las = 1, cex.axis = 0.7, at = -1 + yr01 * 2,
         line = -1, tick = F)
  }
}

if(PDF){dev.off()}

#===============================================
# Plot without dates and with edges only if
# statistically distinguishable from zero with transparency
# proportional to rhat
#===============================================
if(PDF){pdf('../Plots/Graphs.pdf', height = 8, width = 8)}
par(mfrow = c(1,1), mar = c(0,3,0,0), family = 'serif', bg = 'white')
# Create matrixA_related = construct_adj_matrix(Result = mle_CIs, Entry = 'r2.5.')
A_related = construct_adj_matrix(Result = mle_CIs, Entry = 'rhat')
A_related_lci = construct_adj_matrix(Result = mle_CIs, Entry = 'r2.5.')
A_related[A_related_lci < eps] = 0 # Edit s.t. only not stat diff from related have weight

# Create graphs 
G_related = graph_from_adjacency_matrix(A_related, mode='upper', diag=F, weighted=T)
V(G_related)$site = SNPData[V(G_related)$name, 'City']

set.seed(1)
SHAPES = c('circle', 'square') # Different shapes distinguish different sites 

for(i in 1:nrow(geo_dist_info$pairwise_site_distance)){
  
  # Extract cities
  s1 = as.character(geo_dist_info$pairwise_site_distance[i,'X1'])
  s2 = as.character(geo_dist_info$pairwise_site_distance[i,'X2'])
  
  # Filter graph by city
  G = extract_city_graph(x = G_related, city1 = s1, city2 = s2, 'black')
  
  # Plot with legend
  set.seed(1)
  plot.igraph(G, vertex.size = 3, vertex.label = NA, 
              vertex.shape = SHAPES[as.numeric(V(G)$site == s1)+1])
  legend('left', pch = 23, bty = 'n', inset = -0.07, 
         pt.bg = Cols[Cols %in% unique(V(G)$color[V(G)$color!="#FFFFFF"])], cex = 0.7, pt.cex = 1, 
         legend = C_names[as.character(names(Cols[Cols %in% unique(V(G)$color[V(G)$color!="#FFFFFF"])]))])
  legend('bottom', pch = c(0,1), pt.bg = 'white', legend = c(s1,s2), bty = 'n')
}


#===========================================================
# Looking specifically at clonal components found across 
# Buenaventura and Tumaco (some of also found in Guapi)
# Edge transparency proportion to rhat if statistically
# distinguishable from zero (since based on Comp_G)
#===========================================================
par(bg = 'white')
G_BT = induced_subgraph(Comp_G, vids = which(grepl("Tu", V(Comp_G)$label) & grepl("Bu", V(Comp_G)$label)))
A_BT = round(as_adjacency_matrix(G_BT, type = "lower", attr = 'weight'),3)
rep_samples = c('sid319', 'sid235', 'sid321', 'sid138', 'sid310') # From A_BT
unique_M_BT = M_high[rep_samples] # Extract the membership numbers

# Extract all those with said membership (includes some from guapi)
G_clonal_comp_BT = induced_subgraph(All_G$Unfiltered, vids = which(M_high %in% unique_M_BT)) 

# Plot all those with said membership
set.seed(2)
plot(G_clonal_comp_BT, vertex.size = 3, vertex.label = NA, 
     vertex.color = M_cols[V(G_clonal_comp_BT)$name], 
     edge.color = sapply(E(G_clonal_comp_BT)$weight, adjustcolor, col = 'black')) # Make transparency a function of rhat)

# Extract edge vertex names as matrix
EV_names = do.call(rbind, strsplit(attributes(E(G_clonal_comp_BT))$vname, split = "\\|"))

# Extract average relatedness between clonal components (includes some guapi)
Mean_sd_matrix = array(NA, dim = c(length(unique_M_BT), length(unique_M_BT), 4), 
                       dimnames = list(C_names[as.character(unique_M_BT)], C_names[as.character(unique_M_BT)], c('mu','k','r2.5.','r97.5.')))
for(j in 1:(length(unique_M_BT)-1)){
  for(i in (j+1):length(unique_M_BT)){
    members_i = names(which(M_high == unique_M_BT[i])) 
    members_j = names(which(M_high == unique_M_BT[j])) 
    inds = (EV_names[,1] %in% members_i & EV_names[,2] %in% members_j) | EV_names[,1] %in% members_j & EV_names[,2] %in% members_i 
    Mean_sd_matrix[i,j,'mu'] = mean(E(G_clonal_comp_BT)$weight[inds])
    
    # Extract the mean lower confidence interval 
    inds_mle = mle_CIs$sample_comp %in% apply(EV_names[inds,], 1, function(x)paste(x[1],x[2],sep='_')) | mle_CIs$sample_comp %in% apply(EV_names[inds,], 1, function(x)paste(x[2],x[1],sep='_'))
    if(Mean_sd_matrix[i,j,'mu'] == mean(mle_CIs$rhat[inds_mle])){ # check
      Mean_sd_matrix[i,j,'r2.5.'] = max(mle_CIs$`r2.5.`[inds_mle]) # To see highest lower bound
      Mean_sd_matrix[i,j,'r97.5.'] = max(mle_CIs$`r97.5.`[inds_mle]) # To see highest lower bound
      Mean_sd_matrix[i,j,'k'] = max(mle_CIs$khat[inds_mle])
    } 
  }
}
round(Mean_sd_matrix, 3)



#########################################################
# 3) In this section we generate and plot PCA and 
# compare with gravity network
#########################################################
par(bg = 'white', mar = c(5,4,1,1), mfrow = c(2,1))
A_pca = A_related
diag(A_pca) = 1
A_pca[lower.tri(A_pca)] = A_pca[upper.tri(A_pca)]
pcaIBD = prcomp(A_pca, center = TRUE, scale. = TRUE)

# Plot PCA
site_cols = brewer.pal(5, 'Dark2')
names(site_cols) = unique(SNPData[,'City'])
# barplot(pcaIBD$sdev[1:10]/sum(pcaIBD$sdev)*100,
#         ylab = expression('Percent'~sigma),
#         names.arg = paste('PC', 1:10), las = 2)

plot(y = pcaIBD$rotation[,1], x = -pcaIBD$rotation[,2], 
     ylab = sprintf('PCA 1 (%s)', round(pcaIBD$sdev[1]/sum(pcaIBD$sdev),3)), 
     xlab = sprintf('-PCA 2 (%s)', round(pcaIBD$sdev[2]/sum(pcaIBD$sdev),3)), 
     pch = 20, col = site_cols[SNPData[rownames(pcaIBD$x), 'City']])
legend('bottomleft', pch = 20, col = site_cols, legend = names(site_cols), inset = 0.01)

# PCA appears to be similar to the gravity network 
net <- graph.data.frame(geo_dist_info$pairwise_site_distance[,1:2], directed = FALSE)
E(net)$weight = geo_dist_info$pairwise_site_distance$gravity_estimate/100000
as_adjacency_matrix(net, attr = 'weight')
set.seed(5)
plot(net, edge.width = E(net)$weight, vertex.frame.color = site_cols[V(net)$name],
     vertex.color = site_cols[V(net)$name], vertex.size = 15, vertex.label = NA)


if(PDF){dev.off()}