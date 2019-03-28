##############################################################
# Script to explore trends in IBD with or without removing 
# highly related samples that are liable to amplify evidence 
# of relatedness between sites. 

# Two options: 
# 1) remove edges and see if trend remains (almost all samples remain)
# 2) remove vertices and see if trend remains (only 87 samples remain)
# Second option amounts to a increasingly more severe definition of a within
# site clique, meaning that lines between sites are increasingly independent
# To-do: results are uninterpretable, without error bars, mle_df, 
# remoteness vs transmission: confounding. 
# show plots to Dan and Caroline
##############################################################

rm(list = ls())
require(igraph) # To make graph, construct components etc.
require(RColorBrewer)
# Load SNP data for cities
load('~/Documents/BroadLaptop/ColombianBarcode/RData/Archive/Result.RData')
load('~/Documents/BroadLaptop/ColombianBarcode/RData/SNPData.RData')
load('~/Documents/BroadLaptop/ColombianBarcode/RData/geo_dist_info.RData')
PDF = F # Set to TRUE to plot graphs to pdf

#==========================================================
# Function to create colours
#==========================================================
cols = colorRampPalette(brewer.pal(12, "Paired")) 


#==========================================================
# Function to create an adjacency matrix
#=========================================================
construct_adj_matrix = function(Result){
  sample_names = unique(c(as.character(Result$sample1), as.character(Result$sample2)))
  sample_count = length(sample_names)
  adj_matrix = array(data = NA, dim = c(sample_count, sample_count), 
                     dimnames = list(sample_names, sample_names))
  for(i in 1:nrow(Result)){
    indi = as.character(Result$sample1[i])
    indj = as.character(Result$sample2[i])
    adj_matrix[indi, indj] = Result$IBD[i]
  }
  return(adj_matrix)
}


#==========================================================
# Function to filter highly related samples within sites
#==========================================================
rm_highly_related_within = function(Threshold, Result, Cities, Edge){
  
  if(Edge){
    # Simply remove edges
    Result_filtered = Result[Result$IBD < Threshold,]
  } else {
    # First construct an adjaceny matrix of IBD
    adj_matrix = construct_adj_matrix(Result)
    sample_names = rownames(adj_matrix)
    
    # Second, modify s.t. only edges >= threshold are non-zero
    adj_matrix_high = adj_matrix
    adj_matrix_high[adj_matrix < Threshold] = 0
    
    # Third, create a graph and list its components
    G = graph_from_adjacency_matrix(adjmatrix = adj_matrix_high, 
                                    mode='upper', 
                                    diag=F, weighted = 'w', 
                                    add.colnames = T)
    S = components(graph = G) # Lists components
    
    # Forth, draw a single sample per site within a component 
    rep_samples = unlist(lapply(1:S$no, function(s){
      samples_ids = which(S$membership==s)
      if(length(samples_ids) == 1){
        rep_samples = sample_names[samples_ids] # If only one sample per component, keep
      } else { # Keep one per city
        cities = Cities[sample_names[samples_ids]]
        rep_samples = sapply(unique(cities), function(city){
          # Pick sample deterministically (for comparison across plots)
          rep_sample = paste0('sid',max(as.numeric(gsub('sid','',names(cities)[cities == city])))) 
        })
      }
      return(rep_samples)
    }))
    
    # Fifth, restructure Result s.t. only one sample per city per component
    inds = Result$sample1 %in% rep_samples & Result$sample2 %in% rep_samples
    Result_filtered = Result[inds,]
  }
  return(Result_filtered)
}

#==========================================================
# Function to return city subgraph
#==========================================================
extract_city_graph = function(x, city1, city2){
  inds = V(x)$site == city1|V(x)$site == city2 # indices for cities
  z = induced_subgraph(x, vids = which(inds)) # extract subgraph
  Y = as.numeric(V(z)$site == city1) + Jitter[V(z)$name] # vertical layout
  attr(z, 'layout') = cbind(X = V(z)$date, Y) # Full layout
  V(z)$color = M_cols[V(z)$name] # Colour by clonal component
  E(z)$weight[is.na(E(z)$weight)] = 0 # NAs introduced when filter edges 
  E(z)$color = sapply(E(z)$weight, adjustcolor, col = 'darkgray')
  E(z)$color[E(z)$weight == 1] = 'black'
  E(z)$width <- E(z)$weight
  return(z)}

#===========================================================
# Analyses 
#===========================================================

# Thresholds to rm edges or vertices and cities of samples
Thresholds = c(1.1,1,0.75,0.5) 
Cities = SNPData$City; names(Cities) = row.names(SNPData) 

# Filter vertices and name
All_results_v = lapply(Thresholds, rm_highly_related_within, Result, Cities, Edge = F)
All_results_e = lapply(Thresholds, rm_highly_related_within, Result, Cities, Edge = T)
names(All_results_v) = Thresholds
names(All_results_e) = Thresholds

# Extract summaries
sapply(All_results_v, function(x){c('Edge count' = nrow(x), 
                                    'Vertex count' = length(unique(c(x$sample1,x$sample1))))})
sapply(All_results_e, function(x){c('Edge count' = nrow(x), 
                                    'Vertex count' = length(unique(c(x$sample1,x$sample1))))})

# Create adjacency matrices 
All_adj_matrix_v = lapply(All_results_v, construct_adj_matrix)
All_adj_matrix_e = lapply(All_results_e, construct_adj_matrix)

# Create graphs 
All_G_v = lapply(All_adj_matrix_v, graph_from_adjacency_matrix, mode='upper', diag=F, weighted=T)
All_G_e = lapply(All_adj_matrix_e, graph_from_adjacency_matrix, mode='upper', diag=F, weighted=T)

# Add meta data
All_G_v = lapply(All_G_v, function(x){
  V(x)$site = SNPData[V(x)$name, 'City']
  V(x)$date = SNPData[V(x)$name, 'COLLECTION.DATE']
  return(x)
})
All_G_e = lapply(All_G_e, function(x){
  V(x)$site = SNPData[V(x)$name, 'City']
  V(x)$date = SNPData[V(x)$name, 'COLLECTION.DATE']
  return(x)
})

#---------------------------------------------------------------------
# Based on unfiltered results, jitter, clonal component membership
#--------------------------------------------------------------------- 
# Adjency matrix unfiltered (same as All_adj_matrix_v[[1]] and All_adj_matrix_e[[1]])
A_high = construct_adj_matrix(Result)
A_high[A_high < 1] = 0 # Edit s.t. only clonal edges have non-zero weight
G_high = graph_from_adjacency_matrix(A_high, mode='upper', diag=F, weighted=T) # Construct graph 
C_high = components(G_high) # Extract components from graph
M_high = C_high$membership # Extract membership of vertices
Jitter = sapply(rownames(A_high), function(x)rnorm(1,0,0.05)) # For graph layout

# Create clonal component colours
Cols = cols(sum(C_high$csize > 1)) # Enumerate colours
names(Cols) = unique((1:C_high$no)[C_high$csize > 1])
M_cols = sapply(M_high, function(x){ # Return white for singleton 
  ifelse(C_high$csize[x]==1,'#FFFFFF',Cols[as.character(x)])})
names(M_cols) = names(M_high)





#-----------------------------------------------------
# For each site comparison vizualise effect of filter
#-----------------------------------------------------
All_G = list('Filter by edge' = All_G_e, 'Filter by vertex' = All_G_v)

# Generate results for all city pairs
for(i in 1:nrow(geo_dist_info$pairwise_site_distance)){
  
  for(l in 1:length(All_G)){
    X = All_G[[l]] # For each strategy
    
    # Extract cities
    s1 = geo_dist_info$pairwise_site_distance[i,'X1']
    s2 = geo_dist_info$pairwise_site_distance[i,'X2']
    
    # Filter graphs by city
    All_G12 = lapply(X, extract_city_graph, city1 = s1, city2 = s2)
    
    # Plot the graph corresponding to each filtered data set
    if(PDF){pdf('../Plots/Graphs.pdf')}
    par(mfrow = c(2,2), mar = c(1,3,1,0.2), family = 'serif')
    for(j in 1:length(All_G12)){
      G = All_G12[[j]]
      plot.igraph(G, layout = attributes(G)$layout, vertex.size = 3, vertex.label = NA)
      legend('bottomleft', pch = 21, bty = 'n', inset = -0.05, 
             pt.bg = unique(V(G)$color), cex = 0.75,
             legend = rep('', length(unique(V(G)$color))))
      mtext(side = 2, line = 1.5, adj = 0, 
            sprintf('%s with threshold = %s', names(All_G)[l], names(All_G12)[j]))
      mtext(side = 1, sprintf('%s', s2), line = -1)
      mtext(side = 3, sprintf('%s', s1), line = -1)
      box(which = 'figure')
    }
    if(PDF){dev.off()}
  }
}


#-----------------------------------------------
# Fractions highly related - to finish
#-----------------------------------------------
par(mfrow = c(4,4))
sorted_site_comp = names(sort(geo_dist_info$pairwise_site_distance_all[unique(Result$site_comp)]))

for(i in rev(1:length(All_results))){
  for(threshold in c(1,0.75,0.5,0.25)){
    X = All_results[[i]]
    prop_highly_related = sapply(sorted_site_comp, function(x){
      inds = X$site_comp == x
      mean(X$IBD[inds] >= threshold)
    })
    barplot(prop_highly_related, names.arg = sorted_site_comp, cex.names = 0.5,
            las = 2, main = threshold, 
            density = rep(c(100,25), c(5, 10)))
  }}




