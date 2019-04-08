##############################################################
# Script to explore trends in IBD with or without removing 
# highly related samples that are liable to amplify evidence 
# of relatedness between sites. 

# Two options: 
# 1) remove edges and see if trend remains (almost all samples remain)
# 2) remove vertices and see if trend remains (samples deplete)

# To-do: 
# remoteness vs transmission: confounding. 
# I think this is the "way forward"
# clear up plots, show to Dan and Caroline
# bootstap by sample rather than sample comp
##############################################################
rm(list = ls())
require(igraph) # To make graph, construct components etc.
require(RColorBrewer)
# Load SNP data for cities
load('../RData/mle_CIs.RData')
load('../RData/SNPData.RData')
load('../RData/geo_dist_info.RData')
PDF = T # Set to TRUE to plot graphs to pdf
PLOT_GRAPHS = T # Set to TRUE to plot graphs 
eps = 0.01
Thresholds = c(eps, 0.25, 0.5)

#++++++++++++++++++++++++++
# mle_CIs = mle_CIs[(mle_CIs$City1 == "Tado" | mle_CIs$City2 == "Tado"), ] 
# mle_CIs[mle_CIs$City2 == "Tado", 'City2'] = 'Quibdo'
# mle_CIs$City12 = apply(mle_CIs[,c('City1','City2')], 1, function(x)paste(sort(x),collapse="_"))
# geo_dist_info$pairwise_site_distance <- geo_dist_info$pairwise_site_distance[-c(1,3,7,9),]
#++++++++++++++++++++++++++
mle_CIs$site_comp = apply(mle_CIs[,c('City1','City2')], 1, function(x)paste(sort(x),collapse="_"))

#==========================================================
# Function to create colours
#==========================================================
cols = colorRampPalette(brewer.pal(12, "Paired")) 

#==========================================================
# Function to create an adjacency matrix
#=========================================================
construct_adj_matrix = function(Result, Entry){
  
  sample_names = unique(c(as.character(Result$individual1), as.character(Result$individual2)))
  sample_count = length(sample_names)
  adj_matrix = array(data = NA, dim = c(sample_count, sample_count), 
                     dimnames = list(sample_names, sample_names))
  for(i in 1:nrow(Result)){
    indi = as.character(Result$individual1[i])
    indj = as.character(Result$individual2[i])
    adj_matrix[indi, indj] = Result[i, Entry]
  }
  return(adj_matrix)
}


#==========================================================
# Function to filter highly related samples within sites
#==========================================================
rm_highly_related_within = function(Result, Cities, Edge){
  
  if(Edge){
    # Keep edges statistically diff. from one
    Result_filtered = Result[Result$`97.5%` < 1-eps,] 
  } else {
    # First construct an adjaceny matrix of 97.5% CI IBD
    adj_matrix = construct_adj_matrix(Result, Entry = '97.5%')
    sample_names = rownames(adj_matrix)
    
    # Second, modify s.t. only edges NOT statistically diff. from one
    adj_matrix_high = adj_matrix
    adj_matrix_high[adj_matrix < 1-eps] = 0
    
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
        # Pick sample with earliest data (deterministically for comparison across plots)
        rep_samples = sapply(unique(cities), function(city){
          samples_to_choose_from = names(cities)[cities == city]
          rep_sample = samples_to_choose_from[which.min(SNPData[samples_to_choose_from, 'COLLECTION.DATE'])]
        })
      }
      return(rep_samples)
    }))
    
    # Fifth, restructure Result s.t. only one sample per city per component
    inds = Result$individual1 %in% rep_samples & Result$individual2 %in% rep_samples
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
  Jitter_sign = rep(1, length(V(z)$site ))
  Jitter_sign[V(z)$site == city2] = -1
  # vertical layout
  Y = as.numeric(V(z)$site == city1) + (Jitter[V(z)$name] * Jitter_sign)
  attr(z, 'layout') = cbind(X = V(z)$date, Y) # Full layout
  V(z)$color = M_cols[V(z)$name] # Colour by clonal component
  E(z)$weight[is.na(E(z)$weight)] = 0 # NAs introduced when filter edges 
  E(z)$width <- E(z)$weight
  E(z)$color = sapply(E(z)$weight, adjustcolor, col = 'white') 
  v_names = do.call(rbind, strsplit(attributes(E(z))$vnames, split = "\\|"))
  edge_col_ind = which(M_cols[v_names[,1]] == M_cols[v_names[,2]] & M_cols[v_names[,1]] != "#FFFFFF")
  if(any(edge_col_ind)){
    E(z)$color[edge_col_ind] = M_cols[v_names[,1]][edge_col_ind]
  }
  return(z)}

#===========================================================
# Plot of CIs 
#===========================================================
Ordered_r = sort.int(mle_CIs$rhat, index.return = T) # Order estimates

# NULL plot
plot(NULL, ylim = c(0,1), xlim = c(1,length(mle_CIs$rhat)), 
     ylab = expression('Relatedness estimate,'~hat(italic(r))), 
     xlab = expression("Sample comparison index ranked by"~hat(italic(r))), 
     bty = 'n', las = 1, panel.first = grid())

# Make transparency a function of lower CI
cols_inds = apply(sapply(1:length(Thresholds), function(i){
  a = rep(1, nrow(mle_CIs))
  ind = mle_CIs$`2.5%` > Thresholds[i]
  a[ind] <- (i+1)
  return(a)}), 1, max)

# Calculate "colour"
cols_CIs = brewer.pal(4, "GnBu")[cols_inds]

# CIs by segment
segments(x0 = 1:length(mle_CIs$rhat), x1 = 1:length(mle_CIs$rhat),
         y0 = mle_CIs$`2.5%`[Ordered_r$ix], y1 = mle_CIs$`97.5%`[Ordered_r$ix],
         col = cols_CIs[Ordered_r$ix])

# Add mles
lines(Ordered_r$x)


#===========================================================
# Analyses 1b) Filter results
#===========================================================
Cities = SNPData$City; names(Cities) = row.names(SNPData) 

# Filter vertices and name
All_results = lapply(c(F,T), rm_highly_related_within, Result = mle_CIs, Cities = Cities)
All_results[[3]] = mle_CIs
names(All_results) = c('Filter by vertex', 'Filter by edge', 'Unfiltered')
save(All_results, file = '../RData/All_results.RData')


# Extract summaries
sapply(All_results, function(x){c('Edge count' = nrow(x), 
                                  'Vertex count' = length(unique(c(x$individual1,x$individual2))))})



#===========================================================
# Analyses 2) Visualise result of filtering
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

#---------------------------------------------------------------------
# Based on unfiltered results, jitter, clonal component membership
#--------------------------------------------------------------------- 
# Adjency matrix unfiltered (same as All_adj_matrix_v[[1]] and All_adj_matrix_e[[1]]
A_high = construct_adj_matrix(mle_CIs, Entry = '97.5%')
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

# Creat a vector of vertex colours
M_cols = sapply(M_high, function(x){ # Return white for singleton 
  ifelse(C_high$csize[x]==1,'#FFFFFF80',Cols[as.character(x)])})
names(M_cols) = names(M_high)

# Extract colours of comparisons for histogram
v_names = do.call(rbind, strsplit(attributes(E(G_high))$vnames, split = "\\|")) # Extract vertices per edge
edge_col_ind = M_cols[v_names[,1]] == M_cols[v_names[,2]] & M_cols[v_names[,1]] != "#FFFFFF"
edge_cols = array(M_cols[v_names[,1]][edge_col_ind], dim = sum(edge_col_ind), # Creat a vector grays
                  dimnames = list(apply(v_names[edge_col_ind, ], 1, function(x)paste(sort(x), collapse = '_'))))
save(edge_cols, file = '../RData/edge_cols.RData')

#-----------------------------------------------------
# Plot the relatedness between components
#-----------------------------------------------------
R_comp = rm_highly_related_within(Result = mle_CIs, Edge = F,
                         Cities = sapply(row.names(SNPData), function(x)return("City")))
A_comp = construct_adj_matrix(R_comp, Entry = 'rhat')
G_comp = graph_from_adjacency_matrix(A_comp, mode='upper', diag=F, weighted=T)
V(G_comp)$color = M_cols[V(G_comp)$name]
Comp_G = induced_subgraph(G_comp, vids = which(V(G_comp)$color!="#FFFFFF80"))
E(Comp_G)$width <- E(Comp_G)$weight
E(Comp_G)$colour <- 'black'
V(Comp_G)$date <- SNPData[V(Comp_G)$name, 'COLLECTION.DATE']
set.seed(150)

# Work out where each clone is seen
cites_per_clone = lapply(V(Comp_G)$color, function(COL){
  inds = mle_CIs$sample_comp %in% names(edge_cols)[edge_cols == COL]
  unique(c(mle_CIs$City1[inds], mle_CIs$City1[inds]))
})

clone_site_no = sapply(cites_per_clone, length) # Number of sites observed
V(Comp_G)$label = sapply(cites_per_clone, function(x){
  paste(sapply(strsplit(x, split = ''), function(y) paste(y[1:2], collapse = '')), collapse = '+')})
attr(Comp_G, 'layout') = cbind(X = V(Comp_G)$date, 
                               Y = clone_site_no + rnorm(length(V(Comp_G)), 0, 0.2)) # Full layout
plot(Comp_G, vertex.size = 10, layout = attributes(Comp_G)$layout, vertex.label.cex = 0.7)
unique_yrs = as.numeric(unique(SNPData[V(Comp_G)$name, "Year"]))
min_yr = min(unique_yrs)
max_yr = max(unique_yrs)
yr01 = (unique_yrs-min_yr)/(max_yr - min_yr)
axis(side = 1, labels = unique_yrs, las = 1, cex.axis = 0.7, 
     at = -1 + yr01 * 2, line = -1, tick = F)



#-----------------------------------------------------
# For each site comparison vizualise effect of filter
# This does not take into account undercertainty in rhat
#-----------------------------------------------------
if(PLOT_GRAPHS){
  if(PDF){pdf('../Plots/Graphs_CIs.pdf', height = 10, width = 8)}
  par(mfrow = c(2,1), mar = c(3,2,3,2), family = 'serif', bg = 'white')
  for(i in 1:nrow(geo_dist_info$pairwise_site_distance)){
    for(l in c(3,1)){
      
      X = All_G[[l]] # For each strategy
      
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
      legend('bottomleft', pch = 21, bty = 'n', 
             pt.bg = sort(unique(V(G)$color[V(G)$color!="#FFFFFF80"])), cex = 1, y.intersp = 0.5,
             legend = rep('', length(unique(V(G)$color))-1))
      mtext(side = 3, line = 0.2, adj = 0, sprintf('%s', names(All_G)[l]), cex = 1.5)
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
}



