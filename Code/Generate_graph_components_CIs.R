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
PLOT_GRAPHS = F # Set to TRUE to plot graphs 
eps = 0.01
Thresholds = c(eps, 0.25, 0.5)

# Add site comps
mle_CIs$site_comp = apply(cbind(SNPData$City[mle_CIs$individual1], SNPData$City[mle_CIs$individual2]), 1, function(x){
  paste(sort(x), collapse = "_")
})


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
        rep_samples = sapply(unique(cities), function(city){
          # Pick sample (deterministically for comparison across plots)
          rep_sample = paste0('sid',max(as.numeric(gsub('sid','',names(cities)[cities == city])))) 
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
  Y = as.numeric(V(z)$site == city1) + Jitter[V(z)$name] # vertical layout
  attr(z, 'layout') = cbind(X = V(z)$date, Y) # Full layout
  V(z)$color = M_cols[V(z)$name] # Colour by clonal component
  E(z)$weight[is.na(E(z)$weight)] = 0 # NAs introduced when filter edges 
  E(z)$color = sapply(E(z)$weight, adjustcolor, col = 'white') 
  E(z)$color = 'white' # First assign all to be white
  v_names = do.call(rbind, strsplit(attributes(E(z))$vnames, split = "\\|"))
  edge_col_ind = which(M_cols[v_names[,1]] == M_cols[v_names[,2]] & M_cols[v_names[,1]] != "#FFFFFF")
  if(any(edge_col_ind)){
    E(z)$color[edge_col_ind] = M_cols[v_names[,1]][edge_col_ind]
  }
  E(z)$width <- E(z)$weight
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

# Extract summaries
sapply(All_results, function(x){c('Edge count' = nrow(x), 
                                  'Vertex count' = length(unique(c(x$individual1,x$individual2))))})



#===========================================================
# Analyses 2) Visualise result of filtering
#===========================================================
if(PLOT_GRAPHS){
  
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
  Jitter = sapply(rownames(A_high), function(x)rnorm(1,0,0.05)) # For graph layout
  
  # Create clonal component colours
  Cols = cols(sum(C_high$csize > 1)) # Enumerate colours
  names(Cols) = unique((1:C_high$no)[C_high$csize > 1])
  M_cols = sapply(M_high, function(x){ # Return white for singleton 
    ifelse(C_high$csize[x]==1,'#FFFFFF',Cols[as.character(x)])})
  names(M_cols) = names(M_high)
  
  #-----------------------------------------------------
  # For each site comparison vizualise effect of filter
  # This does not take into account undercertainty in rhat
  #-----------------------------------------------------
  if(PDF){pdf('./Graphs_CIs.pdf', height = 10, width = 9)}
  par(mfrow = c(3,3), mar = c(0.5,1,1.5,1), family = 'serif', bg = 'white')
  for(i in 1:nrow(geo_dist_info$pairwise_site_distance)){
    
    for(l in length(All_G):1){
      
      X = All_G[[l]] # For each strategy
      
      # Extract cities
      s1 = as.character(geo_dist_info$pairwise_site_distance[i,'X1'])
      s2 = as.character(geo_dist_info$pairwise_site_distance[i,'X2'])
      
      # Filter graph by city
      G = extract_city_graph(x = X, city1 = s1, city2 = s2)
      
      # Plot graph
      plot.igraph(G, layout = attributes(G)$layout, vertex.size = 3, vertex.label = NA)
      rect(xleft = par("usr")[1], ybottom = par("usr")[3], xright = par("usr")[2], ytop = par("usr")[4], col = "lightgray")
      plot.igraph(G, layout = attributes(G)$layout, vertex.size = 3, vertex.label = NA, add = T)
      
      # Annotate
      legend('bottomleft', pch = 21, bty = 'n', 
             pt.bg = unique(V(G)$color[V(G)$color!="#FFFFFF"]), cex = 0.75, y.intersp = 0.5,
             legend = rep('', length(unique(V(G)$color))-1))
      mtext(side = 3, line = 0.2, adj = 0, sprintf('%s', names(All_G)[l]))
      mtext(side = 1, sprintf('%s', s2), line = -1.2)
      mtext(side = 3, sprintf('%s', s1), line = -1.2)
    }
  }
  if(PDF){dev.off()}
}



#-----------------------------------------------
# Fractions highly related - to finish
#-----------------------------------------------
par(mfrow = c(3,3), bg = 'white', pty = 's', mar = c(4,4,1,1))
sorted_site_comp = names(sort(geo_dist_info$pairwise_site_distance_all[unique(mle_CIs$site_comp)]))
set.seed(1) # For reproducibility

for(j in length(All_results):1){
  
  X = All_results[[j]]
  
  for(i in 1:length(Thresholds)){
  # Calculate proportion based on mle
  prop_highly_related = sapply(sorted_site_comp, function(x){
    inds = X$site_comp == x
    mean(X$`2.5%`[inds] >= Thresholds[i])
  })
  
  # Bootstrap to get CIs due to different sample counts per site 
  CIs_site_comp = sapply(sorted_site_comp, function(x){
    No_site_comp = sum(X$site_comp == x) # No. of samples to re-draw 
    prop_bootstrapped = sapply(1:100, function(b){
      booti_ind = sample(which(X$site_comp == x), No_site_comp, replace = T)
      mean(X$`2.5%`[booti_ind] >= Thresholds[i])})
    CIs = quantile(prop_bootstrapped, probs = c(0.025, 0.975))
    return(CIs)
  })
  
  # Barplot
  Midpoints = barplot(prop_highly_related, names.arg = sorted_site_comp, cex.names = 0.5,
                      las = 2, density = rep(c(100,25), c(5, 10)), 
                      ylim = c(0,max(CIs_site_comp, prop_highly_related)), 
                      main = paste0(names(All_results)[j], Thresholds[i]))
  # CIs
  segments(x0 = Midpoints[,1], x1 = Midpoints[,1],
           y0 = CIs_site_comp[1, ],
           y1 = CIs_site_comp[2, ], lty = 1)}
}




