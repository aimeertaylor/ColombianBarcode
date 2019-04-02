##############################################################
# Script to explore trends in IBD with or without removing 
# highly related samples that are liable to amplify evidence 
# of relatedness between sites.
#
# In this scipt I explored propogating uncertainty in rhat 
# via parametric bootstrap realisations of rhat. This
# lead to some proportions based on mles < than proportions 
# based on 2.5% quantiles. This is likely due to the fact that 
# most bootstrap rhats exceed mles (with probability 1 if we
# filter mles > threshold). It appears odd, but, intuition 
# regarding bias is notoriously difficult. Since
# I have now moved on to leveraging our knowledge about 
# uncertainty in rhat by exploring edges that are 
# "not statistically different from 0 or 1, I am archiving this
# code (1st April 2019).
# 
# Two options: 
# 1) remove edges and see if trend remains (almost all samples remain)
# 2) remove vertices and see if trend remains (samples deplete)
##############################################################
rm(list = ls())
require(igraph) # To make graph, construct components etc.
require(RColorBrewer)
# Load SNP data for cities
load('../RData/mle_boot.RData')
load('../RData/SNPData.RData')
load('../RData/geo_dist_info.RData')
PDF = F # Set to TRUE to plot graphs to pdf
PLOT_GRAPHS = F # Set to TRUE to plot graphs 

# Add site comps
mle_boot$site_comp = apply(cbind(SNPData$City[mle_boot$individual1], SNPData$City[mle_boot$individual2]), 1, function(x){
  paste(sort(x), collapse = "_")
})


#==========================================================
# Function to create colours
#==========================================================
cols = colorRampPalette(brewer.pal(12, "Paired")) 


#==========================================================
# Function to create an adjacency matrix
#=========================================================
construct_adj_matrix = function(Result){
  sample_names = unique(c(as.character(Result$individual1), as.character(Result$individual2)))
  sample_count = length(sample_names)
  adj_matrix = array(data = NA, dim = c(sample_count, sample_count), 
                     dimnames = list(sample_names, sample_names))
  for(i in 1:nrow(Result)){
    indi = as.character(Result$individual1[i])
    indj = as.character(Result$individual2[i])
    adj_matrix[indi, indj] = Result$rhat[i]
  }
  return(adj_matrix)
}


#==========================================================
# Function to filter highly related samples within sites
#==========================================================
rm_highly_related_within = function(Threshold, Result, Cities, Edge){
  
  if(Edge){
    # Simply remove edges
    Result_filtered = Result[Result$rhat < Threshold,] 
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
    inds = Result$individual1 %in% rep_samples & Result$individual2 %in% rep_samples
    Result_filtered = Result[inds,]
  }
  return(Result_filtered)
}

#==========================================================
# Function to return city subgraph
#==========================================================
extract_city_graph = function(x, city1, city2, high_threshold){
  inds = V(x)$site == city1|V(x)$site == city2 # indices for cities
  z = induced_subgraph(x, vids = which(inds)) # extract subgraph
  Y = as.numeric(V(z)$site == city1) + Jitter[V(z)$name] # vertical layout
  attr(z, 'layout') = cbind(X = V(z)$date, Y) # Full layout
  V(z)$color = M_cols[V(z)$name] # Colour by clonal component
  E(z)$weight[is.na(E(z)$weight)] = 0 # NAs introduced when filter edges 
  E(z)$color = sapply(E(z)$weight, adjustcolor, col = 'white') 
  v_names_high = attributes(E(z))$vnames[E(z)$weight > high_threshold] 
  e_col_high = M_cols[do.call(rbind,strsplit(v_names_high, split = "\\|"))[,1]]
  E(z)$color[E(z)$weight > high_threshold] = e_col_high  # Change to membership
  E(z)$width <- E(z)$weight
  return(z)}

#===========================================================
# Analyses 1a) Plot results
#===========================================================
any(mle_boot$rhat == 1)
par(mfrow = c(2,1), pty = 'm')
hist(mle_boot$rhat, breaks = 101, col = 'red', main = '', freq = F, 
     xlab = expression('Relatedness estimate,'~hat(italic(r))), las = 1)
high_threshold = 0.8 # from histogram
rhat_CIs = t(apply(mle_boot, 1, function(x)quantile(x$rhats_boot, probs=c(0.025, 0.975)))) # Few seconds
Ordered_r = sort.int(mle_boot$rhat, index.return = T)
plot(NULL, ylim = c(0,1), xlim = c(1,length(mle_boot$rhat)), 
     ylab = expression('Relatedness estimate,'~hat(italic(r))), 
     xlab = expression("Sample comparison index ranked by"~hat(italic(r))), 
     bty = 'n', las = 1)

# CIs by segment
segments(x0 = 1:length(mle_boot$rhat), x1 = 1:length(mle_boot$rhat),
         y0 = rhat_CIs[Ordered_r$ix,1], y1 = rhat_CIs[Ordered_r$ix,2],
         col = 'lightgray')
lines(Ordered_r$x)

# # Bootstrap points (huge plot)
# count = 1
# for(s in Ordered_r$ix){
#   points(y = mle_boot[s, 'rhats_boot'][[1]], x = rep(count, 100), 
#          pch = ".", col = rainbow(100))
#   count = count + 1}
# lines(Ordered_r$x)




#===========================================================
# Analyses 1b) Filter results
#===========================================================
# Thresholds to rm edges or vertices and cities of samples
Thresholds = c(1,high_threshold) 
Cities = SNPData$City; names(Cities) = row.names(SNPData) 

# Filter vertices and name
All_results_v = lapply(Thresholds, rm_highly_related_within, Result = mle_boot, Cities, Edge = F)
All_results_e = lapply(Thresholds, rm_highly_related_within, Result = mle_boot, Cities, Edge = T)
names(All_results_v) = Thresholds
names(All_results_e) = Thresholds

# Extract summaries
sapply(All_results_v, function(x){c('Edge count' = nrow(x), 
                                    'Vertex count' = length(unique(c(x$individual1,x$individual2))))})
sapply(All_results_e, function(x){c('Edge count' = nrow(x), 
                                    'Vertex count' = length(unique(c(x$individual1,x$individual2))))})

#===========================================================
# Analyses 2) Visualise result of filtering
#===========================================================
if(PLOT_GRAPHS){
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
  # Adjency matrix unfiltered (same as All_adj_matrix_v[[1]] and All_adj_matrix_e[[1]]
  A_high = construct_adj_matrix(mle_boot)
  A_high[A_high < high_threshold] = 0 # Edit s.t. only clonal edges have non-zero weight
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
  All_G = list('Unfiltered' = All_G_e[1], 
               'Filter by edge' = All_G_e[2], 
               'Filter by vertex' = All_G_v[2])
  
  # Generate results for all city pairs
  if(PDF){pdf('./Graphs.pdf', height = 10, width = 9)}
  par(mfrow = c(3,3), mar = c(0.5,1,1.5,1), family = 'serif', bg = 'white')
  for(i in 1:nrow(geo_dist_info$pairwise_site_distance)){
    
    for(l in 1:length(All_G)){
      
      X = All_G[[l]] # For each strategy
      
      # Extract cities
      s1 = geo_dist_info$pairwise_site_distance[i,'X1']
      s2 = geo_dist_info$pairwise_site_distance[i,'X2']
      
      # Filter graphs by city
      All_G12 = lapply(X, extract_city_graph, city1 = s1, city2 = s2, high_threshold)
      
      # Plot the graph corresponding to each filtered data set
      for(j in 1:length(All_G12)){
        G = All_G12[[j]]
        
        plot.igraph(G, layout = attributes(G)$layout, vertex.size = 3, vertex.label = NA)
        rect(xleft = par("usr")[1], ybottom = par("usr")[3], xright = par("usr")[2], ytop = par("usr")[4], col = "lightgray")
        plot.igraph(G, layout = attributes(G)$layout, vertex.size = 3, vertex.label = NA, add = T)
        
        legend('bottomleft', pch = 21, bty = 'n', 
               pt.bg = unique(V(G)$color[V(G)$color!="#FFFFFF"]), cex = 0.75, y.intersp = 0.5,
               legend = rep('', length(unique(V(G)$color))-1))
        mtext(side = 3, line = 0.2, adj = 0, sprintf('%s', names(All_G)[l]))
        mtext(side = 1, sprintf('%s', s2), line = -1.2)
        mtext(side = 3, sprintf('%s', s1), line = -1.2)
      }
    }
  }
  if(PDF){dev.off()}
}


#-----------------------------------------------
# Fractions highly related - to finish
#-----------------------------------------------
par(mfrow = c(3,3), bg = 'white', pty = 's', mar = c(4,4,1,1))
sorted_site_comp = names(sort(geo_dist_info$pairwise_site_distance_all[unique(mle_boot$site_comp)]))
List = list(All_results_e[[1]], All_results_e[[2]], All_results_v[[2]]) 
set.seed(1) # For reproducibility

for(i in 1:length(List)){
  
  X = List[[i]]
  
  for(threshold in c(0.8,0.5,0.25)){
    
    # Calculate proportion based on mle
    prop_highly_related = sapply(sorted_site_comp, function(x){
      inds = X$site_comp == x
      mean(X$rhat[inds] >= threshold)
    })
    
    # Bootstrap to get CIs due to different sample counts per site (use mles)
    CIs_site_comp = sapply(sorted_site_comp, function(x){
      No_site_comp = sum(X$site_comp == x) # No. of samples to re-draw 
      prop_bootstrapped = sapply(1:100, function(b){
        booti_ind = sample(which(X$site_comp == x), No_site_comp, replace = T)
        mean(X$rhat[booti_ind] >= threshold)})
      CIs = quantile(prop_bootstrapped, probs = c(0.025, 0.975))
      return(CIs)
    })
    
    # Bootstrap to get CIs due to rhat uncertainity (use all samples)
    CIs_rhat = sapply(sorted_site_comp, function(x){
      rhat_boots = X[X$site_comp == x, "rhats_boot"]
      prop_bootstrapped = sapply(1:100, function(b){
        bootb_rhat = sapply(rhat_boots, function(x) x[b])
        mean(bootb_rhat >= threshold)})
      CIs = quantile(prop_bootstrapped, probs = c(0.025, 0.975))
      return(CIs) 
    })
    
    #---------------------------------------------
    # Plot 
    #---------------------------------------------
    YLIM = c(0,max(CIs_rhat, CIs_site_comp, prop_highly_related))
    
    # Barplot
    Midpoints = barplot(prop_highly_related, names.arg = sorted_site_comp, cex.names = 0.5,
                        las = 2, main = threshold, density = rep(c(100,25), c(5, 10)), 
                        ylim = YLIM)
    
    # Uncertainity due to rhats not always less than that due to sample counts
    print(apply(cbind(CIs_site_comp[2, ]-CIs_site_comp[1, ], CIs_rhat[2, ]-CIs_rhat[1, ]), 1, 
                which.min))
    
    # CIs
    segments(x0 = Midpoints[,1]-0.1, x1 = Midpoints[,1]-0.1,
             y0 = CIs_site_comp[1, ],
             y1 = CIs_site_comp[2, ], col = 'blue', lty = 1)
    segments(x0 = Midpoints[,1]+0.1, x1 = Midpoints[,1]+0.1, 
             y0 = CIs_rhat[1, ], 
             y1 = CIs_rhat[2, ], col = 'red', lty = 1)
  }
}




