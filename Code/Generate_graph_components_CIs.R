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
PDF = F # Set to TRUE to plot graphs to pdf
PLOT_GRAPHS = T # Set to TRUE to plot graphs 
eps = 0.01
Thresholds = c(eps, 0.25, 0.5)

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


plot(G_high)

# Layout within sites
Jitter = M_high*10^-3 # For graph layout (function on M)
Jitter[M_high %in% which(C_high$csize < 2)] = -0.15 
Jitter = Jitter + rnorm(length(Jitter), 0 , 0.05)

# Create clonal component colours
Cols = cols(sum(C_high$csize > 1)) # Enumerate colours
names(Cols) = unique((1:C_high$no)[C_high$csize > 1])

# Creat a vector of vertex colours
M_cols = sapply(M_high, function(x){ # Return white for singleton 
  ifelse(C_high$csize[x]==1,'#FFFFFF',Cols[as.character(x)])})
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
par(mfrow = c(1,1))
R_comp = rm_highly_related_within(Result = mle_CIs, Edge = F,
                                  Cities = sapply(row.names(SNPData), function(x)return("City")))
A_comp = construct_adj_matrix(R_comp, Entry = 'rhat')
G_comp = graph_from_adjacency_matrix(A_comp, mode='upper', diag=F, weighted=T)
V(G_comp)$color = M_cols[V(G_comp)$name]
Comp_G = induced_subgraph(G_comp, vids = which(V(G_comp)$color!="#FFFFFF"))
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
plot(Comp_G, vertex.size = 10, 
     vertex.label.cex = 0.5, 
     #layout = attributes(Comp_G)$layout, 
     vertex.label = paste0(SNPData[V(Comp_G)$name, 'City'], '\n', M_high[V(Comp_G)$name], '\n', SNPData[V(Comp_G)$name, 'COLLECTION.DATE']), 
     vertex.label.color = 'black')

unique_yrs = as.numeric(unique(SNPData[V(Comp_G)$name, "Year"]))
min_yr = min(unique_yrs)
max_yr = max(unique_yrs)
yr01 = (unique_yrs-min_yr)/(max_yr - min_yr)
axis(side = 1, labels = unique_yrs, las = 1, cex.axis = 0.7, 
     at = -1 + yr01 * 2, line = -1, tick = F)

# Add dates somehow
legend('bottomleft', pch = 21, bty = 'n', cex = 0.7, pt.bg = Cols, y.intersp = 0.7,
       legend = names(Cols), inset = -0.05)




#-----------------------------------------------------
# Looking specifically at clonal components found across 
# Buenaventura and Tumaco
# What are the expected segment lengths for the three that have relatedness close to 0.5? 
#-----------------------------------------------------
Buenaventura_Tumaco = induced_subgraph(Comp_G, vids = which(grepl("Tu", V(Comp_G)$label) & grepl("Bu", V(Comp_G)$label)))
plot(Buenaventura_Tumaco)
A_Buenaventura_Tumaco = round(as_adjacency_matrix(Buenaventura_Tumaco, type = "lower", attr = 'weight'),3)

# Vizualise the data: 
rep_samples = c('sid319', 'sid235', 'sid321', 'sid138', 'sid310')
unique_M_BT = M_high[rep_samples]
G_clonal_comp_BT = induced_subgraph(All_G$Unfiltered, vids = which(M_high %in% unique_M_BT))
set.seed(2)
plot(G_clonal_comp_BT, vertex.size = 3, vertex.label = NA, vertex.color = M_cols[V(G_clonal_comp_BT)$name], 
     edge.width = E(G_clonal_comp_BT)$weight)
EV_names = do.call(rbind, strsplit(attributes(E(G_clonal_comp_BT))$vname, split = "\\|"))

# Extract average relatedness between clonal components
Mean_sd_matrix = array(NA, dim = c(length(unique_M_BT), length(unique_M_BT), 4), 
                       dimnames = list(unique_M_BT, unique_M_BT, c('mu','k','2.5%','97.5%')))
for(i in 1:(length(unique_M_BT)-1)){
  for(j in (i+1):length(unique_M_BT)){
    members_i = names(which(M_high == unique_M_BT[i])) 
    members_j = names(which(M_high == unique_M_BT[j])) 
    inds = (EV_names[,1] %in% members_i & EV_names[,2] %in% members_j) | EV_names[,1] %in% members_j & EV_names[,2] %in% members_i 
    Mean_sd_matrix[i,j,'mu'] = mean(E(G_clonal_comp_BT)$weight[inds])
    
    # Extract the mean lower confidence interval 
    inds_mle = mle_CIs$sample_comp %in% apply(EV_names[inds,], 1, function(x)paste(x[1],x[2],sep='_')) | mle_CIs$sample_comp %in% apply(EV_names[inds,], 1, function(x)paste(x[2],x[1],sep='_'))
    if(Mean_sd_matrix[i,j,'mu'] == mean(mle_CIs$rhat[inds_mle])){ # check
      Mean_sd_matrix[i,j,'2.5%'] = max(mle_CIs$`2.5%`[inds_mle]) # To see highest lower bound
      Mean_sd_matrix[i,j,'97.5%'] = max(mle_CIs$`97.5%`[inds_mle]) # To see highest lower bound
      Mean_sd_matrix[i,j,'k'] = max(mle_CIs$khat[inds_mle])
    } 
  }
}
round(Mean_sd_matrix, 3)

#-----------------------------------------------------
# Looking specifically at related components found across 
# Buenaventura and Tumaco
#-----------------------------------------------------
# Adjency matrix unfiltered (same as All_adj_matrix_v[[1]] and All_adj_matrix_e[[1]]
inds = mle_CIs$City12 %in% c("Buenaventura_Tumaco")
A_related_BT = construct_adj_matrix(Result = mle_CIs[inds,], Entry = '2.5%')
A_related_BT[is.na(A_related_BT)] = 0 # Edit s.t. only not stat diff from related have weight
A_related_BT[A_related_BT < eps] = 0 # Edit s.t. only not stat diff from related have weight
G_related_BT = graph_from_adjacency_matrix(A_related_BT, mode='upper', diag=F, weighted=T) # Construct graph 
C_related_BT = components(G_related_BT) # Extract components from graph
C_related_BT$no

V(G_related_BT)$site = SNPData[V(G_related_BT)$name, 'City']
V(G_related_BT)$data = SNPData[V(G_related_BT)$name, 'COLLECTION.DATE']

plot(G_related_BT, vertex.size = 3, vertex.label = NA, 
     #vertex.color = cols(11)[C_related_BT$membership], 
     vertex.color = M_cols[V(G_related_BT)$name])

plot(G_related_BT, vertex.size = 3, vertex.label = NA, 
     vertex.color = cols(11)[C_related_BT$membership], 
     layout = cbind(V(G_related_BT)$data, 
                    as.numeric(V(G_related_BT)$site == 'Buenaventura') + Jitter[(V(G_related_BT)$name)]))

#-----------------------------------------------------
# For each site comparison vizualise effect of filter
# This does not take into account uncertainty in rhat
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
      legend('bottomleft', pch = 21, bty = 'n', cex = 0.7, 
             pt.bg = Cols[Cols %in% unique(V(G)$color[V(G)$color!="#FFFFFF"])], y.intersp = 0.7,
             legend = names(Cols[Cols %in% unique(V(G)$color[V(G)$color!="#FFFFFF"])]))
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
  
  
  
  #===============================================
  # Plot without dates and with edges only if
  # statistically distinguishable from zero
  #===============================================
  # Create matrixA_related = construct_adj_matrix(Result = mle_CIs, Entry = '2.5%')
  A_related = construct_adj_matrix(Result = mle_CIs, Entry = 'rhat')
  A_related_lci = construct_adj_matrix(Result = mle_CIs, Entry = '2.5%')
  A_related[A_related_lci < eps] = 0 # Edit s.t. only not stat diff from related have weight
  
  # Create graphs 
  G_related = graph_from_adjacency_matrix(A_related, mode='upper', diag=F, weighted=T)
  V(G_related)$site = SNPData[V(G_related)$name, 'City']
  
  set.seed(1)
  SHAPES = c('circle', 'square')
  
  par(mfrow = c(1,1), mar = c(3,2,3,2), family = 'serif', bg = 'gray')
  for(i in 1:nrow(geo_dist_info$pairwise_site_distance)){
    
    # Extract cities
    s1 = as.character(geo_dist_info$pairwise_site_distance[i,'X1'])
    s2 = as.character(geo_dist_info$pairwise_site_distance[i,'X2'])
    
    # Filter graph by city
    G = extract_city_graph(x = G_related, city1 = s1, city2 = s2)
    E(G)$color[grepl("#FFFFFF", E(G)$color)] = "#FFFFFF" # make opaque
    
    # Plot with legend
    plot.igraph(G, vertex.size = 3, vertex.label = NA, 
                vertex.shape = SHAPES[as.numeric(V(G)$site == s1)+1])
    legend('bottomleft', pch = 23, bty = 'n', inset = -0.05, 
           pt.bg = Cols[Cols %in% unique(V(G)$color[V(G)$color!="#FFFFFF"])], cex = 0.75, y.intersp = 0.7,
           legend = names(Cols[Cols %in% unique(V(G)$color[V(G)$color!="#FFFFFF"])]))
    legend('top', pch = c(0,1), pt.bg = 'white', legend = c(s1,s2), bty = 'n')
    
  }
  if(PDF){dev.off()}
}


# PCA
A_pca = A_related
diag(A_pca) = 1
A_pca[lower.tri(A_pca)] = A_pca[upper.tri(A_pca)]
pcaIBD = prcomp(A_pca, center = TRUE, scale. = TRUE)

par(mfrow = c(1,1))
site_cols = brewer.pal(5, 'Dark2')
names(site_cols) = unique(SNPData[,'City'])
barplot(pcaIBD$sdev[1:10]/sum(pcaIBD$sdev)*100,
        ylab = expression('Percent'~sigma),
        names.arg = paste('PC', 1:10), las = 2)

plot(y = pcaIBD$rotation[,1], x = pcaIBD$rotation[,2], pch = 20, col = site_cols[SNPData[rownames(pcaIBD$x), 'City']])
legend('bottomleft', pch = 20, col = site_cols, legend = names(site_cols), inset = 0.01, cex = 0.7)

# Gravity network 
geo_dist_info$pairwise_site_distance[]
net <- graph.data.frame(geo_dist_info$pairwise_site_distance[,1:2], directed = FALSE)
E(net)$weight = geo_dist_info$pairwise_site_distance$gravity_estimate/100000
as_adjacency_matrix(net, attr = 'weight')
set.seed(5)
plot(net, edge.width = E(net)$weight, vertex.color = site_cols[V(net)$name], vertex.size = 15, vertex.label = NA)
