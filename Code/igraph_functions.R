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
    Result_filtered = Result[Result$`r97.5.` < 1-eps,] 
  } else {
    # First construct an adjaceny matrix of 97.5% CI IBD
    adj_matrix = construct_adj_matrix(Result, Entry = 'r97.5.')
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
extract_city_graph = function(x, city1, city2, col_edge = 'white', 
                              a = 0, b = 1, c = 2){
  inds = V(x)$site == city1|V(x)$site == city2 # indices for cities
  z = induced_subgraph(x, vids = which(inds)) # extract subgraph
  Jitter_sign = rep(1, length(V(z)$site ))
  Jitter_sign[V(z)$site == city2] = -1
  # vertical layout
  Y = as.numeric(V(z)$site == city1) + (Jitter[V(z)$name] * Jitter_sign)
  attr(z, 'layout') = cbind(X = V(z)$date, Y) # Full layout
  V(z)$color = M_cols[V(z)$name] # Colour by clonal component
  E(z)$weight[is.na(E(z)$weight)] = 0 # NAs introduced when filter edges 
  
  # Scale edge widths (not weights) and transparancy to range a, b, c
  weights_rescaled = a + (b-a) * (E(z)$weight - min(E(z)$weight)) / (max(E(z)$weight) - min(E(z)$weight))
  E(z)$width = weights_rescaled * c
  E(z)$color = sapply(weights_rescaled, function(x) adjustcolor(col_edge, alpha.f = x)) # Make transparency a function of rhat
  v_names = do.call(rbind, strsplit(attributes(E(z))$vnames, split = "\\|"))
  edge_col_ind = which(M_cols[v_names[,1]] == M_cols[v_names[,2]] & M_cols[v_names[,1]] != "#FFFFFF")
  if(any(edge_col_ind)){
    E(z)$color[edge_col_ind] = M_cols[v_names[,1]][edge_col_ind]
  }
  return(z)}
