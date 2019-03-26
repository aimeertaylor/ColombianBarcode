##############################################################
# Script to explore trends in IBD with or without removing 
# highly related samples that are liable to amplify evidence 
# of relatedness between sites. 

# Two options: 
# 1) remove edges and see if trend remains (almost all samples remain)
# 2) remove vertices and see if trend remains (only 87 samples remain)
# Second option amounts to a increasingly more severe definition of a within
# site clique, meaning that lines between sites are increasingly independent
# To-do: results are uninterpretable, without error bars 
##############################################################
rm(list = ls())
require(igraph) # To make graph, construct components etc.
# Load SNP data for cities
load('~/Documents/BroadLaptop/ColombianBarcode/RData/Archive/Result.RData')
load('~/Documents/BroadLaptop/ColombianBarcode/RData/SNPData.RData')
load('~/Documents/BroadLaptop/ColombianBarcode/RData/geo_dist_info.RData')
Cities = SNPData$City; names(Cities) = row.names(SNPData)

#==========================================================
# Mechanism to filter highly related samples within sites
#==========================================================
rm_highly_related_within = function(Threshold, Result, Cities, Edge){
  
  if(Edge){
    # First construct an adjaceny matrix of IBD
    sample_names = unique(c(as.character(Result$sample1), as.character(Result$sample2)))
    sample_count = length(sample_names)
    adj_matrix = array(data = NA, dim = c(sample_count, sample_count), 
                       dimnames = list(sample_names, sample_names))
    for(i in 1:nrow(Result)){
      indi = as.character(Result$sample1[i])
      indj = as.character(Result$sample2[i])
      adj_matrix[indi, indj] = Result$IBD[i]
    }
    
    # Second, modify the adjaceny matrix s.t. 
    # only those greater than or equal to the threshold have a non-zero value
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
          rep_sample = sample(names(cities)[cities == city], 1)
        })
      }
      return(rep_samples)
    }))
    
    # Fifth, restructure Result s.t. only one sample per city per component
    inds = Result$sample1 %in% rep_samples & Result$sample2 %in% rep_samples
    Result_filtered = Result[inds,]
  } else {
    Result_filtered = Result[Result$IBD < Threshold,]
  }
  return(Result_filtered)
}



# Filter results
All_results = lapply(c(1.1,1,0.75,0.5), rm_highly_related_within, 
                     Result, Cities, Edge = T )
sapply(All_results, nrow)
sapply(All_results, function(x){length(unique(c(x$sample1,x$sample1)))})

# Negative trend in distance with and across sites remains
for(i in 1:length(All_results)){
  X = All_results[[i]]
  print(summary(lm(IBD ~ geo_dist, data = X)))
  print(summary(lm(IBD ~ geo_dist, data = X[X$geo_dist > 0,])))
 
}

# Fractions highly related
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




