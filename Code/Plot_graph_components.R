##############################################################
# Script to plot igraph results 
#
# Re related edge viz (unrelated edges not plotted): 
# Variation in edge widths between 0 and 1 do not show on pdfs 
# Edge transparancy between 0 and 1 does not show when printed
# A.s. use both: e.g. map related edge transparency from a to b 
# and plot width as transparency * c
#
# Not A_related multiply (4 times) defined as sometimes but not 
# always based on filtered and sometimes but not always unrelated 
# edges set to zero - could name differently but done as of yet
##############################################################

rm(list = ls())
library(igraph) # To make graph, construct components etc.
library(RColorBrewer) # For colours
library(kableExtra) # Needed for kable
load('../RData/All_results.RData') # Load All_results
load('../RData/SNPData.RData') # Load SNP data for cities
load('../RData/geo_dist_info_cities.RData')
source('./igraph_functions.R')
eps <- 0.01 # Below which LCI considered close to zero
cols <- colorRampPalette(brewer.pal(12, "Paired")) # Function to create colours
PDF <- T # Set to TRUE to plot graphs to pdf
a = 0.05; b = 1; c = 2 # Scaling parameters for edge width and transparancy 

# 5 x 1 vector of city colours named by city
cities <- unique(SNPData$City)
cols_cities <- brewer.pal(length(cities), 'Spectral') 
names(cols_cities) <- rev(cities) 

# Add accents
cols_cities_text <- gsub('Quibdo', 'Quibdó', 
                         gsub('Tado', 'Tadó', names(cols_cities)))



#===========================================================
# Get clonal component (CC) membership for all 325 samples
#===========================================================
#-----------------------------------------------------------
# Create graph of all samples to extract CCs
A_high <- construct_adj_matrix(All_results$Unfiltered, Entry = 'r97.5.') # Adj. matrix unfiltered using 'r97.5.' 
A_high[A_high < 1-eps] <- 0 # Edit s.t. only clonal have weight
G_high <- graph_from_adjacency_matrix(A_high, mode='upper', diag=F, weighted=T) # Construct graph 
C_high <- components(G_high) # Extract CCs from graph
M_high <- C_high$membership # Extract CC membership per sample (numeric name of each CC)
#----------------------------------------------------------

#----------------------------------------------------------
# Create a vector of all CCs with 2+ samples and assign each a unique colour 
CCs <- cols(sum(C_high$csize > 1)) # First get unique colours 
names(CCs) <- (1:C_high$no)[C_high$csize > 1] # Name CCs using numeric CC name given igraph
#----------------------------------------------------------

#----------------------------------------------------------
# Create CC character names and order by earliest detected sample per CC
CC_chr_names <- paste0('CC',1:sum(C_high$csize > 1)) 
# Map CC character names to dates of earliest detected sample per CC
# Remove all duplicate CCs regardless of city using a dummy city name
Rhat_filtCC = rm_highly_related_within(Result = All_results$Unfiltered, Edge = F, 
                                       Cities = sapply(row.names(SNPData), function(x)"DummyCity"))
# Extract remaining sample ids
sids_noCCdupl = unique(c(as.character(Rhat_filtCC$individual1), 
                         as.character(Rhat_filtCC$individual2))) # removes 125 samples

# Keep only those in a CC with one or more sample
sids_1stperCC = sids_noCCdupl[M_high[sids_noCCdupl] %in% names(CCs)]

# Extract dates for filtered sample ids 
sample_dates <- SNPData[sids_1stperCC, 'COLLECTION.DATE'] # extract dates 
order_of_sample_dates <- sort.int(sample_dates, index.return = T)$ix # extract order 
sids_1stperCC_ordered <- sids_1stperCC[order_of_sample_dates] # reorder sample ID by date 
names(CC_chr_names) = as.character(M_high[sids_1stperCC_ordered]) # Ensure CC names are ordered as memberships
#--------------------------------------------------------

#----------------------------------------------------------
# Create a list with CC members for each cc and assess cliques
# Added Jan 2021
Clonal_components <- lapply(names(CC_chr_names), function(cc){
  sids <- names(which(M_high == as.numeric(cc)))
  as.character(SNPData[sids, "SAMPLE.CODE"])
})

# Rename Clonal_components
names(Clonal_components) <- CC_chr_names

# Reorder Clonal_components and save
save(Clonal_components, file = "../RData/Clonal_components.RData")

# Check to see how many cliques there are per clonal_component
clique_count_per_component <- sapply(names(CC_chr_names), function(cc){
  sids <- names(which(M_high == as.numeric(cc)))
  subG <- subgraph(G_high, sids)
  x <- maximal.cliques.count(subG, min = 2)
  # Plot, remembering this is after setting edges to zero if not > 0.99
  if(x != 1) plot(subG, main = CC_chr_names[cc])
  return(x)
})
names(clique_count_per_component) <- CC_chr_names

# Three components are not cliques: 2, 15, 18
which(clique_count_per_component != 1)
#----------------------------------------------------------


#--------------------------------------------------------
# Extract table of site and date of earliest sample per CC 
supp_table = data.frame(CC = CC_chr_names, 
                        Earliest_sample = sids_1stperCC_ordered, 
                        Date_earliest_sample = as.character(SNPData[sids_1stperCC_ordered, 'COLLECTION.DATE']), 
                        Site_earliest_sample = as.character(SNPData[sids_1stperCC_ordered, 'City']))
rownames(supp_table) <- names(CC_chr_names)
#--------------------------------------------------------


#--------------------------------------------------------
# In response to reviewer #1's comment about the longevity of CCs from Tado and Quibdo
# Extract the longevity of each CC (and all cities detected) to add to supplementary table 
# And to plot by city

# Start with a vectors whose cc names are the numeric values assigned by igraph
CC_longevities <- array(dim = length(CC_chr_names),  
                        dimnames = list(names(CC_chr_names)))
CC_cities <- array(dim = length(CC_chr_names),  
                   dimnames = list(names(CC_chr_names)))
CC_sizes <- array(dim = length(CC_chr_names),  
                  dimnames = list(names(CC_chr_names)))

for(cc_no_chr in names(CC_chr_names)){ # For each CC with 2 or more samples
  sids <- names(M_high[M_high == as.numeric(cc_no_chr)]) # Extract samples that belong to said cc
  if (supp_table[cc_no_chr,"Earliest_sample"] != sids[which.min(SNPData[sids, 'COLLECTION.DATE'])]) {
    stop("First sample per CC mismatch")
  } else {
    date_range <- range(SNPData[sids, 'COLLECTION.DATE']) # Extract min and max date
    CC_longevities[cc_no_chr] <- diff(date_range) # Calculate the difference in days
    CC_cities[cc_no_chr] <- paste(sort(unique(SNPData[sids, 'City'])), collapse = " ")
    CC_sizes[cc_no_chr] <- length(unique(sids))
  }
}

supp_table$All_cities_detected <- CC_cities[rownames(supp_table)]
supp_table$Sample_count <- CC_sizes[rownames(supp_table)]
supp_table$Longevity <- CC_longevities[rownames(supp_table)]
#--------------------------------------------------------


#--------------------------------------------------------
# # Plot longevities for different cities (exc. CCs detected in multiple cities)
# # - not a great way to visualise
# par(mfrow = c(c(5,1)), mar = c(4,4,0,0))
# for(city in cities){
#   inds <- supp_table$All_cities_detected == city
#   longevities <- supp_table[inds, "Longevity"]
#   barplot(sort(longevities, decreasing = T), col = cols_cities[city], 
#           ylim = range(supp_table$Longevity))
# }

# Plot longevities on log scale and inc. CCs detected in multiple cities
longevities_sorted <- sort.int( supp_table$Longevity, decreasing = T, index.return = T)
cols <- cols_cities[supp_table$All_cities_detected[longevities_sorted$ix]]
cols[is.na(cols)] <- "#BEBEBEFF" 
x <- barplot(height = supp_table$Longevity[longevities_sorted$ix], 
             col = cols, las = 1, log = "y", yaxp = c(1,1,1), 
             ylab = "Longevity (days)", xlab = "Clonal components", xaxt = "n")
text(labels = supp_table$CC[longevities_sorted$ix], pos = 1, offset = 1, 
     x = x-0.8, y = 1, xpd = TRUE, srt = 50, cex = 0.5)
legend("topright", fill = c("#BEBEBEFF",cols_cities),
       legend = c("Multiple cities",cols_cities_text), bty = 'n')

# Plot longevities versus sample count
plot(y = supp_table$Longevity[longevities_sorted$ix], 
     x = supp_table$Sample_count[longevities_sorted$ix], 
     col = cols, pch = 20, 
     ylab = "Longevity (days)", xlab = "Clonal component sample count")

# Positively correlated
cor.test(supp_table$Longevity, 
         supp_table$Sample_count)

# Compute average longevity per city / multiple cities
Average_longevities <- data.frame(sapply(cities, function(city){
  inds <- supp_table$All_cities_detected == city
  c(mean(supp_table$Longevity[inds]), # regardless of sample count
    mean(supp_table$Longevity[inds]/supp_table$Sample_count[inds])) # per sample
}))

# Add averaged for multiple cities
inds <- !supp_table$All_cities_detected %in% cities
Average_longevities$`Muliple cities` <- c(mean(supp_table$Longevity[inds]), 
                                          mean(supp_table$Longevity[inds]/supp_table$Sample_count[inds]))

# Either way, CCs in Tumaco, Buenaventura and Multiple cities have larger longevities 
# than Tado and Quido: 
sort(Average_longevities[1,]) # Average longevity of CC
sort(Average_longevities[2,]) # Average longevity of CC per sample
#--------------------------------------------------------

#--------------------------------------------------------
# Print supplementary table (add accents in tex)
cols_to_inc <- c("CC", "Sample_count", "Longevity", "Date_earliest_sample", "Site_earliest_sample")
kableExtra::kable(supp_table[,cols_to_inc], row.names = FALSE, format = "latex")
#--------------------------------------------------------

#--------------------------------------------------------
# Create a vector of vertex colours by CC - needed for city pair plots
M_cols = sapply(M_high, function(x){ # Return white for singleton 
  ifelse(C_high$csize[x]==1,'#FFFFFF',CCs[as.character(x)])})
#--------------------------------------------------------

#--------------------------------------------------------
# Adjacency over all related with rhat statistically distinguishable from zero
# (no filter)
A_related = construct_adj_matrix(All_results$Unfiltered, Entry = 'rhat')
A_low = construct_adj_matrix(All_results$Unfiltered, Entry = 'r2.5.')
A_related[A_low < eps] = 0 # Edit s.t. only not stat diff from clonal have weight
#--------------------------------------------------------


#==================================================================================
# CC pie-graph whose edges are weighted according to average relatedness 
# averaged of those that are statistically distinguishable from zero
#==================================================================================.
if(PDF){pdf(sprintf('../Plots/All_CCs.pdf'), height = 8, width = 8)}
par(mfrow = c(1,1), family = 'serif')

for(Remove_singletons in c(TRUE,FALSE)){ # If true, remove CCs with only one parasite sample 
  
  # For each pair of CCs with 2 or more samples, extract samples and calculate average relatedness
  CCs_inc <- if(Remove_singletons){names(CCs)}else{as.character(1:C_high$no)}
  CCpairs <- gtools::combinations(n = length(CCs_inc), r = 2, v = CCs_inc) # Get CC pairs
  colnames(CCpairs) = c('individual1','individual2')
  
  # Average over relatedness estimates that are statistically distinguishable from zero
  R_comp <- data.frame(CCpairs, 'rhat_av' = NA, stringsAsFactors = F) 
  for(i in 1:nrow(R_comp)){
    
    # Extract samples per CC
    samples_CC1 <- names(M_high)[as.character(M_high) == R_comp$individual1[i]]
    samples_CC2 <- names(M_high)[as.character(M_high) == R_comp$individual2[i]]
    
    # Construct sample pairs that bridge CCs
    sample_pairs_between_CCs <- expand.grid(samples_CC1, samples_CC2)
    
    # Extract rhats 
    rhats = apply(sample_pairs_between_CCs, 1, function(sids){
      rhat = A_related[sids[1], sids[2]]
      # Thus if rhat is na try swapping rows and columns:
      if(is.na(rhat)){rhat <- A_related[sids[2], sids[1]]}
      return(rhat)
    })
    
    # Average rhats
    R_comp[i, 'rhat_av'] = mean(rhats) 
  }
  
  # Create adjacency matrix over clonal components
  Comp_A = construct_adj_matrix(R_comp, Entry = 'rhat_av')
  # Create graph
  Comp_G = graph_from_adjacency_matrix(Comp_A, mode='upper', diag=F, weighted=T)
  
  # Aside: related components among CCs 
  if(Remove_singletons){
    writeLines(sprintf("Among the clonal components, there are %s related components size %s", 
                       components(Comp_G)$no, 
                       paste(components(Comp_G)$csize, collapse = " and ")))
    
    # Find and inspect parasite members of the anamolous CC
    anomaly_CC_sID1 = names(which(components(Comp_G)$membership != 1)) # The sample ID retained in Comp_G
    anomaly_CC_sIDs = names(which(M_high == anomaly_CC_sID1)) # All sample IDs
    
    
    # Extract and print to screen their data (check with Diego )
    anomaly_CC_SNPData = SNPData[anomaly_CC_sIDs,]
    writeLines(sprintf("The anomalous CC, %s, contains %s samples:", 
                       CC_chr_names[anomaly_CC_sID1], 
                       length(anomaly_CC_sIDs)))
    print(anomaly_CC_SNPData[,1:5])
  }
  
  # For each CC with 2 or more samples, extract number of samples per city (order consistently)
  sample_count_per_city_per_CC <- lapply(CCs_inc, function(CC){
    
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
  names(sample_count_per_city_per_CC) = CCs_inc
  
  # Print non-zero edge weights before rescaling
  writeLines(sprintf("Relatedness values of edges range from %s to %s", 
                     round(min(E(Comp_G)$weight),3), round(max(E(Comp_G)$weight),3)))
  
  # Scale edge widths and transparancy to range a, b, c
  weights_rescaled = a + (b-a) * (E(Comp_G)$weight - min(E(Comp_G)$weight)) / (max(E(Comp_G)$weight) - min(E(Comp_G)$weight))
  vertex_size = sapply(sample_count_per_city_per_CC[V(Comp_G)$name], sum)
  vertex_size[vertex_size > 1] <- vertex_size[vertex_size > 1] + 5
  vertex_size[vertex_size == 1] <- vertex_size[vertex_size == 1] + 2
  
  # Plot as pie (Tim's suggestion)
  set.seed(1)
  plot(Comp_G, 
       layout = layout_with_fr, 
       vertex.shape = "pie", 
       vertex.pie = sample_count_per_city_per_CC[V(Comp_G)$name],
       vertex.pie.color = list(cols_cities[cities]),
       vertex.size = vertex_size, 
       vertex.label = CC_chr_names[V(Comp_G)$name], 
       vertex.label.cex = 0.35, 
       vertex.label.color = 'black',
       vertex.frame.color = NA,
       edge.width = weights_rescaled * c, 
       edge.color = sapply(weights_rescaled, function(x)adjustcolor('black', alpha.f = x)))
  
  # Legend of cities
  legend('bottomleft', pch = 16, bty = 'n', cex = 0.7, pt.cex = 1.5, col = cols_cities, 
         legend = cols_cities_text, inset = 0.1)
}
if(PDF){dev.off()}






#############################################################
# City pair plots
#############################################################

# Add jitter for layout within cities
Jitter = M_high*10^-3 # For graph layout (function on M)
Jitter[M_high %in% which(C_high$csize < 2)] = -0.15 
Jitter = Jitter + rnorm(length(Jitter), 0 , 0.05)

#===========================================================
# For each site comparison vizualise effect of filter with
# dates and edges proportional to related rhats 
#===========================================================
if(PDF){pdf('../Plots/Graphs_timespace.pdf', height = 5, width = 10)}
par(mfrow = c(1,1), mar = c(3,4,3,0.1), family = 'serif', bg = 'white')
filter_name = c("Unfiltered" = "Unfiltered", "Filter by vertex" = "Filter CCs")
SHAPES = c('circle', 'square') # Different shapes distinguish different sites 
for(l in c("Unfiltered","Filter by vertex")){
  
  # Reconstruct adjacency since using both unfiltered and filtered
  X = All_results[[l]]
  A_related_lci = construct_adj_matrix(Result = X, Entry = 'r2.5.')
  A_related = construct_adj_matrix(Result = X, Entry = 'rhat')
  A_related[A_related_lci < eps] = 0 # Edit s.t. unrelated not plotted
  G_related = graph_from_adjacency_matrix(A_related, mode='upper', diag=F, weighted=T)
  V(G_related)$site = SNPData[V(G_related)$name, 'City'] # Add meta data
  V(G_related)$date = SNPData[V(G_related)$name, 'COLLECTION.DATE'] # Add meta data
  
  for(i in 1:nrow(geo_dist_info$pairwise_site_distance)){
    
    # Extract cities
    s1 = as.character(geo_dist_info$pairwise_site_distance[i,'X1'])
    s2 = as.character(geo_dist_info$pairwise_site_distance[i,'X2'])
    
    # Filter graph by city
    G = extract_city_graph(x = G_related, city1 = s1, city2 = s2, col_edge = 'black', 
                           a = a, b = b, c = c)
    
    # Calculate years
    colldates = SNPData[V(G)$name, "COLLECTION.DATE"]
    min_dt = min(colldates)
    max_dt = max(colldates)
    years1stjan = seq(as.Date(round.POSIXt(min_dt, units = "years")), max_dt, by = 'year') 
    yr01 = as.numeric(years1stjan-min_dt)/as.numeric(max_dt - min_dt)
    
    # Plot graph
    plot.igraph(G, layout = attributes(G)$layout, 
                vertex.size = 0, 
                vertex.label = NA, 
                edge.color = NA,
                asp = 0) # For layout
    abline(v = -1 + yr01 * 2, lty = 'dotted', col = 'lightgray') # Add grid line
    plot.igraph(G, layout = attributes(G)$layout, add = T, 
                vertex.shape = SHAPES[as.numeric(V(G)$site == s1)+1], 
                vertex.size = 3, 
                vertex.label = NA) # For layout
    
    # Overlay coloured only
    E(G)$color[grepl("#000000", E(G)$color)] = NA
    plot.igraph(G, layout = attributes(G)$layout, add = T,
                vertex.shape = SHAPES[as.numeric(V(G)$site == s1)+1], 
                vertex.size = 3, 
                vertex.label = NA)
    
    # Annotate
    mtext(side = 3, line = 2, adj = 0, sprintf('%s', filter_name[l]), cex = 1)
    legend('left', pch = 23, bty = 'n', inset = -0.07,
           pt.bg = CCs[CCs %in% unique(V(G)$color[V(G)$color!="#FFFFFF"])],
           cex = 0.7, pt.cex = 1,
           legend = CC_chr_names[as.character(names(CCs[CCs %in% unique(V(G)$color[V(G)$color!="#FFFFFF"])]))])
    legend('top', pch = c(0,1), pt.bg = 'white', legend = gsub('Quibdo', 'Quibdó',
                                                               gsub('Tado', 'Tadó', 
                                                                    c(s1,s2))),
           bty = 'n', inset = -0.14)
    axis(side = 1, labels = years1stjan, las = 2, cex.axis = 0.7, at = yr01 * 2 - 1,
         line = -1.2, tick = F) # Years
  }
}

if(PDF){dev.off()}



#===============================================
# Plot without dates (not filtered)
#===============================================
if(PDF){pdf('../Plots/Graphs.pdf', height = 5, width = 10)}
par(mfrow = c(1,1), mar = c(0,4,0,0), family = 'serif', bg = 'white', pty = 'm')
# Re-create matrix A_related since need it based on unfiltered
A_related = construct_adj_matrix(Result = All_results$Unfiltered, Entry = 'rhat')
A_related_lci = construct_adj_matrix(Result = All_results$Unfiltered, Entry = 'r2.5.')
A_related[A_related_lci < eps] = 0 # Edit s.t. only not stat diff from related have weight

# Create graphs 
G_related = graph_from_adjacency_matrix(A_related, mode='upper', diag=F, weighted=T)
V(G_related)$site = SNPData[V(G_related)$name, 'City']

# ---------------------------------------------------------------------
# Aside: note that there are two related components across the entire graph
# The second consists of the DD2 clonal component plus a loosely related 
# sample, "sid95". This is the outlier in the Buenaventura Tumaco plot 
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
         legend = CC_chr_names[as.character(names(CCs[CCs %in% unique(V(G)$color[V(G)$color!="#FFFFFF"])]))])
  legend('topleft', pch = c(0,1), pt.bg = 'white', legend = gsub('Quibdo', 'Quibdó',
                                                                 gsub('Tado', 'Tadó', 
                                                                      c(s1,s2))), 
         bty = 'n', inset = 0.7)
}

if(PDF){dev.off()}


#===========================================================
# Looking specifically at clonal components found across 
# Buenaventura and Tumaco (some of also found in Guapi)
# Edge transparency proportion to rhat if statistically
# distinguishable from zero (since based on Comp_G)
#===========================================================
par(bg = 'white')
# Reconstruct graph s.t. unrelated not zero
A_related = construct_adj_matrix(All_results$Unfiltered, Entry = 'rhat')
G_related = graph_from_adjacency_matrix(A_related, mode='upper', diag=F, weighted=T)
V(G_related)$site = SNPData[V(G_related)$name, 'City']

G_B = induced_subgraph(G_related, vids = which(grepl("Buenaventura", V(G_related)$site)))
G_T = induced_subgraph(G_related, vids = which(grepl("Tumaco", V(G_related)$site)))

# Some NAs due to only components with two or more samples having names
BT_CCs = intersect(CC_chr_names[as.character(M_high[V(G_B)$name])], CC_chr_names[as.character(M_high[V(G_T)$name])])
BT_CCs = BT_CCs[!is.na(BT_CCs)]
unique_M_BT = as.numeric(names(CC_chr_names[CC_chr_names %in% BT_CCs]))

# Extract all those with said membership (includes some from guapi)
G_clonal_comp_BT = induced_subgraph(G_related, vids = which(M_high %in% unique_M_BT)) 

# Plot all those with said membership
set.seed(2)
plot.igraph(G_clonal_comp_BT, vertex.size = 3, vertex.label = NA, 
            vertex.color = M_cols[V(G_clonal_comp_BT)$name], 
            edge.color = sapply(E(G_clonal_comp_BT)$weight, adjustcolor, col = 'black')) # Make transparency a function of rhat)

# Extract edge vertex names as matrix
EV_names = do.call(rbind, strsplit(attributes(E(G_clonal_comp_BT))$vname, split = "\\|"))

# Extract average relatedness between clonal components (includes some guapi)
Mean_sd_matrix = array(NA, dim = c(length(unique_M_BT), length(unique_M_BT), 6), 
                       dimnames = list(CC_chr_names[as.character(unique_M_BT)], CC_chr_names[as.character(unique_M_BT)], 
                                       c('r mean','r 2.5%','r 97.5%',
                                         'k mean','k 2.5%','k 97.5%')))
for(j in 1:(length(unique_M_BT)-1)){
  for(i in (j+1):length(unique_M_BT)){
    members_i = names(which(M_high == unique_M_BT[i])) 
    members_j = names(which(M_high == unique_M_BT[j])) 
    inds = (EV_names[,1] %in% members_i & EV_names[,2] %in% members_j) | EV_names[,1] %in% members_j & EV_names[,2] %in% members_i 
    Mean_sd_matrix[i,j,'r mean'] = mean(E(G_clonal_comp_BT)$weight[inds])
    
    # Extract the mean lower confidence interval 
    inds_mle = (as.character(All_results$Unfiltered$individual1) %in% members_i & as.character(All_results$Unfiltered$individual2) %in% members_j) | 
      as.character(All_results$Unfiltered$individual1) %in% members_j & as.character(All_results$Unfiltered$individual2) %in% members_i 
    if (Mean_sd_matrix[i,j,'r mean'] == mean(All_results$Unfiltered$rhat[inds_mle])){ # check
      Mean_sd_matrix[i,j,'r 2.5%'] = max(All_results$Unfiltered$`r2.5.`[inds_mle]) # To see highest lower bound
      Mean_sd_matrix[i,j,'r 97.5%'] = max(All_results$Unfiltered$`r97.5.`[inds_mle]) 
      Mean_sd_matrix[i,j,'k mean'] = max(All_results$Unfiltered$khat[inds_mle])
      Mean_sd_matrix[i,j,'k 2.5%'] = max(All_results$Unfiltered$`k2.5.`[inds_mle]) # To see highest lower bound
      Mean_sd_matrix[i,j,'k 97.5%'] = max(All_results$Unfiltered$`k97.5.`[inds_mle]) 
    } 
  }
}

# Order matrices by CC name before printing
order_by_CC_name = sort.int(as.numeric(gsub('CC', '',CC_chr_names[as.character(unique_M_BT)])), index.return = T)$ix

# Print to screen
round(Mean_sd_matrix[order_by_CC_name,order_by_CC_name,], 3)

# Print for overleaf
require(kableExtra)
kable(round(Mean_sd_matrix[order_by_CC_name,order_by_CC_name,'r mean'],3), format = 'latex') 
round(Mean_sd_matrix[order_by_CC_name,order_by_CC_name, 'r 2.5%'], 3)

# Aside: relatedness between CC20, CC15 and CC14: valid example of recombination? Likely via some intermediates
earliest_sids_CC20_15_14 <- sids_1stperCC_ordered[c(20,15,14)] # For c("CC20", "CC15", "CC14")
G_earliest_sids_CC20_15_14 <- induced_subgraph(G_related, vids = earliest_sids_CC20_15_14)
plot(G_earliest_sids_CC20_15_14, 
     edge.width = E(G_earliest_sids_CC20_15_14)$weight*2, 
     vertex.label = CC_chr_names[as.character(M_high[earliest_sids_CC20_15_14])])
round(E(G_earliest_sids_CC20_15_14)$weight,2)
