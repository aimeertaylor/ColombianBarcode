#################################################################################
#' Script to generate and plot a graph of relatedness between samples
#################################################################################
rm(list = ls())
library(igraph) 
library(RColorBrewer)
source('../igraph_functions.R') # for construct_adj_matrix
eps <- 0.01 # threshold below which LCI is considered zero in Taylor et al. 2019
a = 0; b = 0.5; c = 0.25 # Scaling parameters for plotting edge width and transparancy
PDF <- TRUE

# Load relatedness results and metadata
freqs_used <- "Taylor2020"
load(sprintf('../../RData/mles_CIs_extended_freqs%s_meta.RData', freqs_used)) 
load(sprintf('../../RData/metadata_extended.RData'))
load("../../RData/sids_to_remove_from_graphs.RData")

# Define city cols
cities <- unique(metadata$City)
cols_cities <- array(c(rev(brewer.pal(5, 'Spectral')), 
                       brewer.pal(length(cities)-5, 'Dark2')), 
                     dimnames = list(cities))

# Remove samples with one or more NA estimate
sids_to_keep <- setdiff(metadata$SAMPLE.CODE, sids_remv) 
inds_to_keep <- mle_CIs$individual1 %in% sids_to_keep & mle_CIs$individual2 %in% sids_to_keep
mles_to_keep <- mle_CIs[inds_to_keep,]

A_est <- construct_adj_matrix(mles_to_keep, Entry = 'rhat')

if(PDF) pdf("../../Plots/Relatedness_graphs.pdf")

for(high_relatedness_threshold in c(0,0.25)){
  A_est[A_est < high_relatedness_threshold] <- 0
  fields::image.plot(A_est)
  
  # Create graph 
  G_est <- graph_from_adjacency_matrix(A_est, mode='upper', diag=F, weighted = T) 
  art_pnts <- names(articulation.points(G_est))
  ccompnts <- components(G_est)
  
  # Add metadata
  V(G_est)$marker_count <- metadata[V(G_est)$name,  "snp_count"]
  V(G_est)$city <- metadata[V(G_est)$name, "City"]
  V(G_est)$color <- cols_cities[V(G_est)$city] 
  
  weights_rescaled = a + (b-a) * (E(G_est)$weight - min(E(G_est)$weight)) / (max(E(G_est)$weight) - min(E(G_est)$weight))
  
  # Plot graph 
  plot(G_est, 
       main = sprintf("High relatedness threshold: %s", high_relatedness_threshold), 
       cex.main = 0.75, 
       vertex.label.color = "black", 
       vertex.label.cex = 1, 
       vertex.frame.color = NA,
       vertex.label = NA, 
       vertex.shape = c("square","circle")[metadata[V(G_est)$name, "PloSGen2020"] + 1],
       vertex.size = 0.5*log(metadata[V(G_est)$name, "snp_count"]),
       edge.width = weights_rescaled * c, 
       edge.color = sapply(weights_rescaled, function(x)adjustcolor('black', alpha.f = x)))
  
  # Legend
  cities_ <- unique(V(G_est)$city)
  legend('bottomleft', pch = 16, bty = 'n', cex = 0.5, pt.cex = 1, 
         inset = -0.1, col = cols_cities[cities_], legend = cities_)
  
  # R estimate range
  range(A_est, na.rm = T)
}

if(PDF) dev.off()
