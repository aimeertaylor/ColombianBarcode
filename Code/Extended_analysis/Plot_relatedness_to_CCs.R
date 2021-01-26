##############################################################
# Script to plot matrices of relatedness between 
# "new" samples and "old" CCs (CCs of Taylor et al.)
# Need to plot heatmaps using ggplot to specify NA values
##############################################################
rm(list = ls())
library(RColorBrewer) # For colours
library(igraph)
eps <- 0.01 # threshold below which LCI is considered zero in Taylor et al. 2019
PDF <- T

# Load metadata
load(sprintf('../../RData/metadata_extended.RData'))
load(file = "../../RData/relatedness_to_CCs.RData")
load("../../RData/Clonal_components_extended_FSVC_sid.RData")
CC_extended <- Clonal_components
load("../../RData/Clonal_components.RData")
CC_original <- Clonal_components

# CC counts and sizes
cc_original_sizes <- sapply(CC_original, length)
cc_extended_sizes <- sapply(CC_extended, length)
cc_no_original <- length(cc_original_sizes)
cc_no_extended <- length(cc_extended_sizes)

# Concatenate all CCs
CC_all <- c(CC_original, CC_extended)
CC_sizes <- c(cc_original_sizes,cc_extended_sizes)

# Extract city data
cities <- unique(metadata$City)
cols_cities <-c(rev(brewer.pal(5, 'Spectral')),
                brewer.pal(length(cities)-5, 'Dark2'))
names(cols_cities) <- as.character(cities)

# Sample counts per city for pie charts
sample_count_per_city_per_CC <- array(0, dim = c(length(CC_all),
                                                 length(cities)), 
                                      dimnames = list(names(CC_all), cities))

sample_count_per_city_per_CC <- lapply(CC_all, function(cc){
  # Initialise empty vector of city counts
  city_count_per_CC_inc_zeros = array(0, dim = length(cities), dimnames = list(cities))
  # Extract number of samples per CC per city
  city_count_per_CC_exc_zeros <-  table(metadata[cc,"City"])
  # Populate initial empty vector of city counts
  city_count_per_CC_inc_zeros[names(city_count_per_CC_exc_zeros)] <- city_count_per_CC_exc_zeros
  return(city_count_per_CC_inc_zeros)
})

# Average relatedness between CCs
CCs_av_rhat <- sapply(cc_extended_to_cc_original, function(x) x["av_rhat", ])

# Re-order extended clonal components by relatedness 
# (otherwise ordered by date first detected)
CCs_extended_names_ordered <- names(sort(colSums(CCs_av_rhat)))
cc_extended_to_cc_original <- cc_extended_to_cc_original[CCs_extended_names_ordered]
CCs_av_rhat <- CCs_av_rhat[,names(sort(colSums(CCs_av_rhat)))]

# Order sids by cc followed by samples that were removed from extended clonal components
FSVC_sid <- names(sid_extended_to_cc_original)
FSVC_sid_ordered <- c(unlist(CC_extended[CCs_extended_names_ordered]), 
                      FSVC_sid[!FSVC_sid %in% unlist(CC_extended)])
names(FSVC_sid_ordered) <- FSVC_sid_ordered
sid_extended_to_cc_original <- sid_extended_to_cc_original[FSVC_sid_ordered]

if (PDF) pdf("../../Plots/Relatedness_to_CCs.pdf", )

#================================================
# Heatplots between CCs
#================================================
for(x in c("av_rhat", "av_r2.5", "av_r97.5", "mn_r2.5", "mx_r97.5")){
  to_plot <- sapply(cc_extended_to_cc_original, function(y) y[x, ])
  fields::image.plot(to_plot, xaxt = "n", yaxt = "n", main = x,
                     breaks = sort(c(eps, seq(0,1,length.out = 9), 1-eps)),  
                     col = c("black",grey.colors(9)), border = "white", 
                     lab.breaks = sort(c(eps, seq(0,1,length.out = 9), 1-eps)))
  axis(side = 1, at = seq(0,1,length.out = nrow(to_plot)), 
       labels = rownames(to_plot), cex.axis = 0.5, las = 2)
  axis(side = 2, at = seq(0,1,length.out = ncol(to_plot)), 
       labels = colnames(to_plot), cex.axis = 0.5, las = 1)
}

#================================================
# Heatplots between sids and CCs
# Order clonal components and then
# FSVC_sids by clonal components and then add
# the samples that were removed due to high 
# connectedness but low data
#================================================
for(x in c("av_rhat", "av_r2.5", "av_r97.5", "mn_r2.5", "mx_r97.5")){
  to_plot <- sapply(sid_extended_to_cc_original, function(y) y[x, ])
  fields::image.plot(to_plot, xaxt = "n", yaxt = "n", main = x,
                     breaks = sort(c(eps, seq(0,1,length.out = 9), 1-eps)),  
                     col = c("black",grey.colors(9)), 
                     lab.breaks = sort(c(eps, seq(0,1,length.out = 9), 1-eps)))
  axis(side = 1, at = seq(0,1,length.out = nrow(to_plot)), 
       labels = rownames(to_plot), cex.axis = 0.5, las = 2)
  axis(side = 2, at = seq(0,1,length.out = ncol(to_plot)), 
       labels = colnames(to_plot), cex.axis = 0.25, las = 1)
}



#============================================
# Graph CCs 
# Different ways of accessing graph aspects
# attributes(E(BiG), "weight")
# vertex_attr(BiG, "name")
# attr(V(BiG), "name")
#============================================
BiG <- graph_from_incidence_matrix(CCs_av_rhat, weighted = TRUE)
vertex_sizes =  CC_sizes[attr(V(BiG), "names")]


# Space vertices according to size
vertex_spacing_original <- rep(1,cc_no_original)
for(i in 2:cc_no_original){
  vertex_spacing_original[i] <- vertex_spacing_original[i-1] + 
    (vertex_sizes[rownames(CCs_av_rhat)][i-1])/2  + (vertex_sizes[rownames(CCs_av_rhat)][i])/2 
}

vertex_spacing_extended <- rep(1,cc_no_extended)
for(i in 2:cc_no_extended){
  vertex_spacing_extended[i] <- vertex_spacing_extended[i-1] + 
    (vertex_sizes[colnames(CCs_av_rhat)][i-1])/2  + (vertex_sizes[colnames(CCs_av_rhat)][i])/2 
}

# Name
names(vertex_spacing_original) <- rownames(CCs_av_rhat)
names(vertex_spacing_extended) <- colnames(CCs_av_rhat)

# Transform to -1,1 igraph range; see locator(2)
# ( (old_value - old_min) / (old_max - old_min) ) * (new_max - new_min) + new_min
vertex_spacing_original <- ((vertex_spacing_original-min(vertex_spacing_original))/diff(range(vertex_spacing_original)))*2-1
vertex_spacing_extended <- ((vertex_spacing_extended-min(vertex_spacing_extended))/diff(range(vertex_spacing_extended)))*2-1

# By default the rows are drawn first in a bipartite graph
vertex_spacing <- c(vertex_spacing_original, vertex_spacing_extended)
range(vertex_spacing) 

my_bi_partite_layout <- cbind(rep(0:1, c(nrow(CCs_av_rhat), 
                                         ncol(CCs_av_rhat))),
                              vertex_spacing[V(BiG)$name])

#======================================
# Plot
#' CCs original (left), plotted in order of first (bottom) to last (top) detected. 
#' CCs extended (right), plotted in order of average relatedness to CCs original from 
#' least related (bottom) to most related (top).
#======================================
par(mar = c(1,1,1,1))
writeLines("Average relatedness estimates range from and to:")
round(range(edge_attr(BiG, "weight")), 3)

plot(BiG,
     layout = my_bi_partite_layout, 
     vertex.shape = "pie", 
     vertex.pie = sample_count_per_city_per_CC[V(BiG)$name], 
     vertex.pie.color = list(cols_cities[cities]), 
     vertex.pie.lwd = 0.25, 
     vertex.frame.color = NA, #NB, colours pi outline 
     vertex.size = vertex_sizes, 
     # vertex.pie.lty = "dashed", 
     # vertex.pie.border = list("hotpink"), # Doesn't work: always black
     # vertex.label.cex = 0.35, 
     # vertex.label.color = 'black', 
     vertex.label = NA,
     edge.width = edge_attr(BiG, "weight"),
     edge.color = sapply(edge_attr(BiG, "weight"), 
                         function(x)adjustcolor('black', alpha.f = x)))

# # Add CC outline
# plot(BiG, add = T,
#      layout = my_bi_partite_layout,
#      vertex.size = vertex_sizes,
#      vertex.frame.color = "white", #NB, coulours pi outline
#      vertex.frame.lwd = 0.5,
#      vertex.color = NA,
#      vertex.label = NA,
#      edge.color = NA)

axis(side = 2, at = vertex_spacing_original, 
     labels = names(vertex_spacing_original), line = -2, 
     las = 1, cex.axis = 0.25, tick = F)

axis(side = 4, at = vertex_spacing_extended, 
     labels = names(vertex_spacing_extended), line = -2, 
     las = 1, cex.axis = 0.25, tick = F)

# Legend
city_counts <- sapply(sample_count_per_city_per_CC[V(BiG)$name], function(x){x})
cities_ <- names(which(rowSums(city_counts) > 0))
legend('bottom', pch = 16, bty = 'n', cex = 0.5, pt.cex = 1, 
       col = cols_cities[cities_], legend = cities_)

if (PDF) dev.off()