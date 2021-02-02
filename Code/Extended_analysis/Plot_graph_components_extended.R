##############################################################
# Script to plot generate graphs for extended analysis
# 
# To-do list
# - Other allele frequencies
# - Cex value of heatmap legends
# - NA value of heatmap legends
# - Add inter-relatedess & collection date of first (as annotation)
# - Another graph but with samples instead of ccs_extended
# - Re-assess cliqual clonal components original Taylor et al. 
# - Update Angela's Ven diagram (e.g. enumerate number of new samples that belong to old CCs)
#
# Action points:
# - modularise "Clarifying relatedmess"
# - Generate clonal components all together
# - Fraction highly related (Ecuador Colombia connection)
# - Name samples removed
# - Sort under/over flow
##############################################################
rm(list = ls())
library(igraph) # To make graph, construct components etc.
library(RColorBrewer)
source('../igraph_functions.R') # construct_adj_matrix
eps <- 0.01 # threshold below which LCI is considered zero in Taylor et al. 2019
source("../Extended_analysis/get_vertex_clique_matrix.R")

# Load relatedness results
freqs_used <- "Taylor2020"
load(sprintf('../../RData/mles_CIs_extended_freqs%s_meta.RData', freqs_used)) 

# Load metadata
load(sprintf('../../RData/metadata_extended.RData'))

# Define city cols
cities <- unique(metadata$City)
cols_cities <- array(c(rev(brewer.pal(5, 'Spectral')), 
                       brewer.pal(length(cities)-5, 'Dark2')), 
                     dimnames = list(cities))

# Extract sids ordered by date
DE_sid <- dplyr::arrange(metadata[metadata$PloSGen2020, ], COLLECTION.DATE)$SAMPLE.CODE
FSVC_sid <- dplyr::arrange(metadata[!metadata$PloSGen2020, ], COLLECTION.DATE)$SAMPLE.CODE
DE_pair_ind <- (mle_CIs$individual1 %in% DE_sid) & (mle_CIs$individual2 %in% DE_sid)
FSVC_pair_ind <- (mle_CIs$individual1 %in% FSVC_sid) & (mle_CIs$individual2 %in% FSVC_sid)

# Using definition of clone where UCI "touches" one and LCI > 0.5
r_est = mle_CIs[FSVC_pair_ind,]
A_low <- construct_adj_matrix(r_est, Entry = 'r2.5.')
A_high <- construct_adj_matrix(r_est, Entry = 'r97.5.')
A_clonal <- array(NA, dim = dim(A_high), dimnames = dimnames(A_high))
A_clonal[(A_high > (1-eps) & A_low > 0.5)] <- 1 
A_clonal[!(A_high > (1-eps) & A_low > 0.5)] <- 0 
G_clonal <- graph_from_adjacency_matrix(A_clonal, mode='upper', diag=F) 

# Add metadata
V(G_clonal)$marker_count <- metadata[V(G_clonal)$name,  "snp_count"]
V(G_clonal)$city <- metadata[V(G_clonal)$name, "City"]
V(G_clonal)$color <- cols_cities[V(G_clonal)$city] 

# Plot graph 
plot(G_clonal, vertex.label = NA, vertex.size = 3)
articulation.points(G_clonal)

# Extract number of cliques that each sample belongs to and sort
clique_count_per_sample <- colSums(get_vertex_clique_matrix(G_clonal))

# Make symetric and plot to check
A_clonal[lower.tri(A_clonal)] <- t(A_clonal)[lower.tri(A_clonal)]
diag(A_clonal) <- 1
image(A_clonal)  

# Compute number of NA comparisons per sample
na_count_per_sample <- rowSums(is.na(A_clonal))

par(mfrow = c(2,2))

# Plot clique versus marker count for inputG
plot(y = clique_count_per_sample, 
     x = metadata[names(clique_count_per_sample), "snp_count"], 
     ylab = "Per-sample clique membership count", 
     xlab = "Per-sample marker data count",
     bty = "n", pch = 20, cex.main = 1) 

# Plot clique versus marker count for inputG
plot(y = clique_count_per_sample, 
     x = na_count_per_sample[names(clique_count_per_sample)], 
     ylab = "Per-sample clique membership count", 
     xlab = "Per-sample NA relatedness count",
     bty = "n", pch = 20, cex.main = 1) 

# Plot clique versus marker count for inputG
plot(y = metadata[names(na_count_per_sample), "snp_count"], 
     x = na_count_per_sample, 
     xlab = "Per-sample marker data count", 
     ylab = "Per-sample NA relatedness count",
     bty = "n", pch = 20, cex.main = 1) 

#-----------------------------------------------------------
# Aside: could use community detection 
# as alternative to more strictly defined components and cliques
#-----------------------------------------------------------
wc <- cluster_walktrap(G_clonal) 
plot(wc, G_clonal, vertex.label = NA, vertex.size = 3, vetex.color = V(G_clonal)$color)


#-------------------------------------------------------
# Decided upon subgraph
# Check sids_removed with Angela and Manuela
#-------------------------------------------------------
snp_count_threshold <- 42 # Based on plots above
sids_to_keep <- names(clique_count_per_sample)[!metadata[names(clique_count_per_sample), "snp_count"] <= snp_count_threshold]
sids_removed <- setdiff(names(clique_count_per_sample), sids_to_keep)
subG <- subgraph(G_clonal, sids_to_keep)
art_pt <- articulation.points(subG)
CsubG <- components(subG)

Clonal_components <- lapply(1:length(CsubG$csize), function(cc){
  names(which(CsubG$membership == as.numeric(cc)))
})

# Extract the sid of the earliest detected sample per cc
sids_1stperCC <- t(sapply(Clonal_components, function(cc){
  if (all(is.na(metadata[cc,"COLLECTION.DATE"]))) {
    earliest_date <- min(metadata[cc,"Year"])
    sids <- cc[which(metadata[cc,"Year"] == earliest_date)]
  } else {
    date <- min(metadata[cc,"COLLECTION.DATE"], na.rm = T)
    year <- format(as.Date(min(metadata[cc,"Year"]), format = "%Y"), "%Y")
    if (year >= format(date, "%Y")) {
      earliest_date <- date
      sids <- cc[which(metadata[cc,"COLLECTION.DATE"] == earliest_date)]
    } else {
      earliest_date <- as.numeric(as.character(year))
      sids <- cc[which(metadata[cc,"Year"] == earliest_date)]
    }
  }
  sid <- sids[1] # Take the first if there is more than one
}))

# Extract metadata, ordered by date, for the each of the earliest detected sample per cc
sids_1stperCC_metadata <- dplyr::arrange(metadata[sids_1stperCC,], Year, COLLECTION.DATE)

# Create "CC#" character name order by the earliest detected sample per CC
CC_chr_names <- paste0("cc_", 1:length(CsubG$csize)) # Name differently 
names(CC_chr_names) = as.character(CsubG$membership[sids_1stperCC_metadata$SAMPLE.CODE]) # Ensure CC names are ordered as memberships

# Rename and reorder Clonal_components
names(Clonal_components) <- 1:length(Clonal_components)
names(Clonal_components) <- CC_chr_names[names(Clonal_components)]
Clonal_components <- Clonal_components[CC_chr_names]

# save
save(Clonal_components, 
     file = "../../RData/Clonal_components_extended_FSVC_sid.RData")

dev.off()



