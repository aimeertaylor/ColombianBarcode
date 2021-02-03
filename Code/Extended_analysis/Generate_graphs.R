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
# - Generate clonal components all together (plot bipartite alltogether)
# - Fraction highly related (Ecuador Colombia connection)
# - Sort under/over flow
# - Update Angela and Manuela
#
#' Conclusion, stop searching for "ideal" components among new samples. -
#' Instead, use new clonal components as a tool only e.g. to organise samples on
#' heat map.
##############################################################
rm(list = ls())
library(igraph) 
library(RColorBrewer)
source('../igraph_functions.R') # construct_adj_matrix
eps <- 0.01 # threshold below which LCI is considered zero in Taylor et al. 2019
source("../Extended_analysis/get_vertex_clique_matrix.R")

# Load relatedness results and metadata
freqs_used <- "Taylor2020"
load(sprintf('../../RData/mles_CIs_extended_freqs%s_meta.RData', freqs_used)) 
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

# Using definition of clone where UCI "touches" one and LCI > threshold
A_low <- construct_adj_matrix(mle_CIs[FSVC_pair_ind,], Entry = 'r2.5.')
A_high <- construct_adj_matrix(mle_CIs[FSVC_pair_ind,], Entry = 'r97.5.')
A_est <- construct_adj_matrix(mle_CIs[FSVC_pair_ind,], Entry = 'rhat')

# Note that we either need to do the one/zero dichotomy in order to avoid
# weighted = T in graph_from_adjacency_matrix (we need to avoid weighted = T in
# graph_from_adjacency_matrix because some edges have NA weights).
# Alternatively, we can set NAs to zero and use weighted = T (see commented
# below)
# A_clonal <- A_est
# A_clonal[A_high < (1-eps) | A_low < 0.95] <- 0
# A_clonal[is.na(A_clonal)] <- 0
A_clonal <- array(NA, dim = dim(A_high), dimnames = dimnames(A_high))
A_clonal[A_high >= (1-eps) & A_low > 0.95] <- 1
A_clonal[A_high < (1-eps) | A_low < 0.95] <- 0

# Create graph and add metadata
G_clonal <- graph_from_adjacency_matrix(A_clonal, mode='upper', diag=F, weighted = T) 
V(G_clonal)$marker_count <- metadata[V(G_clonal)$name,  "snp_count"]
V(G_clonal)$city <- metadata[V(G_clonal)$name, "City"]
V(G_clonal)$color <- cols_cities[V(G_clonal)$city] 

# Plot graph 
plot(G_clonal, vertex.label = NA, 
     vertex.size = log(metadata[V(G_clonal)$name, "snp_count"]))

art_pnts <- articulation.points(G_clonal)
ccompnts <- components(G_clonal)

# Extract 
Clonal_components <- lapply(1:length(ccompnts$csize), function(cc){
  names(which(ccompnts$membership == as.numeric(cc)))
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
CC_chr_names <- paste0("cc_", 1:length(ccompnts$csize)) # Name differently
names(CC_chr_names) = as.character(ccompnts$membership[sids_1stperCC_metadata$SAMPLE.CODE]) # Ensure CC names are ordered as memberships

# Rename and reorder Clonal_components
names(Clonal_components) <- 1:length(Clonal_components)
names(Clonal_components) <- CC_chr_names[names(Clonal_components)]
Clonal_components <- Clonal_components[CC_chr_names]

# save
save(Clonal_components, file = "../../RData/Clonal_components_extended_FSVC_sid.RData")



