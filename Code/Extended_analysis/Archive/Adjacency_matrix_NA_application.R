#################################################################################
#' Obsolte, since graph_from_adjacency_matrix understood
#'
#' Script to demonstrate the impact of igraph's treatment of NAs on the graphs
#' generated using the extended relatedness estimates. Compare the number of
#' cliques each sample belongs to to the number of NA comparisons each sample
#' has. If G_clonal was created using graph_from_adjacency_matrix with weighted
#' = NULL (default), we expect fewer than average cliques among samples with
#' many NA comparisons. Otherwise, if G_clonal was created with weighted = T, we
#' expect more than average clique counts among samples with many NA
#' comparisons.
#################################################################################
rm(list = ls())
library(igraph) 
source('../igraph_functions.R') # construct_adj_matrix
source("../Extended_analysis/get_vertex_clique_matrix.R")
eps <- 0.01 # threshold below which LCI is considered zero in Taylor et al. 2019

# Load relatedness results and metadata
load('../../RData/mles_CIs_extended_freqsTaylor2020_meta.RData') 
load('../../RData/metadata_extended.RData')

# Extract sids ordered by date
FSVC_sid <- dplyr::arrange(metadata[!metadata$PloSGen2020, ], COLLECTION.DATE)$SAMPLE.CODE
FSVC_pair_ind <- (mle_CIs$individual1 %in% FSVC_sid) & (mle_CIs$individual2 %in% FSVC_sid)

# Using definition of clone where UCI "touches" one and LCI > threshold
r_est = mle_CIs[FSVC_pair_ind,]
A_low <- construct_adj_matrix(r_est, Entry = 'r2.5.')
A_high <- construct_adj_matrix(r_est, Entry = 'r97.5.')
A_clonal <- array(NA, dim = dim(A_high), dimnames = dimnames(A_high))
A_clonal[(A_high > (1-eps) & A_low > 0.95)] <- 1
A_clonal[!(A_high > (1-eps) & A_low > 0.95)] <- 0

# Creat graph and add metadata
G_clonal <- graph_from_adjacency_matrix(A_clonal, mode='upper', diag=F, 
                                        weighted = T) 

# Extract number of cliques that each sample belongs to and sort
clique_count_per_sample <- colSums(get_vertex_clique_matrix(G_clonal))

# Make symetric and plot to check
A_clonal_full <- A_clonal
A_clonal_full[lower.tri(A_clonal_full)] <- t(A_clonal_full)[lower.tri(A_clonal_full)]
diag(A_clonal_full) <- 1
image(A_clonal)
image(A_clonal_full) 

# Compute number of NA comparisons per sample
na_count_per_sample <- rowSums(is.na(A_clonal_full))

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




