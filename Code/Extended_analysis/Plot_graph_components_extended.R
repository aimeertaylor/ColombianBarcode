##############################################################
# Script to plot igraph results using cliques to define components. 
# 
# To-do list
# - Other allele frequencies
# - Convert files into markdown (format_mles_CIs..., this script, Generate_and_plot_relatedness.. )
# - Cex value of heatmap legends
# - NA value of heatmap legends
# - What is the relatedness with Guyana samples? 
# - Why those cluster names? 
# - Dive into Ecuador Colombia connection
# - Add inter-relatedess & collection date of first (as annotation)
# - Another graph but with samples instead of ccs_extended
# - Re-assess cliqual clonal components original
# - Update Angela's Ven diagram (e.g. enumerate number of new samples that belong to old CCs)
# 
#' Doesn't comply to transitivity 
#' Can recover CCs by progressively removing highly related samples 
#' Question: which clonal definition makes most sense? 
#' Question: where to stop with the filtering of samples? 
#' Question: are these samples the "low quality" ones? 
#' Does cc_7 cluster with CCXX when we progressively look for structure among all
#' Re-add clique analysis for the Taylor et al. cliques
##############################################################
rm(list = ls())
library(igraph) # To make graph, construct components etc.
library(RColorBrewer) # For colours
source('../igraph_functions.R') # construct_adj_matrix
eps <- 0.01 # threshold below which LCI is considered zero in Taylor et al. 2019

# Load relatedness results
freqs_used <- "Taylor2020"
load(sprintf('../../RData/mles_CIs_extended_freqs%s_meta.RData', freqs_used)) 

# Load metadata
load(sprintf('../../RData/metadata_extended.RData'))

# Extract city data
cities <- unique(c(mle_CIs$City1, mle_CIs$City2))
cols_cities <-c(rev(brewer.pal(5, 'Spectral')),
                1:(length(cities)-5))#brewer.pal(length(cities)-5, 'PiYG'))
names(cols_cities) <- as.character(cities)

# Extract sids ordered by date
FSVC_sid <- dplyr::arrange(metadata[!metadata$PloSGen2020, ], COLLECTION.DATE)$SAMPLE.CODE
FSVC_pair_ind <- (mle_CIs$individual1 %in% FSVC_sid) & (mle_CIs$individual2 %in% FSVC_sid)


#----------------------------------------------------------------
# Matrix and graph of all samples inc. those in Taylor et al. 2020
#----------------------------------------------------------------
# Using definition of clone where UCI "touches" one and LCI > 0.5
# (this definition is based on plots of relatedness estimates)
A_low <- construct_adj_matrix(mle_CIs, Entry = 'r2.5.')
A_high <- construct_adj_matrix(mle_CIs, Entry = 'r97.5.')  
A_high[!(A_high >= (1-eps) & A_low > 0.5)] <- 0 # Edit s.t. only clonal have weight
A_high[(A_high >= (1-eps) & A_low > 0.5)] <- 1 # Edit s.t. only clonal have weight
image(A_high)
G_clonal <- graph_from_adjacency_matrix(A_high, mode='upper', diag=F, weighted=T) 
V(G_clonal)$color <- cols_cities[metadata[attr(V(G_clonal), "names"), "City"]]
plot(G_clonal, vetex.color = V(G_clonal)$color, vertex.label = NA, vertex.size = 3)
cities_ <- unique(metadata[attr(V(G_clonal), "names"), "City"])
legend('topleft', pch = 16, bty = 'n', cex = 0.5, pt.cex = 1, 
       col = cols_cities[cities_], legend = cities_, inset = -0.1)
writeLines(sprintf("%s clonal component/s with two or more samples", 
                   sum(components(G_clonal)$csize > 1))) 



# ----------------------------------------------------------------
# Matrix of samples not in Taylor et al. 2020
# Using definition of clone where UCI "touches" one and LCI > 0.5
# (this definition is based on plots of relatedness estimates)
# ----------------------------------------------------------------
A_low <- construct_adj_matrix(mle_CIs[FSVC_pair_ind,], Entry = 'r2.5.')
A_high <- construct_adj_matrix(mle_CIs[FSVC_pair_ind,], Entry = 'r97.5.') 
A_high[!(A_high >= (1-eps) & A_low > 0.5)] <- 0 # Edit s.t. only clonal have weight
A_high[(A_high >= (1-eps) & A_low > 0.5)] <- 1 # Edit s.t. only clonal have weight
image(A_high)


# ----------------------------------------------------------------
# Matrix of samples not in Taylor et al. 2020
# Using definition of clone where rhat > 0.99 
# (this definition also seems reasonable given relatedness estimate plots)
# ----------------------------------------------------------------
A_est <- construct_adj_matrix(mle_CIs[FSVC_pair_ind,], Entry = 'rhat') 
A_est[A_est < 1-eps] <- 0 # Edit s.t. only clonal have weight
A_est[A_est >= (1-eps)] <- 1 # Edit s.t. only clonal have weight


# ----------------------------------------------------------------
# How different are matrices defined using rhat versus CI?
# ----------------------------------------------------------------
identical(A_est,A_high) # Not the same
fields::image.plot(A_est!=A_high, col = c("gray", "black"))
sum(A_est!=A_high, na.rm = T) # Only 153 differences
sum(A_est,na.rm=T) # More conservative 
sum(A_high,na.rm=T) # Classifies more comparisons as clonal


# ----------------------------------------------------------------
# Graph using A_est (might change to A_high later)
# One big hairball
# ----------------------------------------------------------------
G_clonal <- graph_from_adjacency_matrix(A_est, mode='upper', diag=F, weighted=T) 
V(G_clonal)$color <- cols_cities[metadata[attr(V(G_clonal), "names"), "City"]]
plot(G_clonal, vetex.color = V(G_clonal)$color, vertex.label = NA, vertex.size = 3)
cities_ <- unique(metadata[attr(V(G_clonal), "names"), "City"])
legend('topleft', pch = 16, bty = 'n', cex = 0.5, pt.cex = 1, 
       col = cols_cities[cities_], legend = cities_, inset = -0.1)
writeLines(sprintf("%s clonal component/s with two or more samples", 
                   sum(components(G_clonal)$csize > 1))) 



#=============================================================
# Exploration of structure with G_clonal
#=============================================================
writeLines("Samples whose removal increases the number of graph components:") 
articulation_points(G_clonal)

# Cliques
Cliques <- max_cliques(G_clonal, min = 2) 
clique_no <- length(Cliques)
writeLines(sprintf("Number of clonal cliques: %s", clique_no))

# Clique overlap
clique_names <- paste0("Clique", 1:clique_no)
names(Cliques) <- clique_names
Clique_intersects <- list()
for(i in 1:(clique_no-1)){
  for(j in (i+1):clique_no){
    i_ <- clique_names[i]
    j_ <- clique_names[j]
    x <- intersect(names(Cliques[[i_]]),
                   names(Cliques[[j_]])) # names returns characters (otherwise returns vertex count)
    if(length(x) > 0) {
      Clique_intersects[[i_]][[j_]] <- x 
    } else {
      next()
    }
  }
}

# Clique per sample
Clique_sample_matrix <- sapply(unique(unlist(sapply(Cliques, names))), function(sid){
  sapply(Cliques, function(Clique){
    as.character(sid) %in% names(Clique)
  })
})


# Sort intersecting samples by the number of cliques they belong to
clique_count_per_sample <- colSums(Clique_sample_matrix)

# Are the highly connected points those with less data? 
plot(y = as.numeric(clique_count_per_sample), 
     x = metadata[names(clique_count_per_sample), "snp_count"], pch = 20, 
     ylab = "Number of clique memberships per sample", 
     xlab = "Available data per sample (SNPs)")
abline(h = 30)

# Zoom in: 
plot(y = as.numeric(clique_count_per_sample), 
     x = metadata[names(clique_count_per_sample), "snp_count"], pch = 20, 
     ylab = "Number of clique memberships per sample", 
     xlab = "Available data per sample (SNPs)", 
     xlim = c(0,100))
abline(h = 30)

# Progressively removing samples that have few SNPs and are highly connected
# each time using a more stringent snp_count_threshold shows: 
# 1) the emergence and then dissappearance of clonal components
# 2) the emergence and then disappearence of attrition points
# I consider the optimum to be highest clonal components without attrition points
# since, thereafter, one is just removing samples from within the clonal components
snp_count_thresholds <- seq(0,100,1) # Inspect limits by hand
art_pt_counts <- array(NA, dim = length(snp_count_thresholds), dimnames = list(snp_count_thresholds))
cc_counts <- array(NA, dim = length(snp_count_thresholds), dimnames = list(snp_count_thresholds))
  
for(snp_count_threshold in snp_count_thresholds){ 
  
  sids_to_keep <- names(clique_count_per_sample)[!metadata[names(clique_count_per_sample), "snp_count"] <= snp_count_threshold]
  subG <- subgraph(G_clonal, sids_to_keep)
  art_pt <- articulation.points(subG)
  cc_count <- sum(components(subG)$csize > 1)
  cc_counts[as.character(snp_count_threshold)] <- cc_count
  art_pt_counts[as.character(snp_count_threshold)] <- length(art_pt)
  
  # Visualise results every 5th increment
  if ((snp_count_threshold %% 5)==0) {
    print(snp_count_threshold)
    writeLines(sprintf("%s clonal component/s with two or more samples", cc_count)) 
    print(cbind(art_pt, snp_count = metadata[art_pt, "snp_count"]))
    plot(subG, vertex.label = NA, vertex.size = 3)
  }
}

#-------------------------------------------------------
# Use plots of clonal component counts and articulation 
# points to select the best snp_count_threshold for filtering 
#-------------------------------------------------------
par(mfrow = c(2,1))
plot(y = cc_counts, x = snp_count_thresholds, type = "l")
abline(v = as.numeric(names(which.max(cc_counts))))
plot(y = art_pt_counts, x = snp_count_thresholds, type = "l")

#-----------------------------------------------------------
# Aside: could use community detection 
# as alternative to more strictly defined components and cliques
#-----------------------------------------------------------
wc <- cluster_walktrap(subG) 
plot(wc, subG, vertex.label = NA, vertex.size = 3, vetex.color = V(G_clonal)$color)


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





