##############################################################
# Script to plot igraph results 
#' Remove low quality samples if CCs very sensitive or 
#' use cliques instead of components. 
#' 
#' To-do list
#' Finish plot relatedness estimates, 
#' CC heatmap (using R_comp), 
#' other allele frequencies
#' Does it make sense to plot between as agraph? 
#' Summary, doesn't comply to transitive 
#' Can more-or-less recover the transitive nature by removing 50 ish samples 
#' Could these samples be the "low quality" ones? Need to talk to Angela and Manuela
##############################################################
rm(list = ls())
library(igraph) # To make graph, construct components etc.
library(RColorBrewer) # For colours
source('../igraph_functions.R') # construct_adj_matrix
eps <- 0.01 # threshold below which LCI is considered zero in Taylor et al. 2019

# Load relatedness results
freqs_used <- "Taylor2020"
load(sprintf('../../RData/mles_CIs_extended_freqs%s_meta.RData', freqs_used)) 
mle_CIs <- mle_CIs[!is.na(mle_CIs$rhat), ] # remove NAs 
mle_CIs <- mle_CIs[!(mle_CIs$r2.5. < eps & mle_CIs$r97.5. > (1-eps)), ] # Remove uninformative

# Load metadata
load(sprintf('../../RData/metadata_extended.RData'))
rownames(metadata) <- metadata$SAMPLE.CODE

# Load SNP data to get snp counts per sample
load(file = "../../RData/snpdata_extended.RData")
sample_snp_counts <- colSums(!is.na(snpdata[,-(1:2)]))

# Extract city data
cities <- unique(c(mle_CIs$City1, mle_CIs$City2))
cols_cities <-c(rev(brewer.pal(5, 'Spectral')),
                1:(length(cities)-5))#brewer.pal(length(cities)-5, 'PiYG'))
names(cols_cities) <- as.character(cities)

# Extract sids ordered by date
DE_sid <- dplyr::arrange(metadata[metadata$PloSGen2020, ], COLLECTION.DATE)$SAMPLE.CODE
FSVC_sid <- dplyr::arrange(metadata[!metadata$PloSGen2020, ], COLLECTION.DATE)$SAMPLE.CODE
DE_pair_ind <- (mle_CIs$individual1 %in% DE_sid & mle_CIs$individual2 %in% DE_sid)
FSVC_pair_ind <- (mle_CIs$individual1 %in% FSVC_sid) & (mle_CIs$individual2 %in% FSVC_sid)

#=================================
# Some quick graphs
#=================================
# Using definition of clone where UCI "touches" one and LCI > 0.5
# (this definition is based on plots of relatedness estimates)
A_low <- construct_adj_matrix(mle_CIs, Entry = 'r2.5.')
A_high <- construct_adj_matrix(mle_CIs, Entry = 'r97.5.') # Adj. matrix unfiltered using 'r97.5.' 
A_high[!(A_high >= (1-eps) & A_low > 0.5)] <- 0 # Edit s.t. only clonal have weight
A_high[(A_high >= (1-eps) & A_low > 0.5)] <- 1 # Edit s.t. only clonal have weight
G_high <- graph_from_adjacency_matrix(A_high, mode='upper', diag=F, weighted=T) # Construct graph 
V(G_high)$color <- cols_cities[metadata[attr(V(G_high), "names"), "City"]]
plot(G_high, vetex.color = V(G_high)$color, vertex.label = NA, vertex.size = 3)
cities_ <- unique(metadata[attr(V(G_high), "names"), "City"])
legend('topleft', pch = 16, bty = 'n', cex = 0.5, pt.cex = 1, 
       col = cols_cities[cities_], legend = cities_, inset = -0.1)
writeLines(sprintf("%s clonal component/s with two or more samples", 
                   sum(components(G_high)$csize > 1))) 
image(A_high)

# Using definition of clone where UCI "touches" one as in Taylor et al. 2020
A_high <- construct_adj_matrix(mle_CIs[DE_pair_ind,], Entry = 'r97.5.') # Adj. matrix unfiltered using 'r97.5.' 
A_high[A_high < 1-eps] <- 0 # Edit s.t. only clonal have weight
A_high[A_high >= 1-eps] <- 1 # Edit s.t. only clonal have weight
G_high <- graph_from_adjacency_matrix(A_high, mode='upper', diag=F, weighted=T) # Construct graph 
V(G_high)$color <- cols_cities[metadata[attr(V(G_high), "names"), "City"]] # see ?attributes(V(G_high))
plot(G_high, vetex.color = V(G_high)$color, vertex.label = NA, vertex.size = 3)
Cliques <- max_cliques(G_high, min = 2) # Clonal component "11" splits into two cliques; see below
cities_ <- unique(metadata[attr(V(G_high), "names"), "City"])
legend('topleft', pch = 16, bty = 'n', cex = 0.5, pt.cex = 1, 
       col = cols_cities[cities_], legend = cities_, inset = -0.1)
writeLines(sprintf("%s clonal components with two or more samples", 
           sum(components(G_high)$csize > 1))) 
image(A_high)

# Using definition of clone where rhat > 0.99 
# (this definition also seems reasonable given relatedness estimate plots)
A_high <- construct_adj_matrix(mle_CIs[FSVC_pair_ind,], Entry = 'rhat') # Adj. matrix unfiltered using 'r97.5.' 
A_high[A_high < 1-eps] <- 0 # Edit s.t. only clonal have weight
A_high[A_high > (1-eps)] <- 1 # Edit s.t. only clonal have weight
G_high <- graph_from_adjacency_matrix(A_high, mode='upper', diag=F, weighted=T) # Construct graph 
V(G_high)$color <- cols_cities[metadata[attr(V(G_high), "names"), "City"]]
cities_ <- unique(metadata[attr(V(G_high), "names"), "City"])
legend('topleft', pch = 16, bty = 'n', cex = 0.5, pt.cex = 1, 
       col = cols_cities[cities_], legend = cities_, inset = -0.1)
writeLines(sprintf("%s clonal component/s with two or more samples", 
                   sum(components(G_high)$csize > 1))) 

# Using definition of clone where UCI "touches" one and LCI > 0.5
# (this definition is based on plots of relatedness estimates)
A_low <- construct_adj_matrix(mle_CIs[FSVC_pair_ind,], Entry = 'r2.5.')
A_high <- construct_adj_matrix(mle_CIs[FSVC_pair_ind,], Entry = 'r97.5.') # Adj. matrix unfiltered using 'r97.5.' 
A_high[!(A_high >= (1-eps) & A_low > 0.5)] <- 0 # Edit s.t. only clonal have weight
A_high[(A_high >= (1-eps) & A_low > 0.5)] <- 1 # Edit s.t. only clonal have weight
G_high <- graph_from_adjacency_matrix(A_high, mode='upper', diag=F, weighted=T) # Construct graph 
V(G_high)$color <- cols_cities[metadata[attr(V(G_high), "names"), "City"]]
plot(G_high, vetex.color = V(G_high)$color, vertex.label = NA, vertex.size = 3)
cities_ <- unique(metadata[attr(V(G_high), "names"), "City"])
legend('topleft', pch = 16, bty = 'n', cex = 0.5, pt.cex = 1, 
       col = cols_cities[cities_], legend = cities_, inset = -0.1)
writeLines(sprintf("%s clonal component/s with two or more samples", 
                   sum(components(G_high)$csize > 1))) 
image(A_high)

#=============================================================
# Some exploration of whichever graph was last generated
#=============================================================
writeLines("Samples whose removal increases the number of graph components:") 
articulation_points(G_high)

# Cliques
Cliques <- max_cliques(G_high, min = 2) 
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

# Are the highly connected points those with less data
# Plot clique_count_per_sample against amonunt of data
# Clearly the sample with 150 memberships is an outlier: "Pf055" 
plot(y = as.numeric(clique_count_per_sample), 
     x = sample_snp_counts[names(clique_count_per_sample)], pch = 20, 
     ylab = "Number of clique memberships per sample", 
     xlab = "Available data per sample (SNPs)")

# Zoom in: seems to be a trend below 50 SNPs   
plot(y = as.numeric(clique_count_per_sample), 
     x = sample_snp_counts[names(clique_count_per_sample)], pch = 20, 
     ylab = "Number of clique memberships per sample", 
     xlab = "Available data per sample (SNPs)", 
     xlim = c(0,60))

# What happens if we remove those samples (20 samples)
problem_sample_ind <- sample_snp_counts[names(clique_count_per_sample)] < 50

# What happens if we remove the most overlapped? 
subG <- subgraph(G_high, names(problem_sample_ind)[!problem_sample_ind])
components(subG)
sum(components(subG)$csize > 1)
articulation.points(subG)
plot(subG, vertex.label = NA, vertex.size = 3)

# Community detection (as alternative to more strictly defined components and cliques)
wc <- cluster_walktrap(subG) 
plot(wc, subG, vertex.label = NA, vertex.size = 3, vetex.color = V(G_high)$color)

# # What happens if we remove the most overlapped? 
# ccps_threshold = 11
# subG <- subgraph(G_high, names(clique_count_per_sample[clique_count_per_sample < ccps_threshold]))
# plot(subG, vertex.label = NA, vertex.size = 3)
# writeLines(sprintf("Removes %s samples", sum(clique_count_per_sample >= ccps_threshold)))
# writeLines(sprintf("%s clonal component/s with two or more samples", 
#                    sum(components(subG)$csize > 1))) 
# # Sample SNP counts of removed samples
# removed_sample_snp_counts <- sample_snp_counts[names(which(clique_count_per_sample >= ccps_threshold))]
# barplot(removed_sample_snp_counts)
# writeLines("The metadata for the removed samples are:")
# metadata[names(removed_sample_snp_counts), ] 






