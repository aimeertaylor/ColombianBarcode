################################################################################
################################################################################
rm(list = ls())
library(igraph) 
library(RColorBrewer)
source('../igraph_functions.R') # construct_adj_matrix
eps <- 0.01 # threshold below which LCI is considered zero in Taylor et al. 2019
a = 0; b = 1; c = 0.5 # Scaling parameters for plotting edge width and transparancy 
PDF <- TRUE
label_art_pnts <- TRUE

# Load relatedness results, metadata and data
freqs_used <- "Taylor2020"
load(sprintf('../../RData/mles_CIs_extended_freqs%s_meta.RData', freqs_used)) 
load(sprintf('../../RData/metadata_extended.RData'))
load(file = "../../RData/snpdata_extended.RData")

# Load component and extract sids
LCI_threshold <- 0.75 # Toggle between 0, 0.75 and 0.95 
load(file = sprintf("../../RData/Clonal_components_extended_all_LCIthrehold_%s.RData", LCI_threshold))
sids <- Clonal_components$cc_7

# Sort sids, inc. those on either spurious edge (see Clonal_components.pdf)
break1 <- c("03014D0", "TC01D0", "SPT26314", "SPT26309","Pf067","Pf030","SPT26313")
break2 <- c("PW0067-C", "PW0080-C")
sids_nucleus <- sids[!sids %in% c(break1, break2)]
sids <- c(sids_nucleus, break1, break2) # Re-group ordered by breaks

# Inspect metadata
table(metadata[unlist(Clonal_components$cc_7), c("City", "Year")]) # all 
table(metadata[break1, c("City", "Year")]) # break one 
table(metadata[break2, c("City", "Year")]) # break two
table(metadata[sids_nucleus, c("City", "Year")]) # nucleus

# Plot data
image(as.matrix(snpdata[,sids])) # SNP data suggest definitely not clones

# Extract relatedness estimates (toggle between nucleus or not)
mles <- mle_CIs[mle_CIs$individual1 %in% sids_nucleus & mle_CIs$individual2 %in% sids_nucleus,]
#mles <- mle_CIs[mle_CIs$individual1 %in% sids & mle_CIs$individual2 %in% sids,]

# Inspect low snp_counts
table(sort(mles$snp_count[mles$snp_count <= 100]))

# Define city cols
cities <- unique(metadata$City)
cols_cities <- array(c(rev(brewer.pal(5, 'Spectral')), 
                       brewer.pal(length(cities)-5, 'Dark2')), 
                     dimnames = list(cities))

# Using definition of clone where UCI "touches" one and LCI > threshold
A_low <- construct_adj_matrix(mles, Entry = 'r2.5.')
A_high <- construct_adj_matrix(mles, Entry = 'r97.5.')
A_est <- construct_adj_matrix(mles, Entry = 'rhat')

A_clonal <- A_est 
A_clonal[A_high < (1-eps) | A_low <= LCI_threshold] <- 0
writeLines(sprintf("Minimum point estimate among those considered clones: %s", 
                   round(sort(unique(as.vector(A_clonal)))[2],2)))

# Create graph 
G_clonal <- graph_from_adjacency_matrix(A_clonal, mode='upper', diag=F, weighted = T) 
art_pnts <- names(articulation.points(G_clonal))
ccompnts <- components(G_clonal)

# Add metadata
V(G_clonal)$marker_count <- metadata[V(G_clonal)$name,  "snp_count"]
V(G_clonal)$city <- metadata[V(G_clonal)$name, "City"]
V(G_clonal)$color <- cols_cities[V(G_clonal)$city] 

# Highlight articulation points if any
if(length(art_pnts) > 0 & label_art_pnts){
  vertex_labels <- V(G_clonal)$name
  vertex_labels[!V(G_clonal)$name %in% art_pnts] <- ""
}

# Rescale edge weights for plotting
if (identical(round(unique(E(G_clonal)$weight)), 1)){ # If all one
  weights_rescaled <- E(G_clonal)$weight
} else { # Otherwise
  weights_rescaled = a + ((b-a) * (E(G_clonal)$weight - min(E(G_clonal)$weight)) / (max(E(G_clonal)$weight) - min(E(G_clonal)$weight)))
}

# Plot graph 
plot(G_clonal,
     vertex.label.color = "black", 
     vertex.label.cex = 1, 
     vertex.label = metadata[V(G_clonal)$name, "Year"],
     vertex.frame.color = NA, 
     vertex.shape = c("square","circle")[metadata[V(G_clonal)$name, "PloSGen2020"] + 1],
     edge.width = weights_rescaled * c, 
     edge.color = sapply(weights_rescaled, function(x) adjustcolor('black', alpha.f = x)),
     vertex.size = (1/c)*log(metadata[V(G_clonal)$name, "snp_count"]))

# Title
title(main = sprintf("Clonal def. UCI >= (1-%s) & LCI > %s", 
                     eps, LCI_threshold), cex.main = 1)

# Legend
cities_ <- unique(V(G_clonal)$city)
legend('bottomleft', pch = 16, bty = 'n', cex = 0.5, pt.cex = 1, 
       inset = -0.1, col = cols_cities[cities_], legend = cities_)

# Extract sids per component
# include components of size one (for Compare_components.R)
Clonal_components <- lapply(1:length(ccompnts$csize), function(cc){
  names(which(ccompnts$membership == as.numeric(cc)))
})
names(Clonal_components) <- paste0("cc_", 1:length(Clonal_components))

# How many non-cliques? 
# Check to see how many cliques there are per clonal_component
clique_count_per_component <- sapply(names(Clonal_components), function(cc){
  sids <- Clonal_components[[cc]]
  if(length(sids) > 1){
    subG <- subgraph(G_clonal, sids)
    x <- maximal.cliques.count(subG, min = 2)
    # Plot, remembering this is after setting edges to zero if not clonal
    if(x != 1) plot(subG, 
                    main = cc, 
                    sub = x,
                    vertex.color = cols_cities[metadata[V(subG)$name, "City"]] , 
                    edge.width = weights_rescaled * c, 
                    edge.color = "black", 
                    vertex.size = (1/c)*log(metadata[V(G_clonal)$name, "snp_count"]), 
                    vertex.label.size = 0.5, 
                    vertex.label = metadata[V(subG)$name, "Year"], #V(subG)$name, #
                    vertex.label.color = "black", 
                    vertex.shape = c("square","circle")[metadata[V(subG)$name, "PloSGen2020"] + 1])
    return(x)
  } else {
    return(0)
  }
})





